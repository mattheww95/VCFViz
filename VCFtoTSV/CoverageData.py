"""
Recreate the '.cache_snv_coverages.json' file output by vcfparser for 
incorporating coverage in the rendered images. The previous version used
pysam, however pysam is quite slow therefore this step will use samtools in
a subprocess instead. 
The json snv file consisted of a dictionary containing the absolute path to
the bam then a dictionary of the values and theyre assigned depths.

2022-05-26: Matthew Wells
"""

import datetime
import os
import glob
from re import S
import subprocess
from typing import List
from VCFlogging import VCFLogger as vlog
import json



class SampleMap:
    """
    Per an individual sample, create some coverage information based
    on a sample name provided. So the class will need to: 
        - recurse a directory known to contain the bam file
        - match the file based on the sample name
        - check if it has an index
            - if not create one
        - run samtools depth, and create a json file of the information
    """

    def __init__(self, sample_name, directory_bam) -> None:
        self.sample_name = sample_name
        self.directory_bam = directory_bam
        self.bam_abs_path = None
        self.bai_abs_path = None
        if os.path.isfile(directory_bam):
            self.initialize_sample()
        else:
            self.find_sample()

    def initialize_sample(self):
        """
        Pair up a sample without a glob path
        """
        sample_name = os.path.basename(self.directory_bam)
        sample_name = sample_name[:sample_name.index(".")]
        bai_path = self.directory_bam + ".bai"
        if os.path.isfile(bai_path):
            self.bai_abs_path = bai_path
            self.bam_abs_path = self.directory_bam
        else:
            vlog.logger.info(f"Creating index for sample {self.sample_name}")
            # calling directory bam as should be single file in this instance
            new_idx = subprocess.Popen(["samtools", "index", "-b", self.directory_bam])
            status = new_idx.wait()
            out, err = new_idx.communicate()
            if status == 1:
                vlog.logger.critical(f"Broken SAM/BAM file: {self.bam_abs_path}")
                vlog.logger.critical(f"Subprocess stdout: {out}")
                vlog.logger.critical(f"Program stderr: {err}")
                raise ValueError(self.sample_name)
            vlog.logger.debug(f"Created index for {self.sample_name}")
            self.bai_abs_path = self.directory_bam + ".bai"
            self.bam_abs_path = self.directory_bam
        return 0

    def find_sample(self):
        """
        search the specified directory for the samples name and its bam
        """
        
        bam_dir = self.directory_bam
        for file in glob.glob(f"{bam_dir}/*.bam"): #TODO clean up string formatting in glob
            i = os.path.basename(file)
            sample_name = i[:i.index(".")] # get name up to file name
            if sample_name == self.sample_name:
                self.bam_abs_path = os.path.abspath(file)

                # samtools creates an index also with the filename and bai tacked on
                bai_path = self.bam_abs_path + ".bai" 
                if not os.path.isfile(bai_path):
                    vlog.logger.info(f"Creating index for sample {self.sample_name}")
                    new_idx = subprocess.Popen(["samtools", "index", "-b", self.bam_abs_path])
                    status = new_idx.wait()
                    out, err = new_idx.communicate()
                    if status == 1:
                        vlog.logger.critical(f"Broken SAM/BAM file: {self.bam_abs_path}")
                        vlog.logger.critical(f"Subprocess stdout: {out}")
                        vlog.logger.critical(f"Program stderr: {err}")
                        raise ValueError(self.sample_name)

                    vlog.logger.debug(f"Created index for {self.sample_name}")
                self.bai_abs_path =  bai_path
                vlog.logger.debug(f"Appending bam index for {self.sample_name}")
                return 0
        vlog.logger.critical(f"Could not find bam file for sample {self.sample_name}")
        raise ValueError(f"Could not find bamfile for sample: {self.sample_name}")

class SamplesCoverage:
    """
    take a list of SampleMap's and calculate their depths.

    All SampleMaps should have a map and index, as SampleMap's
    initializer will throw an error if any sample fails
    """
    samples_coverage = {} # keep all of the coverage data static
    def __init__(self, samples: List[SampleMap]) -> None:
        self.samples = samples
        for sample in self.samples: # initialize samples in dictionary
            self.samples_coverage[sample.sample_name] = {}
        #self.retrieve_coverage() # move this out of init
    
    def retrieve_coverage(self):
        """
        call the samtools depth process on the list of samples
        """
        
        chunks_bam = self.chunk_list(5, self.samples)
        for chunk in chunks_bam:
            self.call_coverage_program(chunk)

    def call_coverage_program(self, samples_list):
        """
        call samtools on a sample list to retrieve coverage and return a dictionary of data
        """
        samtools_call = ["samtools", "depth", "-aa"]
        samtools_call.extend([i.bam_abs_path for i in samples_list])
        vlog.logger.info("Preparing sample depth information, this may be slow")
        vlog.logger.info(f"Samples being processed {[i.sample_name for i in samples_list]}")
        start = datetime.datetime.now()
        depth_task = subprocess.Popen(samtools_call, stdout=subprocess.PIPE)
        
        for i in depth_task.stdout:
            # Passing multiple files to samtools depth returns them in order and as bytes
            # so it is nessecary to encode the bytes to text for processing
            depth_data = i.decode("utf-8", "ignore").strip().split("\t")
            for k, cov in enumerate(depth_data[2:]): # skipping chormosome and postion in output
                # depth data is 1st pos and keeping as string to match original format
                self.samples_coverage[samples_list[k].sample_name][depth_data[1]] = int(cov)
        depth_task.communicate() #testing adding communicate as process kept faileing on the cluster remove if not nesseccary
        status = depth_task.poll()
        end = datetime.datetime.now()
        if status == 0:
            vlog.logger.info(f"Gathered coverage data. Process finished in {end - start} seconds")
        else:
            vlog.logger.critical("Could not compute coverage bam paths provided are as follows:")
            for i in samples_list:
                print("-", i.bam_abs_path)
            vlog.logger.critical(f"Samtools command: {' '.join(samtools_call)}")
            raise RuntimeError("Could not compute coverage, received samtools error.")

    @staticmethod
    def chunk_list(chunk_size, list_chunk):
        """
        Convert a list into chunks for batched processing.
        """
        chunks = []
        for i in range(0, len(list_chunk), chunk_size):
            chunk = list_chunk[i:i+chunk_size]
            chunks.append(chunk)
        return chunks


def create_sample_coverages(samples: List[str], search_dir: str, sample_maps: List[SampleMap] = None):
    """
    From all of the sample sheets specified create the sample map objects
    and pass it off to samples coverage to return the coverage obj
    :param samples: A list of sample_names
    :param search_dir: the directory containing bams
    """
    cache_path = os.path.join(search_dir, ".cache_snv_coverages.json")
    vlog.logger.info(f"Searching {cache_path} for depth cache.")
    if sample_maps is None:
        sample_maps = []
        for i in samples:
            sample_maps.append(SampleMap(i, search_dir))
    
    cov_data = SamplesCoverage(sample_maps)
    if os.path.isfile(cache_path):
        with open(cache_path, 'r') as cov_data_:
            try:
                stage_cov = json.load(cov_data_)
            except json.decoder.JSONDecodeError:
                vlog.logger.warning(f"Could not read SNP cache creating new one")
                cov_data.retrieve_coverage()
                write_cache(cache_path=cache_path, cov_data=cov_data)
                return cov_data

            # check that all samples are in cov info
            samples_to_recall = []
            for samp in sample_maps:
                if stage_cov.get(samp.sample_name) is None:
                    vlog.logger.warning(f"Missing coverage for sample {samp.sample_name}"\
                        f" in data cache, regenerating coverage information for missing sample.")
                    samples_to_recall.append(samp)
                    #cov_data.retrieve_coverage()

            if len(samples_to_recall) != 0:
                chunks = cov_data.chunk_list(5, samples_to_recall)
                for chunk in chunks:
                    vlog.logger.info(f"Updating coverage cache to include samples {[i.sample_name for i in chunk]}")
                    cov_data.call_coverage_program(chunk)
                # load and stage data
                # this needs to be refactored to not be copying this code later
                write_cache(cache_path=cache_path, cov_data=cov_data, mode='a')
                return cov_data

            vlog.logger.info(f"Reusing coverage data from previous program run.")
            for key in stage_cov.keys():
                for vals in stage_cov[key].keys():
                    stage_cov[key][vals] = int(stage_cov[key][vals]) 
            cov_data.samples_coverage = stage_cov
            return cov_data
    else:
        cov_data.retrieve_coverage()
    write_cache(cache_path=cache_path, cov_data=cov_data)
    return cov_data

def write_cache(cache_path, cov_data, mode = 'w'):
    """
    A method to write out the snv cache
    """
    if mode !='a':
        with open(cache_path, mode) as cov_file:
            vlog.logger.info(f"Creating depth data cache in {cache_path}.")
            json.dump(cov_data.samples_coverage, cov_file, indent=2)
    else:
        with open(cache_path, 'r') as cov_file:
            vlog.logger.info("Appending data to existing cache.")
            json_data = json.load(cov_file)
            cov_data.samples_coverage = {**json_data, **cov_data.samples_coverage} # merge the dictionaries
            write_cache(cache_path=cache_path, cov_data=cov_data)
    


if __name__=="__main__":
    test_list = ["1" for _ in range(3)]
    test_list2 = ["2" for _ in range(5)]
    test_list3 = ["3" for _ in range(7)]
    SamplesCoverage.chunk_list(5, test_list)
    SamplesCoverage.chunk_list(5, test_list2)
    SamplesCoverage.chunk_list(5, test_list3)

    t1 = SampleMap("22_WPG17_NE_0426_2", "./tests")
    #t1.find_sample()

    t2 = SampleMap("22_AB16_GP_0414", "./tests")
    #t2.find_sample()

    cc = SamplesCoverage([t1, t2])
    #cc.retrieve_coverage()
    create_sample_coverages(["22_WPG17_NE_0426_2", "22_AB16_GP_0414"], "./tests")

    #test should fail
    #t3 = SampleCoverage("22_AB16_tt_0414", "./tests")
    #t3.find_sample()

    #t3 = SampleMap("22_AB16_GG_0414", "./tests")
    #t3.find_sample()