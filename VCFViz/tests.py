"""
Create some unit tests for this program cause hell yeah

2022-05-19: Matthew Wells
"""


import VCFViz.CommandLineArgs as CommandLineArgs
import unittest
from VCFViz.RenderHTML import VCFDataHTML
from VCFViz.VCFToJson import ReadIvar
from VCFViz.VCFToJson import ReadVCF
from VCFViz.CoverageData import CoverageData
from VCFViz.VCFlogging import VCFLogger as vlog
import logging
import time
import InputOptions
import sys
import copy

vlog.logger.setLevel(logging.CRITICAL)

class TestVCFMethods(unittest.TestCase):    
    def test_initialize_voc_tables(self):
        ivar_data = ReadIvar("tests/22_AB16_GP_0414.tsv")
        ivar_data2 = ReadIvar("tests/22_WPG17_NE_0426_2.tsv")
        ivar_data_list = [ivar_data2, ivar_data]
        vcf_html = VCFDataHTML(ivar_data_list, "tests/VCFParser_20220516.txt", "./tests", 30, "/tmp")
        for key in vcf_html.figure_data.keys():
            for mut in vcf_html.figure_data[key].keys():
                vals = len(vcf_html.figure_data[key][mut])
                self.assertEqual(vals, len(ivar_data_list))

class TestVCFToJson(unittest.TestCase):

    def test_ReadVCF(self):
        vcf_file = ReadVCF("tests/test_vcf.vcf")
        for key in vcf_file.vcf_file.keys():
            self.assertEqual(vcf_file.vcf_file[key].row.POS, key)

    def test_ReadIvar(self):
        ivar_file = ReadIvar("tests/22_WPG17_NE_0426_2.tsv")
        for pos in ivar_file.vcf_info:
            for val in ivar_file.vcf_info[pos]:
                self.assertEqual(pos, val.POS)

class TestVCFRenderHTML(unittest.TestCase):
    """
    Updated functionality of the class broke the test, need to rewrite
    """
    ...
    #def test_MultipleFileReads(self):
    #    print("\n") # make logger output
    #    vlog.logger.critical("Running test with 10 files")
    #    start = time.time()
    #    files = [ReadIvar("tests/22_WPG17_NE_0426_2.tsv") for _ in range(10)]
    #    VCFDataHTML(files, "tests/VCFParser_20220516.txt")
    #    end =  time.time()
    #    vlog.logger.critical(f"Time to run 10 files {round(end - start, 3)} seconds")
    #    
    #    vlog.logger.critical("Running test with 100 files")
    #    start = time.time()
    #    files = [ReadIvar("tests/22_WPG17_NE_0426_2.tsv") for _ in range(100)]
    #    VCFDataHTML(files, "tests/VCFParser_20220516.txt")
    #    end =  time.time()
    #    vlog.logger.critical(f"Time to run 100 files {round(end - start, 3)} seconds")

class TestSampleMap(unittest.TestCase):
    """
    Test the sample Coverage examples
    """
    def test_find_sample(self):
        
        #t3 = CoverageData.SamplesCoverage([t1])
        #self.assertRaises(ValueError, t1.find_sample())
        with self.assertRaises(ValueError):
            CoverageData.SampleMap("22_AB16_tt_0414", "./tests")


class TestSamplesCoverage(unittest.TestCase):
    """
    Test Samplescoverage Fucntions
    """
    ...


class TestInputOptions(unittest.TestCase):
    """
    Unit tests for various submission types
    """

    def test_process_submission_sheet(self):
        test_sub_sheet = "tests/test_vcfparser_subsheet_.txt"
        test_metadata_sheet = "tests/VCFParser_tester.txt"
        cov_thresh = 30
        outdir = "/tmp"
        InputOptions.process_submission_sheet(test_sub_sheet, test_metadata_sheet, cov_thresh, outdir)

class TestCommandLineArgs(unittest.TestCase):
    """
    Tests for the various command line args, they should return type errors
    """
    def test_Argparser_input_file(self):
        args_t1 = copy.deepcopy(sys.argv)
        args_t1.extend("input-file")
        t1 = CommandLineArgs.Argparser(args_t1)
        with self.assertRaises(TypeError):
            t1()
    
    def test_Argparser_input_file(self):
        args_t2 = copy.deepcopy(sys.argv)
        args_t2.extend("directory-glob")
        t2 = CommandLineArgs.Argparser(args_t2)
        with self.assertRaises(TypeError):
            t2()

    def test_Argparser_input_file(self):
        args_t3 = copy.deepcopy(sys.argv)
        args_t3.extend("wastewater-run")
        t3 = CommandLineArgs.Argparser(args_t3)
        with self.assertRaises(TypeError):
            t3()



if __name__ == '__main__':
    unittest.main(verbosity=2)