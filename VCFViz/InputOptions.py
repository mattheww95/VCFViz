"""
Handle different input options being either globs, input sheet or specification of two files

2022-05-30: Matthew Wells
"""



try:
    from VCFViz.VCFlogging import VCFLogger as vlog
    from VCFViz.VCFToJson import ReadIvar, ReadVCF, VariantData
    from VCFViz import RenderHTML
    from VCFViz.RenderHTML import VCFDataHTML
    from VCFViz import CoverageData
    from VCFViz.CreateExcelReports import HTMLToExcel
except ImportError:
    from VCFlogging import VCFLogger as vlog
    from VCFToJson import ReadIvar, ReadVCF, VariantData
    import RenderHTML
    from RenderHTML import VCFDataHTML
    import CoverageData
    from CreateExcelReports import HTMLToExcel

    
from datetime import datetime
import os
import glob

def determine_vcf_format(sample_file):
    """
    Determine which vcf reader method to use
    """
    if sample_file[sample_file.rindex(".")+1:] == "tsv":
        return (ReadIvar, "*.tsv")
    else:
        return (ReadVCF, "*.vcf")

#Submission sheet input (Retain sample order)
def process_submission_sheet(sample_sheet: str, metadata: str, coverage_threshold: int, output_directory: str):
    """
    Process a submission sheet that provides:
        - sample name
        - Ivar sheet path
        - bam path
    """
    samples = []
    var_data = []
    sheet_cov_data = []
    file_reader_method = None
    glob_string = None
    with open(sample_sheet, 'r') as samples_:
        for i in samples_.readlines():
            val = i.strip().split("\t")
            if len(val) != 3:
                vlog.logger.critical(val)
                vlog.logger.critical("Specified sheet does not match needed criteria.")
                vlog.logger.critical("Sheet should be tab delimited and ordered: sample name, vcf path, bam path")
                exit(-1)
            file_reader_method, glob_string = determine_vcf_format(val[1])
            ivar_data = file_reader_method(val[1])
            ivar_data.sample_name = val[0]
            var_data.append(ivar_data)
            cov_data = CoverageData.SampleMap(val[0], val[2])
            cov_data.sample_name = val[0]
            sheet_cov_data.append(cov_data)
            samples.append(val[0])
        
            #samples_process[val[0]] = (ivar_data, cov_data) # 1: ivardata 2: bam path
    vlog.logger.info(f"Creating coverage data for {len(samples)} samples and outputting data to {output_directory}.")
    sample_cov_data = CoverageData.create_sample_coverages(samples, output_directory, sheet_cov_data)
    rendered = RenderHTML.VCFDataHTML(VariantData(var_data), metadata, None, coverage_threshold, output_directory, sample_cov_data)
    rendered.combine_html_plots()

#Glob directories
def glob_directories(ivar_directory:str, bam_directory: str, metadata: str, coverage_threshold: int, output_directory: str):
    """
    The main function to call in prepareing the samples
    """
    file_reader_method = None
    glob_ext = None
    #TODO this is a bottle neck that can be easily eliminated in the future
    if file_reader_method or glob_ext is None:
        file_reader_method, glob_ext = determine_vcf_format(os.path.join(ivar_directory, os.listdir(ivar_directory)[0]))
    
    #ivar_data = [file_reader_method(os.path.join(ivar_directory, i)) for i in os.listdir(ivar_directory) if os.path.splitext(i)[-1].lower() == glob_ext]
    ivar_data = [file_reader_method(i) for i in glob.glob(os.path.join(ivar_directory, glob_ext))]
    if not ivar_data:
        vlog.logger.critical(f"Could not find vcf files in {os.path.join(ivar_directory, glob_ext)}")
    vcf_html = VCFDataHTML(VariantData(ivar_data), metadata, bam_directory, coverage_threshold, output_directory)
    vcf_html.combine_html_plots()

#cmd line sample specification
def wastewater_run(input_directory, metadata, coverage_threshold):
    """
    Run the new vcfparser on the wastewater directories
    """
    start = datetime.now()
    for i in os.listdir(input_directory):
        variants = os.path.join(input_directory, i, "variants")
        bams = os.path.join(input_directory, i, "bam")
        if os.path.isdir(variants) and os.path.isdir(bams):
            out_dir = os.path.join(input_directory, i)
            try:
                glob_directories(variants, bams, metadata, coverage_threshold, out_dir)
            except RuntimeError:
                pass
        else:
            pass
    end = datetime.now()
    vlog.logger.info(f"Finished: {input_directory} in {end - start}")

def create_summary_excel_report(directory_html, output_path):
    """
    Create a summary report of the html run information output into excel
    """
    vlog.logger.info("Creating Excel summary of HTML information.")
    xx = HTMLToExcel(directory_html, output_path)
    xx.html_dict_to_excel()
    vlog.logger.info(f"Excel file output to {output_path}")


if __name__ == "__main__":
    test_sub_sheet = "tests/test_vcfparser_subsheet_.txt"
    test_metadata_sheet = "tests/VCFParser_tester.txt"
    cov_thresh = 30
    outdir = "./tests"
    cov = 30
    parser_sheet = "tests/VCFParser_tester.txt"