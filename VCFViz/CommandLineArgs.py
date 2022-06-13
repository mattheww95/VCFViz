"""
Parse the command line parameters and pass them to the program. Dispatch
the right args to the right function.

2022-05-29: Matthew Wells
"""
import argparse
import VCFViz.InputOptions as InputOptions
import os
from typing import Any
import sys

class Argparser:
    """
    handle the args for the different input options

    I organized some of the parser options in a dumb way but oh well...
    """
    args = []

    functions_call = {
                    "input-file": InputOptions.process_submission_sheet, 
                    "directory-glob": InputOptions.glob_directories, 
                    "wastewater-run": InputOptions.wastewater_run,
                    "summarize-excel": InputOptions.create_summary_excel_report
                    }
    def __call__(self) -> Any:
       
        parser = argparse.ArgumentParser(description="A VCFParser edit with minor changes.")
        subparsers = parser.add_subparsers(help="Pick a run mode for vcfparser.")
        #--- Submission Sheet Entry ---
        parser_1 = subparsers.add_parser("input-file", help="Run vcfparser with an input file.")
        parser_1.add_argument("-s", "--sample-sheet", help="Input file of samples names and paths to use")
        parser_1.add_argument("-o", "--output-directory", help="The output directory to use, default is current directory", 
        default=os.getcwd())
        parser_1.add_argument("-c", "--coverage-threshold", help="Set the minimum depth of coverage for allele prescence, default is 0", default=0, type=int)
        parser_1.add_argument("-m", "--metadata", help="The metadata sheet to use for subsetting VCF files.")

        #--- Directory Glob Entry ---
        parser_2 = subparsers.add_parser("directory-glob", help="Run vcfparser by passing in directories with a glob pattern.")
        parser_2.add_argument("-i", "--ivar-directory", help="The directory containing ivar outputs files.")
        parser_2.add_argument("-b", "--bam-directory", help="The directory containing the bamfiles matching the Ivar filies")
        parser_2.add_argument("-o", "--output-directory", help="The output directory to use, default is current directory", 
        default=os.getcwd())
        parser_2.add_argument("-c", "--coverage-threshold", help="Set the minimum depth of coverage for allele prescence, default is 0", default=0, type=int)
        parser_2.add_argument("-m", "--metadata", help="The metadata sheet to use for subsetting VCF files.")

        #--- Wastewater Directory Run ---
        parser_3 = subparsers.add_parser("wastewater-run", help="Run vcfparser on a reportable directory setup by the wastewater group.")
        parser_3.add_argument("-i", "--input-directory", help="Input of wastewater data configured directory")
        parser_3.add_argument("-c", "--coverage-threshold", help="Set the minimum depth of coverage for allele prescencem default is 30", default=30, type=int)
        parser_3.add_argument("-m", "--metadata", help="The metadata sheet to use for subsetting VCF files.")

        #--- Post run to optionally summarize html reports into a spreadsheet
        parser_4 = subparsers.add_parser("summarize-excel", help="Create a summary excel file of the final html data.")
        parser_4.add_argument("-d", "--directory-html", help="The directory containing the html run data.")
        parser_4.add_argument("-o", "--output-path", help="Place to output summary data as an excel file.")


        if len(self.args) == 0:
            parser.print_help()
            sys.exit(-1)
        else:
            parser_args = parser.parse_args()
            if len(self.args) == 1:
                parser.print_help()
                exit(-1)
            if len(self.args) == 2:
                parser.print_help()
                exit(-1)
            function_to_call = self.args[1]
            self.functions_call[function_to_call](**parser_args.__dict__)
    
    def __init__(self, *args, **kwargs):
        self.args = args[0]

def run_main():
    """
    The main entry point for the program.
    """
    run1 = Argparser(sys.argv)
    run1()

if __name__ == "__main__":
    run_main(sys.argv)
    #test_sub_sheet = "tests/test_vcfparser_subsheet_.txt"
    #test_metadata_sheet = "tests/VCFParser_tester.txt"
    #cov_thresh = 30
    #outdir = "./tests"
    #sys.argv.extend(["input-file", "-s", test_sub_sheet, "-o", outdir, "-c", str(cov_thresh), "-m", test_metadata_sheet])
    #t1 = Argparser(sys.argv)
    #t1()