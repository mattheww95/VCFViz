"""
Read in a VCF File and convert it to a json file for later formatting.

Proper VCF data lends itself to more of a json format after parsing
hence while dataclasses and tuples are nice it likely make moresense to
convert it to a dictionary.

Mandatory fields are howerver:
    1) #CHROM
    2) POS
    3) ID
    4) REF
    5) ALT
    5) QUAL
    7) FILTER
    8) INFO

They must be in that order as well.

2022-05-09: Matthew Wells
"""

from collections import namedtuple
from dataclasses import dataclass
from typing import NamedTuple
import os


class VCFTags:
    """
    A class to contain the format and info tags
    """
    contig: dict = {}
    FORMAT: dict = {}
    INFO: dict = {}

class FlagDescriptors(NamedTuple):
    ID: str
    NUMBER: str
    TYPE: str # specifies what number is as not always a number
    DESCRIPTION: str

class ReadVCF:
    """
    Read in the vcf file and parse the header.
    """
    __slots__ = ['file_name', '__dict__']
    header_info = VCFTags
    SPLIT_CHAR = "$" # Trying to find something ascii yet not used in informatics...
    TABLE_DELIMITER = "\t"
    VCFRow = namedtuple("VCFRow", ["row", "INFO", "FORMAT"])
    vcf_file = {}

    def __init__(self, file_name) -> None:
        self.table_start = None
        self.file_name = file_name
        self.read_vcf_header()  
        self.read_vcffile()  
    

    def read_vcffile(self):
        """
        Read the VCF File in starting at the point defined in the header. this method will
        have to take the vcf table column heads and prepare dictionaries from them of the 
        multiple different pieces of information.
        """
        with open(self.file_name, 'r') as vcf:
            vcf_data = vcf.readlines()
            iter_data = vcf_data[self.table_start:]
        cols = iter_data[0]
        vcf_rows = []
        col_vals = cols.strip().strip("#").split("\t")
        VCFData = namedtuple("VCFRow", col_vals)
        sample_start = col_vals.index("FORMAT") + 1 # format is last tag in standard vcf file, adding one as 0 indexed in cols
        for line in iter_data[1:]:
            vcf_row = VCFData(*line.strip().split("\t"))
            """
            Withing the vcf rows we get a format info field and a format field
            specific to the sample. These correspond to information in the fromat or infor fields
            FORMAT gives the order each samples ifnormation will show up
            INFO is sample specific
            """
            vcf_info_row = self.split_vcf_info_field(vcf_row.INFO)
            form_tags = vcf_row.FORMAT.split(":")
            sample_format_info = {}
            for sample in vcf_row[sample_start:]:
                sample_format_info[sample] = dict(zip(form_tags, sample.split(":")))
            mut_col = vcf_row.POS
            self.vcf_file[mut_col] = self.VCFRow(vcf_row, vcf_info_row, sample_format_info)

            
    def split_vcf_info_field(self, vcf_info_field):
        """
        A method to parse out the vcf info field and so that each tag can be matched to its describing
        information
        :param vcf_info_field: A information string from a vcf field
        """
        line_split = vcf_info_field.split(";")
        information_tags = {}
        for val in line_split:
            tags = val.split("=")
            information_tags[tags[0]] = tags[1]
        return information_tags
        
    def read_vcf_header(self):
        """
        Open the first vcf file and parse the header attributes, all header
        attributes begin with a ##
        """
        headers = []
        with open(self.file_name, 'r') as vcf:
            i = next(vcf)
            track_table_op = 0
            while i[:2] == "##":
                headers.append(i.strip().strip("##"))
                track_table_op += 1
                i = next(vcf)

        self.table_start = track_table_op # set where the information starts for the second read
        for line in headers:
            split_line = line.strip(">").split("<")
            split_line_len = len(split_line)
            if split_line_len > 1:
                type_line = split_line[0].strip("=")
                
                # three commas seperating the values, need to split on 
                # comma but the description tag
                # uses those within in it, so split_char replacement makes it 
                # easy to split the string in the right place
                info_line = split_line[1].replace(",", self.SPLIT_CHAR, 3).split(self.SPLIT_CHAR) 
                if type_line == "contig":
                    #TODO make sure this can be refactored for multiple contigs
                    type_line, self.header_info.__dict__["__annotations__"][type_line] = info_line
                elif split_line_len > 1:
                    type_line = split_line[0].strip("=")
                    vals_format = {}
                    for val in info_line:
                        strlin = val.split("=")
                        vals_format[strlin[0].strip('=').upper()] = strlin[1]
                    getattr(self.header_info, type_line.upper())[vals_format["ID"]] = FlagDescriptors(**vals_format)
                    
class IvarFields(NamedTuple):
        REGION: str
        POS: str
        REF: str
        ALT: str
        REF_DP: str
        REF_RV: str
        REF_QUAL: str
        ALT_DP: str
        ALT_RV: str
        ALT_QUAL: str
        ALT_FREQ: str
        TOTLA_DP: str
        PVAL: str
        PASS: str
        GFF_FEATURE: str
        REF_CODON: str
        REF_AA: str
        ALT_CODON: str
        ALT_AA: str



class ReadIvar:
    """
    Depending on the variant caller used it may not actually be a vcf and is infact a tsv of few fields...

    Ivar does weird things with bed intervals so correcting indexing will need to be performed. However,
    Ivar has a defined set of fields that do not change! and can be specified here.

    For future consideration: https://andersen-lab.github.io/ivar/html/manualpage.html
    Refer to the above manual page.
    """

    def __init__(self, filename):
        self.vcf_info = {} # declaring this with the class creates a shared attribute...
        self.filename = filename
        self.sample_name = self.get_sample_name(self.filename)
        self.read_ivar_file()
    
    @staticmethod
    def get_sample_name(file_name):
        file_only = os.path.basename(file_name)
        file_only = file_only[:file_only.index(".")] # drop extension
        return file_only

    def print_ivar_info(self):
        """
        For dev purposes, just to verify the ivar file is read properly
        """
        for pos in self.vcf_info:
            print(pos, self.vcf_info[pos])
    
    def read_ivar_file(self):
        """
        Read the ivar tsv and return its initialized vcf.
        """
        with open(self.filename, 'r') as vcf:
            vcf_file = vcf.readlines()[1:] # skip columns as specified
            for line in vcf_file:
                lines_split = line.strip().split("\t")
                #ivar_row = self.IvarFields(*lines_split)
                ivar_row = IvarFields(*lines_split)
                # Checking if ivar position is empty as I beleive ivar provided multiple
                # entries for the mutations at the same mutation
                if self.vcf_info.get(ivar_row.POS) is None:
                    self.vcf_info[ivar_row.POS] = []
                    self.vcf_info[ivar_row.POS].append(ivar_row)
                else:
                    self.vcf_info[ivar_row.POS].append(ivar_row)

if __name__=="__main__":
    vcf_file = ReadVCF("tests/test_vcf.vcf")
    #for key in vcf_file.vcf_file.keys():
    #    print(key, vcf_file.vcf_file[key])
    ivar_file = ReadIvar("tests/22_WPG17_NE_0426_2.tsv")
    ivar_file.print_ivar_info()