"""
Parse an arbitrary amount of fields in a VCF file converting it into a table for easier use by other.

2022-03-25: Matthew Wells
"""

from collections import namedtuple
from enum import Enum, auto
import os

vcffile = namedtuple('VCF', [])
# defining VCF line types globally, not putting everything in a data class as no need to increase memory usage
file_format = namedtuple("vcf_format", ["filetype"])
source = namedtuple("source", ["source"])
reference = namedtuple("reference", ["reference"])
file_date = namedtuple("file_date", ["date"])
phasing = namedtuple("phasing", ["phasing"])
command_line = namedtuple("command_line", ["command_line"])
filter_vcf = namedtuple("filter", ["id", "description"])
info = namedtuple("info", ["ID", "number", "type", "description", "source", "version"]) 
format_vcf = namedtuple("format_vcf", ["ID", "Number", "type", "description"]) # dynamic fieldd
alt = namedtuple("alternate", ["ID", "description"])
assembly = namedtuple("assembly", ["url"])
contig = namedtuple("contig", ["ID", "length", "url"]) # can be more then one contig present, fields need to be expanded on
sample = namedtuple("sample", []) # Dynamic field
pedigree = namedtuple("pedigree", []) # Dynamic field
pedigreedb = namedtuple("pedigree_db", ["url"])


vcf_info_format = namedtuple('VCFHeaderInfo', ['ID', 'Number', 'Type', 'Description',
'Source','Version'])

class VCFLines(Enum):
    FILEFORMAT = auto()
    SOURCE = auto()
    REFERENCE = auto()
    FILEDATE = auto()
    PHASING = auto()
    COMMANDLINE = auto()
    FILTER = auto()
    INFO = auto()
    FORMAT = auto()
    ALT = auto()
    ASSEMBLY = auto()
    CONTIG = auto()
    SAMPLE = auto()
    PEDIGREE = auto()
    PEDIGREEDB = auto()


def ParseVCFFields(fp: os.path) -> list:
    """Determine the fields kept within the vcf file
    
    VCF file do follow a standard which can be easily googled, however fields may change within them.
    Additionally, VCF files are often compressed so in the future compression libraries will be required

    :param fp: The vcf file to be read
    :return: Return a list of the column headers within the vcf file
    """
    headers_and_information = []

    with open(fp, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                headers_and_information.append(line)
            else:
                break


def process_vcf_line(vcf_line: str) -> str:
    """Return values format or info ID
    """

    if vcf_line.startswith('##fileformat'):
        split_vcf_line(vcf_line, VCFLines.FILEFORMAT)
    elif vcf_line.startswith('##source'):
        split_vcf_line(vcf_line, VCFLines.SOURCE)
    elif vcf_line.startswith('##reference'):
        split_vcf_line(vcf_line, VCFLines.REFERENCE)
    elif vcf_line.startswith('##filedate'):
        split_vcf_line(vcf_line, VCFLines.FILEDATE)
    elif vcf_line.startswith('##phasing'):
        split_vcf_line(vcf_line, VCFLines.PHASING)
    elif vcf_line.startswith('##commandline'):
        split_vcf_line(vcf_line, VCFLines.COMMANDLINE)
    elif vcf_line.startswith('##FILTER'):
        split_vcf_line(vcf_line, VCFLines.FILTER)
    elif vcf_line.startswith('##INFO'):
        split_vcf_line(vcf_line, VCFLines.INFO)
    elif vcf_line.startswith('##FORMAT'):
        split_vcf_line(vcf_line, VCFLines.FORMAT)
    elif vcf_line.startswith('##ALT'):
        split_vcf_line(vcf_line, VCFLines.ALT)
    elif vcf_line.startswith('##assembly'):
        split_vcf_line(vcf_line, VCFLines.ASSEMBLY)
    elif vcf_line.startswith('##contig'):
        split_vcf_line(vcf_line, VCFLines.CONTIG)
    elif vcf_line.startswith('##SAMPLE'):
        split_vcf_line(vcf_line, VCFLines.SAMPLE)
    elif vcf_line.startswith('##PEDIGREE'):
        split_vcf_line(vcf_line, VCFLines.PEDIGREE)
    elif vcf_line.startswith('##pedigreeDB'):
        split_vcf_line(vcf_line, VCFLines.PEDIGREEDB)

    return '' # remember to delete this!!!


def split_vcf_line(vcf_line: str, line_type: VCFLines) -> vcf_info_format:
    """Split the vcf line to to extract field info
    
    :return: The split vcf line with format fields specified
    """
    ...



