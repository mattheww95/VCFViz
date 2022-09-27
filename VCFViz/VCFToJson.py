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

from collections import namedtuple, defaultdict
from dataclasses import dataclass
from VCFlogging import VCFLogger
from typing import NamedTuple, List
from enum import Enum, auto
import sys
import os
import glob


class MutationTypes(Enum):
    INS = auto()
    DEL = auto()
    SNP = auto()
    MNP = auto()

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
    #header_info = VCFTags
    #SPLIT_CHAR = "$" # Trying to find something ascii yet not used in informatics...
    #TABLE_DELIMITER = "\t"
    #VCFRow = namedtuple("VCFRow", ["row", "INFO", "FORMAT"])
    #vcf_file = defaultdict(lambda: list())

    def __init__(self, file_name) -> None:
        self.header_info = VCFTags
        self.SPLIT_CHAR = "$" # Trying to find something ascii yet not used in informatics...
        self.TABLE_DELIMITER = "\t"
        self.VCFRow = namedtuple("VCFRow", ["row", "INFO", "FORMAT"])
        self.vcf_file = defaultdict(lambda: list())
        self.table_start = None
        self.file_name = file_name
        self.sample_name = ReadIvar.get_sample_name(file_name)
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
        cols = iter_data[0].strip(" #\n").replace(".", "").replace("/", "")
        vcf_rows = []
        col_vals = cols.split("\t")
        
        # Make all sample names valid python variables, this approach is hacky however
        col_vals = [i if not i[0].isdigit() else "SAMPLE_" + i for i in col_vals]

        VCFData = namedtuple("VCFRow", col_vals)
        sample_start = col_vals.index("FORMAT") + 1 # format is last tag in standard vcf file, adding one as 0 indexed in cols
        for line in iter_data[1:]:
            vcf_row = VCFData(*line.strip().split("\t"))
            """
            Withing the vcf rows we get a format info field and a format field
            specific to the sample. These correspond to information in the fromat or infor fields
            FORMAT gives the order each samples information will show up
            INFO is sample specific
            """
            vcf_info_row = self.split_vcf_info_field(vcf_row.INFO)
            form_tags = vcf_row.FORMAT.split(":")
            sample_format_info = {}
            for data, sample in zip(vcf_row[sample_start:], vcf_row._fields[sample_start:]):
                sample_format_info[sample] = dict(zip(form_tags, data.split(":")))
            mut_col = vcf_row.POS
            # Keep position indexing keep multiple rows per entry
            self.vcf_file[mut_col].append(self.VCFRow(vcf_row, vcf_info_row, sample_format_info))

            
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
        REGION: str = None
        POS: str = None
        REF: str = None
        ALT: str = None
        REF_DP: str = None
        REF_RV: str = None
        REF_QUAL: str = None
        ALT_DP: str = None
        ALT_RV: str = None
        ALT_QUAL: str = None
        ALT_FREQ: str = None
        TOTAL_DP: str = None
        PVAL: str = None
        PASS: str = None
        GFF_FEATURE: str = None
        REF_CODON: str = None
        REF_AA: str = None
        ALT_CODON: str = None
        ALT_AA: str = None
        TYPE: str = None



class ReadIvar:
    """
    Depending on the variant caller used it may not actually be a vcf and is infact a tsv of few fields...

    Ivar does weird things with bed intervals so correcting indexing will need to be performed. However,
    Ivar has a defined set of fields that do not change! and can be specified here.

    For future consideration: https://andersen-lab.github.io/ivar/html/manualpage.html
    Refer to the above manual page.

    TODO adding creation of contiguity checking
    TODO Add type field to Ivar info
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

    def determine_ivar_type(self, _alt):
        """
        Determine the type of change ivar is reporting be it, insertion, deletion or snp.
        Ivar does not record MNP changes sadly
        """
        if "+" in _alt:
            return MutationTypes.INS.name
        elif "-" in _alt:
            return MutationTypes.DEL.name
        else:
            return MutationTypes.SNP.name
    
    def read_ivar_file(self):
        """
        Read the ivar tsv and return its initialized vcf.
        """
        with open(self.filename, 'r') as vcf:
            vcf_file = vcf.readlines()[1:] # skip column names
            for line in vcf_file:
                lines_split = line.strip().split("\t")
                """
                Adding type field to entries based on characters in Ivar Field, 
                fields 2 and 3 (zero based) correspond to the needed fields
                """
                lines_split.append(self.determine_ivar_type(lines_split[3]))
                ivar_row = IvarFields(*lines_split)
                # Checking if ivar position is empty as I beleive ivar provided multiple
                # entries for the mutations at the same mutation
                if self.vcf_info.get(ivar_row.POS) is None:
                    self.vcf_info[ivar_row.POS] = list()
                    self.vcf_info[ivar_row.POS].append(ivar_row)
                else:
                    self.vcf_info[ivar_row.POS].append(ivar_row)
    
    def create_contgous_entries(self):
        """
        order contigous mutations so that they can be combined and new ivar entries can be created
        """
        positions = [int(i) for i in self.vcf_info] # get all keys
        positions.sort() # sorts in place
        pos_iter = iter(positions)
        contigous_groups = list()
        contigous_positions = list()
        try:
            contigous_positions.append(next(pos_iter))
        except StopIteration:
            sys.stderr.write("No mutations found in sample.\n")
            exit(-1)
        while True:
            try:
                ps = next(pos_iter)
            except StopIteration:
                break
            else:
                if (contigous_positions[-1] + 1) == ps:
                    contigous_positions.append(ps)
                else:
                    if len(contigous_positions) > 1:
                        contigous_groups.append(contigous_positions)
                    contigous_positions = list()
                    contigous_positions.append(ps)
        self.create_joined_groups(contigous_groups)
    
    def create_joined_groups(self, contig_groups: List[list]):
        """
        Create joined VCFEntries for the MNP's created 
        One issue to keep tabs on is that some ivar rows contain multiple entries per a site
        """
        for group in contig_groups:
            mnp_groups = list()
            for val in group:
                if len(self.vcf_info[str(val)]) >= 2:
                    mnp_groups.append(self.vcf_info[str(val)])
                else:
                    mnp_groups.extend(self.vcf_info[str(val)])
    

class CommonVCF(NamedTuple):
    SampleName: str
    CHROM: str
    TYPE: str
    TOTAL_DEPTH: int
    ALT_DEPTH: int
    POS: int
    ALT: str



class VariantData:
    """
    Class to take either Ivar data or VCF data and create a shared format to pass
    to the html renderer
    2022-09-21
    """
    variant_data = dict()
    def __init__(self, variants_list) -> None:
        self.variants_list = variants_list
        self.translate_vcf()
    
    def translate_vcf(self):
        """
        Pick the right translator to convert the VCF data into a common format for the heatmaps
        """
        #print(self.variants_list.count(type(self.variants_list[0])) == len(self.variants_list))
        out = [type(i) for i in self.variants_list]
        value_type = out[1:] == out[:-1]
        if not value_type:
            VCFLogger.logger.critical("A mixture of VCF files and Ivar files has been recorded. Bailing out")
            exit(-1)

        if out[0] == ReadIvar:
            VCFLogger.logger.info("You are using Ivar files.")
            self.ivar_translation()
        elif out[0] == ReadVCF:
            VCFLogger.logger.info("You are using VCF files.")
            self.vcf_translation()
        else:
            VCFLogger.logger.critical(f"{type(out[0])} is not yet supported for visualization")
            exit(-1)

    def ivar_translation(self):
        """
        Translate the Ivar output into the CommonVCF type format
        """
        for vcf in self.variants_list:
            for pos, values in vcf.vcf_info.items():
                for v in values:
                    #print(v)
                    alt_depth = int(v.ALT_DP) # get reverse read and forward reads
                    total_depth = int(v.TOTAL_DP)
                    mut_type = v.TYPE
                    pos = int(v.POS)
                    chrom = v.REGION
                    sample = vcf.sample_name
                    if self.variant_data.get(sample) is None:
                            self.variant_data[sample] = dict()
                    if self.variant_data[sample].get(pos) is None:
                        self.variant_data[sample][pos] = list()
                    self.variant_data[sample][pos].append(CommonVCF(sample, chrom, mut_type, total_depth, alt_depth, pos, v.ALT))
                    
    
    def vcf_translation(self):
        """
        Translate vcf data into the CommonVCF type format
        """
        for vcf in self.variants_list:
            for key, values in vcf.vcf_file.items():
                for v in values:
                    cigar_data = None
                    cigar_trasformed = None
                    if "," in v.row.ALT:
                        cigar_data = self.split_vcf_record(v)
                    for form in v.FORMAT.values():
                        alt_depth = int(form["AD"].split(",")[-1])
                        total_depth = int(form["DP"])
                        mut_type = v.INFO["TYPE"].upper()
                        pos = int(v.row.POS)
                        chrom = v.row.CHROM
                        sample = vcf.sample_name
                        if cigar_data is None:
                            cigar_data = [v.row.REF, v.row.ALT, alt_depth, v.INFO["CIGAR"], pos]
                        if self.variant_data.get(sample) is None:
                            self.variant_data[sample] = dict()
                        if self.variant_data[sample].get(pos) is None:
                            self.variant_data[sample][pos] = list()     
                        
                        try:
                            if v.INFO["TYPE"].upper() not in {"SNP", "MNP"}:
                                cigar_trasformed = self.cigar_translation(*cigar_data)
                                #print(cigar_trasformed)
                                alt = None # TODO return alt value to fill common VCF format
                                for cigar in cigar_trasformed:
                                    alt = cigar[1]
                                    self.variant_data[sample][pos].append(
                                        CommonVCF(sample, chrom, cigar[2], total_depth, alt_depth, pos, alt))
                            else:
                                self.variant_data[sample][pos].append(CommonVCF(sample, chrom, mut_type, total_depth, alt_depth, pos, v.row.ALT))
                        except ValueError:
                            VCFLogger.logger.critical("Found unhandled situation, multiple variants listed per site")
                            VCFLogger.logger.critical(key)
                            VCFLogger.logger.critical(v)
                            exit(-1)
                        
    
    def split_vcf_record(self, vcf_row):
        """
        Split a vcf record into multiple if multiple entires listed at one site.
        Pos it same, AD format tag differs, Cigar string differs and type
        TODO Currently assuming one sample per a vcf

        Values returned correspond to, alterante, cigar string and depth
        """
        vcf_form = list(vcf_row.FORMAT.keys())
        if len(vcf_form) > 1:
            VCFLogger.logger.critical("More than one sample recorded in VCF file, this behaviour is currently not supported. Bailing out") 
        # 1: for teh AD tag as first value correpsonds to the reference
        data = list(zip([vcf_row.row.REF] * len(vcf_row.row.ALT.split(",")), vcf_row.row.ALT.split(","), vcf_row.FORMAT[vcf_form[0]]["AD"].split(",")[1:], vcf_row.INFO["CIGAR"].split(","), [vcf_row.row.POS] * len(vcf_row.row.ALT.split(","))))
        return data

    def split_cigar_string(self, cigar_op):
        """
        Split a cigar string into an order set of tuples defining operations then values.
        each op is only one character and they come after the value defining them
        """
        number_positions = [k for k, v in enumerate(cigar_op) if v.isalpha()]
        prev = 0
        ops = []
        for v in number_positions:
            ops.append((cigar_op[prev:v + 1][-1], int(cigar_op[prev:v + 1][:-1])))
            prev = v + 1
        return ops

    def apply_cigar_ops(self, cigar, ref, alt, pos):
        """
        Cigar, reference then alterante
        """

        cigar_operations = {
            "M": lambda *args: args[3][args[4]:args[4]+args[1]],
            "X": lambda *args: args[3][args[4]:args[4]+args[1]],
            "I": lambda *args: args[3][args[4]:args[4]+args[1]],
            #"D": lambda *args: "-" * args[1]
            "D": lambda *args: args[2][args[4]:args[4]+args[1]]
        }

        # args1 = cigar op, arg 2 = op length, arg3 ref_pos 
        position_ops = {
            "M": lambda *args: int(args[1]), # no ref increase
            "X": lambda *args: int(args[1]), # increase ref,
            "I": lambda *args: 0, # no position increase
            "D": lambda *args: int(args[1])
        }
        return_values = list()

        new_alt = str()
        running_pos = 0
        ref_pos = int(pos)
        ref_poses = [ref_pos]
        alt_slices = list()
        for i in cigar:
            alt_slices.append(cigar_operations[i[0]](i[0], i[1], ref, alt, running_pos))
            new_alt += alt_slices[-1]
            ref_pos += position_ops[i[0]](i[0], i[1], ref_pos)
            ref_poses.append(ref_pos)
            if i[0] != "D":
                running_pos += i[1]

        data = [i for i in zip(ref_poses, alt_slices, [i[0] for i in cigar])]
        
        variable_tag = {
            "I": ("+", "INS"),
            "X": ("", "SNP"),
            "D": ("-", "DEL"),
        }
        for i in data:
            if i[2] != "M":
                #print(i[0], f"{variable_tag[i[2]]}{i[1]}")
                #POS, IVAR FORMAT
                return_values.append((i[0], f"{variable_tag[i[2]][0]}{i[1]}", variable_tag[i[2]][1]))
        return return_values

        

    def cigar_translation(self, *args, **kwargs):
        """
        Use input cigar string to convert deletions and insertions too reveiling the proper change. This
        only is needed for freebayes deletions.
        """
        # CIGAR Operations: M-match D-Deletion I-Insertion N-Alignment gap X-Mismatch
        # number matching cigar op preceded cigar op, args[4] is pos
        #TODO need to add kwargs and args here
        cigar_transformed = list()
        if type(args[0]) == tuple:
            for arg in args:
                cigar_ops = self.split_cigar_string(arg[3])
                cigar_transformed.extend(self.apply_cigar_ops(cigar_ops, arg[0], arg[1], arg[4]))
                #print(cigar_ops, arg[0])
        else:
            cigar_ops = self.split_cigar_string(args[3])
            cigar_transformed.extend(self.apply_cigar_ops(cigar_ops, args[0], args[1], args[4]))
        return cigar_transformed
              


if __name__=="__main__":
    #for i in glob.glob("tests/*.vcf"):
    #    vcf_file_ = ReadVCF(i)
    #    for f in vcf_file_.vcf_file:
    #        for g in vcf_file_.vcf_file[f]:
    #            for j in g.FORMAT.values():
    #                #print(j)
    #                print(i, int(j["AD"].split(",")[1]) / int(j["DP"]), g.row.ALT, g.row.REF, g.INFO["TYPE"].upper(), g.row.POS, j["AD"].split(",")[1], j["DP"])
    
    #for key in glob.glob("tests/*.tsv"):
    #    ivar_file = ReadIvar(key)
    #    print(ivar_file.vcf_info)
    #ivar_file = ReadIvar("tests/22_WPG17_NE_0426_2.tsv")
    vcf_data = [ReadVCF(key) for key in glob.glob("tests/*.vcf")]
    #vcf_data = [ReadIvar(key) for key in glob.glob("tests/*.tsv")]
    #ivar_data.append(ReadVCF("tests/22_WPG35_SE_0830_1.ivar_trim.sorted.vcf"))
    VariantData(vcf_data)
    
    #ivar_file.print_ivar_info()
    #ivar_file.create_contgous_entries()