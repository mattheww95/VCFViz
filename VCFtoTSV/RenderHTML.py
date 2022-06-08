"""
Take the vcf data and render a heatmap using HTML using CSS.

TODO: add in ability for reversion detection in the future

2022-05-17: Matthew Wells
"""
from collections import namedtuple
from dataclasses import dataclass
from datetime import datetime
import glob
from typing import NamedTuple, List
import os
from VCFtoTSV import CoverageData
from VCFtoTSV.VCFToJson import ReadIvar, ReadVCF, IvarFields
from VCFtoTSV.VCFlogging import VCFLogger as vlog
import statistics




class VCFParserRow(NamedTuple):
        PangoLineage: str
        NextStrainClade: str 
        NucName: str
        AAName: str
        Key: str
        SignatureSNV: str 
        Position: int
        Type: str
        Length: str
        Ref: str
        Alt: str


class DataSheet:
    """
    Read in the vcfparser sheet to decide what mutations to index.
    Columns of used and already know include:
        - VOC (string)
        - Position (Int)
        - Type (Sub, Del, Ins, MNP) Dels and sub, index differently
        - Ref
        - Alt
    The sheet is tab delimited
    """

    def __init__(self, file_name):
        self.file_name = file_name
        self.voc_info = self.group_voc_info()

    def group_voc_info(self):
        """
        From the vcfparser sheet file group the needed data
        """    
        voc_info = {}
        with open(self.file_name, 'r') as metadata:
            lines = metadata.readlines()[1:] # header fields stored in named tuple
            for line in lines:
                line = line.strip().split("\t") # split on tab delimiter remove new lines
                voc = line[0]
                voc_metadata = VCFParserRow(*line[1:])
                if voc_info.get(voc) is None:
                    voc_info[voc] = {}
                
                voc_info[voc][f"{voc_metadata.Ref}{voc_metadata.Position}{voc_metadata.Alt}"] = voc_metadata
        return voc_info

class PlotData(NamedTuple):
    metadata: VCFParserRow
    ivar_row: IvarFields
    sample_name: str
    sample_depth: int = None

class VCFDataHTML:

    """
    Take in the position data and create a table an html table.
    Starting witha single sample as a proof of concept

    2022-05-17: Matthew Wells
    """
    __slots__ = [
        "ivar_data", "__dict__", "indx_samples"]

    html_meta_start = """
    <!DOCTYPE html>
    <html>
    <head>
    <style>
    html {
        zoom: 0.25
    }
    .banner {
        font-family: \"Arial\";
        font-size: 80px;
        height: 100%;
        margin-bottom: 400px;
        padding: 60px;
        text-align: center;
    }
    table, th, tr, td {
      border: 1px solid black;
      border-collapse: collapse;
      text-align: center;
      font-family: \"Arial\";
      font-size: 20px;
      white-space: nowrap;
    }
    h1 {
        font-family: \"Arial\";
        font-size: 80px;
    }
    thead {
        font-size: 40px;
    }
    /*side nave code copy and paster from https://www.w3schools.com/howto/howto_css_fixed_sidebar.asp */
    .sidenav {
      height: 100%; /* Full-height: remove this if you want "auto" height */
      width: 800px; /* Set the width of the sidebar */
      position: fixed; /* Fixed Sidebar (stay in place on scroll) */
      z-index: 1; /*  Stay on top */
      top: 0; /* Stay at the top */
      left: 0;
      background-color: #111; /* Black */
      overflow-x: hidden; /* Disable horizontal scroll */
      font-family: \"Arial\";
      padding-top: 20px;
      margin-right: 500px;
      font-size: 40px;
    }

    /* The navigation menu links */
    .sidenav a {
      padding: 6px 8px 6px 16px;
      text-decoration: none;
      font-size: 80px;
      color: coral;
      display: block;
    }
    .main {
        margin-left: 900px;
    }
    </style>
    </head>
    <body>
    """
    #out_dir = "/mnt/c/Users/mwells/Desktop/NewVCFParser_tests/"
    banner_tags = ["<section class=\"banner\">", "</section>"]
    side_bar_nav = ["<div id=\"sidenav_\" class=\"sidenav\">", "</div>"]
    nav_element = ["<a href=\"#@\">", "</a>"]
    html_meta_end = """
    </body>
    </html>
    """
    final_tag = "test"
    td_tags = ("<td>", "</td>") # data tags
    row_tags = ("<tr>", "</tr>") # row tags wrapping headers and data
    header_tags = ["<th>", "</th>"]
    info_header_tag = ["<h1>", "</h1>"]
    table_tags = ("<table style=\"border:1px solid black;margin-left: 0px\">", "</table>")
    CSS_colours = ["#DEEDCF",
                    "#BFE1B0",
                    "#99D492",
                    "#74C67A",
                    "#56B870",
                    "#39A96B",
                    "#1D9A6C",
                    "#188977",
                    "#137177",
                    "#0E4D64",
                    "#0A2F51",
                    ]
    css_text_colour = "coral"

    def __init__(self, ivar_data: List[ReadIvar], vcf_parser_sheet: str, search_dir: str, cov_thresh: int, out_dir: str, prep_cov_data = None) -> None:
        """
        TODO: have flag for coverage info so that it can run without it
        Can be done better for handing off data, but just to rush out a prototype, e.g. not just ivar specific
        prep_cov_data: is a parameter to be added in the case of preprocessed data is provided
        """
        self.out_dir = out_dir
        self.vcfparser_sheet = vcf_parser_sheet
        self.low_cov_thresh = cov_thresh #TODO make this a param in cmd line
        self.indx_samples = {i.sample_name: i for i in ivar_data}
        if prep_cov_data == None:
            self.cov_info = CoverageData.create_sample_coverages([i.sample_name for i in ivar_data], search_dir)
        else:
            self.cov_info = prep_cov_data
        self.vcf_metadata = DataSheet(vcf_parser_sheet)
        #self.voc_table_data = self.initialize_voc_tables()
        self.ivar_data = ivar_data
        self.figure_data = {}
        for data in self.ivar_data:
            # modifies figure data obj in place, adding in data for figures
            self.figure_data = self.initialize_voc_tables(data, self.figure_data)
        self.create_heatmaps()
    
    def initialize_voc_tables(self, datafile, html_plots: dict):
        """
        From the vcf metadata initialize a dictionary for each voc that can show a queried postion,
        and can be grabbed from the vcf_data

        TODO: methods can be divided and mutations cleaned with gaurd statements
        """

        for key in self.vcf_metadata.voc_info.keys():
            if html_plots.get(key) is None:
                html_plots[key] = {}

            for voc in self.vcf_metadata.voc_info[key]:
                #TODO double value addition is probably happening here
                pos = "".join([i for i in voc if i.isdigit()])
                if html_plots[key].get(voc) is None:
                    html_plots[key][voc] = []
                ivar_data = datafile.vcf_info.get(pos)
                sample_name = datafile.sample_name
                # add in coverage check here
                depth = self.cov_info.samples_coverage[datafile.sample_name][self.vcf_metadata.voc_info[key][voc].Position]
                empty_data = PlotData(self.vcf_metadata.voc_info[key][voc], None, sample_name, int(depth))

                if ivar_data is None:
                    vlog.logger.debug(f"The no data found VOC {key} mutations {voc}")
                    html_plots[key][voc].append(empty_data)
                else:
                    start_len = len(html_plots[key][voc])
                    for i in ivar_data:
                        ret_val = self.append_ivar_info(html_plots, i, key, voc, datafile)
                    if not ret_val and len(html_plots[key][voc]) == start_len:
                        html_plots[key][voc].append(empty_data) # add empty value if data could not be found
        return html_plots
    
    def append_ivar_info(self, html_plots_obj: dict, ivar_data_val, key_val, voic, datafile):
        """
        As ivar data contains a list of values, a for loop is required to run through the data 
        of tuples to identify other postitions. No return value is specified as the dictionary will be mutatated in place
        """
        
        depth = self.cov_info.samples_coverage[datafile.sample_name][self.vcf_metadata.voc_info[key_val][voic].Position]
        plot_data = PlotData(self.vcf_metadata.voc_info[key_val][voic], ivar_data_val, datafile.sample_name, int(depth))
        vcf_data_meta = self.vcf_metadata.voc_info[key_val][voic]
        # to compare indels, vcf parser sheet places ref at front
        if vcf_data_meta.Type != "Sub":
            if vcf_data_meta.Type == "Del":
                # splitting string to rwmove first char as in vcfparser
                # we include the ref codon and ivar includes a starting -
                ivar_del = ivar_data_val.ALT[1:] 
                meta_del = vcf_data_meta.Ref[1:]
                test = ivar_del == meta_del
                if not test:
                    vlog.logger.info(f"Mismatch in deletion from metadata: {meta_del} and VCF deletion {ivar_data_val.ALT}")
                    return False
                else:
                    html_plots_obj[key_val][voic].append(plot_data)
            elif vcf_data_meta.Type == "Ins":
                ivar_ins = ivar_data_val.ALT[1:]
                meta_ins = vcf_data_meta.Alt[1:]
                test = ivar_ins == meta_ins
                if not test:
                    vlog.logger.info(f"Mismatch in insertion from metadata: {meta_ins} and VCF deletion {ivar_data_val.ALT}")
                    return False
                else:
                    html_plots_obj[key_val][voic].append(plot_data)
            elif vcf_data_meta.Type == "Mnp":
                #TODO move cv into static methods
                meta_alt = vcf_data_meta.Alt
                l_meta_alt = len(meta_alt)
                meta_pos = int(vcf_data_meta.Position)
                mnps_add = []
                if ivar_data_val.ALT != meta_alt[0]:
                    return False

                mnps_add.append(ivar_data_val)
                for i in range(1, l_meta_alt):# skip ivar dataval
                    pos_test = datafile.vcf_info.get(str(meta_pos + i))
                    for vcf_val in pos_test: # test one vcfval at a time alter this later
                        if vcf_val.ALT == meta_alt[i]:
                            mnps_add.append(vcf_val)
                
                # Quick lil CV test for the mnps to make sure the occur together
                # NOTE: this is based on my faith and not empircical measures
                alleles = "".join([i.ALT for i in mnps_add])
                if alleles != meta_alt:
                    return False
                depths = [int(i.ALT_DP) for i in mnps_add]
                d_stdev = statistics.pstdev(depths)
                d_avg = statistics.mean(depths)
                depths_cv = (d_stdev / d_avg) * 100

                alt_freqs = [float(i.ALT_FREQ) for i in mnps_add]
                altf_stdev = statistics.pstdev(alt_freqs)
                alt_avg = statistics.mean(alt_freqs)
                alt_cv = (altf_stdev / alt_avg) * 100
                #Change the ivar row in the data to show the average for the mnps
                
                plot_data = plot_data._replace(ivar_row = plot_data.ivar_row._replace(ALT_FREQ = alt_avg))
                plot_data = plot_data._replace(ivar_row = plot_data.ivar_row._replace(ALT = vcf_data_meta.Alt))
                plot_data = plot_data._replace(ivar_row = plot_data.ivar_row._replace(REF = vcf_data_meta.Ref))
                # CV thresholds are set arbitralily, and should reflect the depth
                cov_var = 2.5
                if depths_cv < cov_var and alt_cv < cov_var: # TODO make this calculated based on depth
                    vlog.logger.info(f"Combining {alleles} at position {meta_pos} into MNP")
                    html_plots_obj[key_val][voic].append(plot_data)
                else:
                    vlog.logger.info(f"Could not combine mutations for {voic} due to a Coefficient of Variation greater than {cov_var}.")
                    return False

            else:
                vlog.logger.warning(f"Support not provided for mutation type {vcf_data_meta.Type}")
                return False
        else:
            test_alt = vcf_data_meta.Alt == ivar_data_val.ALT
            test_ref = vcf_data_meta.Ref == ivar_data_val.REF
            query_combo = ivar_data_val.REF + str(ivar_data_val.POS) + ivar_data_val.ALT
            if test_alt and test_ref and self.vcf_metadata.voc_info[key_val].get(query_combo) is not None:
                html_plots_obj[key_val][voic].append(plot_data)
            else:
                vlog.logger.info(f"Mismatch in substitution from metadata: {vcf_data_meta.Ref} at"\
                    f" position {ivar_data_val.POS} and VCF {ivar_data_val.ALT}")
                return False

    def check_alt_prescence(self, combo, voc_key):
        """
        check for prescence of alternate combo in vcfparser data
        """
        val = self.vcf_metadata.voc_info[voc_key].get(combo)
        if val is not None:
            return True
        return False

    def create_heatmaps(self):
        """
        Run the code to create the heatmaps from the intialized data sheets
        """
        figure_data = self.figure_data
        #side_bar_nav = [self.side_bar_nav[0]]
        for data in figure_data.keys():
            html_figure = [self.html_meta_start, f"<h1>{data}</h1>", self.table_tags[0],"<thead>",self.row_tags[0], 
            self.header_tags[0] + "" + self.header_tags[1]]

            voic_data = figure_data[data]
            positions = list(voic_data.keys())

            # replaceing header tags for samples only
            val_headers = [f"<th style=\"transform: rotate(180deg);padding:25px;font-size:50px;writing-mode:vertical-lr;\">{i.sample_name}{self.header_tags[1]}" 
            for i in voic_data[positions[0]]]
            html_figure.extend(val_headers)
            html_figure.append("</thead>")

            for voic in voic_data.keys(): # add figure data
                html_figure.append("<tr style=\"height:200px;\">")
                col_name = voic_data[voic][0].metadata.AAName + "|" + voic_data[voic][0].metadata.NucName # getting first value for row name
                html_figure.append("<td style=\"font-weight: 600;padding: 10px;font-size:50px;\">" + col_name + self.td_tags[1])
                for vcf_row in voic_data[voic]:
                    row_info = []
                    if vcf_row.ivar_row is not None:
                        alt_freq = "LC"
                        colour_css = "#ffffff"
                        if vcf_row.sample_depth >= self.low_cov_thresh:
                            alt_freq = str(round(float(vcf_row.ivar_row.ALT_FREQ), 3)) # get alt freq col
                            colour_css = self.CSS_colours[self.pick_colour(alt_freq=alt_freq)] # get index based on colur scale

                        row_data = f"<td bgcolor={colour_css} style=\"color:{self.css_text_colour};font-weight: 600;padding: 10px;font-size:50px\">"
                        row_info.append(row_data + alt_freq + self.td_tags[1])
                    else:
                        cov = vcf_row.sample_depth
                        #TODO preprocess this to speed it up
                        #TODO make logic prepared for mixture of wildtype and alternate
                        
                        sample_data = self.determine_sample(vcf_row.sample_name)
                        val = sample_data.vcf_info.get(vcf_row.metadata.Position)
                        empt_flag = "WT"
                        if val is not None:
                            empt_flag = "ALT"
                        # = "NC" # get alt freq col
                        if cov == 0:
                            alt_freq = "NC"
                        elif cov <= self.low_cov_thresh:
                            alt_freq = f"{empt_flag}_LC"
                        elif cov > self.low_cov_thresh:
                            alt_freq = f"{empt_flag}"
                        row_data = f"<td style=\"color:{self.css_text_colour};padding: 10px;font-weigth: 600;font-size:50px;\">"
                        row_info.append(row_data + alt_freq + self.td_tags[1])    
                    html_figure.extend(row_info)
                row_info.append(self.row_tags[1])
                       
            html_figure.append(self.table_tags[1])
            html_figure.append(self.html_meta_end)
            html_figure = "\n".join(html_figure)

            with open(os.path.join(self.out_dir, f"{data}_{self.final_tag}.html"), "w") as html_out:
                vlog.logger.info(f"Creating plot for {data}")
                html_out.write(html_figure)
    
    def combine_html_plots(self):
        """
        Create a naviagable webpage of the output plots
        """
        glob_pattern = os.path.join(self.out_dir, f"*{self.final_tag}.html")
        plots = glob.glob(glob_pattern)
        html_pages = []
        headers = []
        header_tag = self.info_header_tag[0]
        vlog.logger.info("Creating combined Report")
        for plot in plots:
            with open(plot, 'r') as cur_plot:
                lines = cur_plot.read()
                #headers.append(lines[lines.index(header_tag)+ len(header_tag):lines.index("</h1>")])
                header = lines[lines.index(header_tag)+ len(header_tag):lines.index("</h1>")]
                headers.append(header)
                html_pages.append(f"<a id=\"{header}\"></a>")
                html_pages.append(lines[lines.index(header_tag):lines.rindex("</body>")].strip())

        headers_formatted = [f"{self.nav_element[0].replace('@', i)}{i}{self.nav_element[1]}" for i in headers]
        banner_formatted = [self.banner_tags[0], 
        f"<h1>VCFParser sheet used: {self.vcfparser_sheet}</h1>", 
        f"<h1>Time combined report created: {datetime.now()}</h1>",
        f"<h1>Depth of Coverage Threshold: {self.low_cov_thresh}</h1>", 
        self.banner_tags[1]]
        side_bar_nav = [self.side_bar_nav[0]]
        side_bar_nav.extend(headers_formatted)
        side_bar_nav.append(self.side_bar_nav[1])
        html_doc = [self.html_meta_start, *banner_formatted, *side_bar_nav, "<div class=\"main\">"]
        html_doc.extend(html_pages)
        html_doc.append("</div>")
        html_doc.append(self.html_meta_end)
        html_fig = "\n".join(html_doc)
        with open(os.path.join(self.out_dir, f"{os.path.basename(self.out_dir)}_doc.html"), "w") as combined:
            combined.write(html_fig)


    def determine_sample(self, sample_name):
        """
        get the ivar information for the sample
        """
        return self.indx_samples.get(sample_name)


    @staticmethod
    def pick_colour(alt_freq: float):
        return int(float(alt_freq) * 10)

    def create_ivar_fig(self):
        """
        This creates heat maps for ivar data (TEST ONLY CURRENTLY)
        """
        header_cols = [self.html_meta_start, self.table_tags[0], self.row_tags[0],
        self.header_tags[0] + "POS" + self.header_tags[1], 
        self.header_tags[0] + self.ivar_data.sample_name + self.header_tags[1],
        self.row_tags[1]]
        row_info = []
        for key in self.ivar_data.vcf_info.keys():
            row_info.append(self.row_tags[0]) # append in row delimeters
            row_info.append(self.td_tags[0] + key + self.td_tags[1])
            # need to assign some logic to tags below
            row_data = self.td_tags[0] # add Make copy to add css to tag
            alt_freq = self.ivar_data.vcf_info[key][0].ALT_FREQ # get alt freq col
            # below colour logic can be altered to scale with colours and mathmatically cleaner
            colour_css = int(float(alt_freq) * 10) # get index based on colur scale
            row_data = f"<td bgcolor={self.CSS_colours[colour_css]} style=\"color:red;\">"
            row_info.append(row_data + alt_freq + self.td_tags[1])
            row_info.append(self.row_tags[1])

        header_cols.extend(row_info)
        header_cols.append(self.table_tags[1])
        header_cols.append(self.html_meta_end)
        html_data = "\n".join(header_cols)
        with open("out_test.html", "w") as html_out:
            html_out.write(html_data)


if __name__ == "__main__":
    ivar_data_ = ReadIvar("tests/22_AB16_GP_0414.tsv")
    ivar_data2 = ReadIvar("tests/22_WPG17_NE_0426_2.tsv")
    #vcf_html = VCFDataHTML([ivar_data2, ivar_data_], "tests/VCFParser_20220516.txt")
    #vcf_html.combine_html_plots()
    parser_sheet = "tests/VCFParser_20220516.txt"
    import InputOptions
    #InputOptions.glob_directories(t_ivar, bam_dir, parser_sheet, cov, out_dir)





