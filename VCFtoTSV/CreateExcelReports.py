"""
Create the excel reports of the vcf data used in creating the output
of the files. This will slow down the code but they have all asked for it.
It will add in dependecies so defining this method function outside of the
class so openpyxl is not always needed.

With great sadness:
2022-05-31: Matthew Wells
"""

# imports I didnt want to bring in, but makes my life easier
# Note the quality of libraries is great, I just didnt want to have to add dependencies
#import openpyxl
import pandas as pd

# Nice python library 
import os
import glob



class HTMLToExcel:
    """
    Read all html tables in a directory and covert to a single excel file
    """
    _glob_word_ = "test" # the word terminating the html plots generated currently
    def __init__(self, directory_recurse) -> None:
        self.html_voc_data = {}
        self.directory_recurse = directory_recurse
        self.directories = self.recurse_directory()
        [self.read_html_file(i) for i in self.directories] # map call wasnt working for somereason

    def read_html_file(self, fp):
        """
        Read html file and 
        """
        #fp = f"./tests/Delta_{self._glob_word_}.html"
        sample_name = os.path.basename(fp)[:os.path.basename(fp).index("_")]
        df = pd.read_html(fp)
        samples = [i for i in list(df[0].columns)[1:]]
        df[0].rename(columns={"Unnamed: 0": "NucName+AAName"}, inplace=True)
        df[0]["AAName"], df[0]["NucName"] = zip(*df[0]["NucName+AAName"].apply(lambda x: x.split("|"))) # create split columns for indexing
        df[0]["VOC"] = sample_name
        df[0]["Position"] = df[0]["NucName"].apply(lambda x: int(''.join([i for i in x if i.isdigit()])))
        df[0] = df[0].reindex(columns=["VOC", "Position", "AAName", "NucName", "NucName+AAName", *samples])
        self.html_voc_data[sample_name] = df[0] # update and maintain a worksheets in place
        print(df[0].columns)
        return True
    
    def recurse_directory(self):
        """
        Find all HTML files for voc data to be read in from.
        """
        glob_str = os.path.join(self.directory_recurse, f"*_{self._glob_word_}.html")
        files = glob.glob(glob_str)
        return files
    
    def html_dict_to_excel(self):
        """
        Take in the
        """
        ...


if __name__=="__main__":
    xx = HTMLToExcel("./tests")
    print(xx.html_voc_data.keys())