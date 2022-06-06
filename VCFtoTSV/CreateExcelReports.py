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



class HTMLToExcel:
    """
    Read all html tables in a directory and covert to a single excel file
    """
    def __init__(self, directory_recurse) -> None:
        self.directory_recurse = directory_recurse
    
    def read_html_file(self):
        """
        Read html file and 
        """
        df = pd.read_html("./tests/Delta_test.html")
        print(type(df))


if __name__=="__main__":
    xx = HTMLToExcel("./tests")
    xx.read_html_file()