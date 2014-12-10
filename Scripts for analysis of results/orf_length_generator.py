"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Generates two separate files with reads below and above a certain ORF length cutoff
Also provides the count for each category
"""

from sys import argv
import argparse
import os.path
import re

class orf_len_calculator(object):
    extra_annot_dict = {}
    
    def __init__(self, extra_annot_file, cutoff):
        self.extra_annot_file = extra_annot_file
        self.cutoff = cutoff
                       
    def parse_extra_annot_file(self):
        fExtraAnnot = open(self.extra_annot_file,"r")
        lines = fExtraAnnot.readlines()
        less_than_cutoff = 0
        more_than_cutoff = 0
        for line in lines:
            lineSplit = line.split("/")
            readName = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
            matchObject = re.findall(r'_+[0-9]*', lineSplit[1])
            value = ''.join(matchObject[0:2])
            valueArr = value.split("_")
            num1 = int(valueArr[1])
            num2 = int(valueArr[2])
            if num1 > num2:
                length = num1 - num2 + 1
            else:
                length = num2 - num1 + 1
            #print length
            if length <= int(self.cutoff):
                less_than_cutoff = less_than_cutoff + 1
            else:
                more_than_cutoff = more_than_cutoff + 1
             
            length = str(length)    
            if self.extra_annot_dict.has_key(readName):
                listItems = self.extra_annot_dict.get(readName)
                if (isinstance(listItems,str)):
                    if not (listItems == length):
                        listItems = []
                        listItems.append(self.extra_annot_dict.get(readName))
                        listItems.append(length)
                        self.extra_annot_dict[readName] = listItems
                else:
                    for l in listItems:
                        if not(l == length):
                           listItems.append(length)
                           self.extra_annot_dict[readName] = listItems
            else:
                self.extra_annot_dict[readName] = length
            
        print "No of reads less than " + str(self.cutoff) + " : ", less_than_cutoff
        print "No of reads more than " + str(self.cutoff) + " : ", more_than_cutoff

    def write_output(self):
        fOutput = open("Annotation_reads_length.csv","w")
        
        for key, value in self.extra_annot_dict.items():
            length = self.extra_annot_dict.get(key)
            if (isinstance(length,str)):
                if int(length) <= int(self.cutoff):
                    flag = 0
                else:
                    flag = 1
                fOutput.write(key + "," + str(length) + "," + str(flag) + "\n")
            else:
                for l in length:
                    if int(l) <= int(self.cutoff):
                        flag = 0
                    else:
                        flag = 1
                    fOutput.write(key + "," + str(l) + "," + str(flag) + "\n")
        
    
def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 1 input file \n1. Extra annotation file of tool\n 2. Cutoff length", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =2 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )

    args = parser.parse_args()
    extra_annot_file = args.filename[0]
    cutoff = args.filename[1]
    olc = orf_len_calculator(extra_annot_file, cutoff)
    olc.parse_extra_annot_file()
    olc.write_output()
    
if __name__ == '__main__':
    main()
