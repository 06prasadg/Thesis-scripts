"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Script to parse ORF caller - MGM results
"""

from sys import argv
import re
import argparse
import os.path

   
class MGM_parser(object):
    
    def __init__(self,fileMGM):
        self.fileMGM = fileMGM

    def analyzing_MetaGeneMark_output(self):
        fMgm = open(self.fileMGM,"r")
        #fileMGMoutput = fileMGM.split (".")[0] + "_output.fasta"
        fileMGMoutput = "MGM_gff_converted_output_with_changes.fasta"
        fMgm_output = open(fileMGMoutput ,"w")
        lines = fMgm.readlines()
        firstNum, secondNum = (0,0)
        for line in lines:
            if line.startswith("gi|403225042|ref|NC_017971.2"):
                lineArr = line.split("\t")
                firstNum = int(lineArr[3])
                secondNum = int(lineArr[4])
                
                ## changes in 9/16/2014
                if(firstNum == 3 and secondNum == 299):
                    firstNum = 1
                    secondNum = 300
                elif(firstNum == 299 and secondNum == 3):
                    firstNum = 300
                    secondNum = 1
                elif(firstNum == 2 and secondNum == 298):
                    firstNum = 1
                    secondNum = 300
                elif(firstNum == 298 and secondNum == 2):
                    firstNum = 300
                    secondNum = 1
                
                ## making changes to make them match to output mapped
                #if(firstNum == 3 and secondNum == 299):
                #    firstNum = firstNum-2
                #elif(firstNum == 2 and secondNum == 298):
                #    firstNum = firstNum-1
                    
                ## doing +2 297_1 => 299_1
                #elif(firstNum == 1 and secondNum == 297):
                #    secondNum = secondNum+2
                ## doing +1 297_1 => 298_1
                elif(firstNum == 1 and secondNum == 297):
                   secondNum = secondNum+1
                    
                readname = lineArr[0]
                if lineArr[6] == "-":
                    newLine = ">" + readname+"_" + str(secondNum) + "_" + str(firstNum) + " Tool=MetaGeneMark\n"
                else:
                    newLine = ">" + readname+"_" + str(firstNum) + "_" + str(secondNum) + " Tool=MetaGeneMark\n"
                fMgm_output.write(newLine)  
        fMgm_output.close()
                               

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 1 input file \n MetaGeneMark fasta", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =1 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )
    
    args = parser.parse_args()
    fileMGM = args.filename[0]
    
    output_directory = args.output_directory
    if output_directory is None:
        path=""
    else:
       path = output_directory
    mgmObject = MGM_parser(fileMGM)
    mgmObject.analyzing_MetaGeneMark_output()
    
if __name__ == '__main__':
    main()
