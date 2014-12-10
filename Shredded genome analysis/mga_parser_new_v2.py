from sys import argv
import re
import argparse
import os.path


"""
Modified on 8/6/2014
This script has correction to invert co-ordinates if there is a - sign.
It wasn't considered earlier.
This is the latest MGA shredded parser script.

9/16/2014
Commented out indicator related logic (-2, -1)
Exact matched results increased a lot.
"""

class mga_parser(object):
    
    virusName = ""
    def __init__(self, mgaFile, path):
        self.mgaFile = mgaFile
        self.path = path
        
    def get_virus_name(self):
        ## extracting the virus name
        #fileName = self.mgaFile.split("//")[-1]
        fileNameArr = self.mgaFile.split(".")
        print fileNameArr
        self.virusName = fileNameArr[0].split("_")[0]
        
    def parsing_mga_shredded_input(self):
        fMgaFile = open(self.mgaFile,"r")
        fileMGAoutput = os.path.join(self.path, self.virusName + "_MGA_shredded_output.fasta")
        fMgaOutput = open(fileMGAoutput,"w")
        lines = fMgaFile.readlines()
        count=0
        for line in lines:
            
            first, second, indicator = (0,0,0)
            if line.startswith("# gi"):
                line = line.strip()
                readName = line.replace("# ",">")
            elif(line.startswith("gene")):
                lineArr = line.split("\t")
                sign = lineArr[3]
                indicator = lineArr[4]
                first = int(lineArr[1])
                second = int(lineArr[2])
                #if indicator == "1":
                #    second = second - 2
                #if indicator == "2":
                #    second = second - 1
                if(sign == "-"):
                    fMgaOutput.write(readName + "_" + str(second) + "_" + str(first) + " Tool=MetaGeneAnnotator" + "\n")
                else:
                    fMgaOutput.write(readName + "_" + str(first) + "_" + str(second) + " Tool=MetaGeneAnnotator" + "\n")
            else:
                pass

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 1 input file \n MGA shredded input", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =1 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )

    args = parser.parse_args()
    mgaFile = args.filename[0]
    
    output_directory = args.output_directory
    if output_directory is None:
        path=""
    else:
       path = output_directory
    
    mp = mga_parser(mgaFile, path)
    mp.get_virus_name()
    mp.parsing_mga_shredded_input()
    
if __name__ == '__main__':
    main()