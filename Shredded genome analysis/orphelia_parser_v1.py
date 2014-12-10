"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Script to parse ORF caller - Orphelia's results
"""

import re
import argparse
import os.path

class orphelia_parser(object):
    
    virusName = ""
    
    def __init__(self, fileOrphelia, path):
        self.fileOrphelia = fileOrphelia
        self.path = path
    
    def get_virus_name(self):
        ## extracting the virus name
        fileName = str(self.fileOrphelia)
        fileNameArr = fileName.split(".")
        self.virusName = fileNameArr[0].split("_")[0]
        
    def analyzing_Orphelia_Output(self):
        fOrphelia = open(self.fileOrphelia,"r")
        fOrpheliaOutput = os.path.join(self.path, self.virusName + "_Orphelia_shredded_output.fasta")
        fOrphelia_output = open(fOrpheliaOutput, "w")
        lines = fOrphelia.readlines()
        for line in lines:
            if "rand" not in line:
                if '-' in line:
                   
                    identifierArray = line.split('-')
                    #print identifierArray
                    midIdentifier = ''.join(identifierArray[0]).split("_")
                    midIdentifier.reverse()
                    #print midIdentifier
                    identifier = '_'.join(midIdentifier[1:-1])
                    #print identifier
                    first = ''.join(list(identifierArray[1])[5:])
                    #print first
                    newLine = ">" + first.strip() + '_' + identifier + " Tool=Orphelia " + '-' + ''.join(list(identifierArray[1])[0:2])+"\n"
                    #newLine1 = newLine.strip()
                    #print newLine1
                    fOrphelia_output.write(newLine)
                    
                else:
        
                    identifierArray = line.split('+')
                    #print identifierArray
                    midIdentifier = ''.join(identifierArray[0]).split("_")
                    #print midIdentifier
                    identifier = '_'.join(midIdentifier[1:-1])
                    #print identifier
                    first = ''.join(list(identifierArray[1])[5:])
                    #print first
                    newLine = ">" + first.strip() + '_' + identifier + " Tool=Orphelia " + '+' + ''.join(list(identifierArray[1])[0:2])+"\n"
                    #print newLine1
                    fOrphelia_output.write(newLine)
            
            else:
                pass
            
        fOrphelia_output.close()
        return

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 1 input file \n Orphelia fasta", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =1 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )
    
    args = parser.parse_args()
    fileOrphelia = args.filename[0]
    output_directory = args.output_directory  
    if output_directory is None:
        path=""
    else:
       path = output_directory
       
    op = orphelia_parser(fileOrphelia, path)
    op.get_virus_name()
    op.analyzing_Orphelia_Output()

if __name__ == '__main__':
    main()
