from sys import argv
import argparse
import os.path
import re

class mga_parser_with_len(object):
    
    def __init__(self, mga_original_file, cutoff):
        self.mga_original_file = mga_original_file
        self.cutoff = cutoff
        
    def parse_original_mga_file(self):
        fMga_original_file = open(self.mga_original_file,"r")
        fileMGAoutput = os.path.join("MGA_shredded_output_with_" + str(self.cutoff) + ".fasta")
        fMgaOutput = open(fileMGAoutput, "w")
        self.cutoff = float(self.cutoff)
        
        lines = fMga_original_file.readlines()
        for line in lines:
            if line.startswith("# gi"):
                line = line.strip()
                readName = line.replace("# ",">")
            elif(line.startswith("gene")):
                lineArr = line.split("\t")
                first = int(lineArr[1])
                second = int(lineArr[2])
                score = float(lineArr[6])
                sign = lineArr[3]
                length = second - first + 1
                length = int(length)
                if score > self.cutoff:
                    if(sign == "-"):
                        fMgaOutput.write(readName + "_" + str(second) + "_" + str(first) + " Tool=MetaGeneAnnotator " + "\n")
                    else:
                        fMgaOutput.write(readName + "_" + str(first) + "_" + str(second) + " Tool=MetaGeneAnnotator " + "\n")
            else:
                pass
        
        fMgaOutput.close()

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 1 input files \n 1. Orginal MGA file \n 2. Cutoff Score", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =2 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )

    args = parser.parse_args()
    mga_original_file = args.filename[0]
    cutoff = args.filename[1]
    mga = mga_parser_with_len(mga_original_file, cutoff)
    mga.parse_original_mga_file()
    
if __name__ == '__main__':
    main()
