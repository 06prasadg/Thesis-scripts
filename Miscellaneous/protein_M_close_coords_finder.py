"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
This script considers co-ordinates of 3 match, 5 match which have their second co-ords different than annotation by  + - 1, + - 2.
Finds proteins of these reads using proteins file and creates 3 separate files - Protein who satisfy above condition and starts with M,
end with M, and a statistics file which gives count of starts, end letters of proteins
"""

from sys import argv
import re
import argparse
import os.path

class protein_finder(object):
    result_dict = {}
    protein_dict = {}
    dna_dict = {}
    startswith_dict = {}
    endswith_dict = {}
    
    def __init__(self, path, resultFile, proteinFile, dnaFile):
        self.path = path
        self.resultFile = resultFile
        self.proteinFile = proteinFile
        self.dnaFile = dnaFile

    def making_result_dict(self):

        fResultFile = open(self.resultFile, "r")
        lines = fResultFile.readlines()
        for line in lines:
            
            annotationArray = line.split(" ")[-1]
            annotation_first = int(annotationArray.split("_")[1])
            annotation_second = int(annotationArray.split("_")[2].strip())
            
            toolArray = line.split(" ")[0].split("/")[1]
            tool_first = int(toolArray.split("_")[1])
            tool_second = int(toolArray.split("_")[2])   
            read = line.split(" ")[0].strip()
            
            if((annotation_first == tool_first + 1 and annotation_second == tool_second) \
                or (annotation_first == tool_first - 1 and annotation_second == tool_second) \
                or (annotation_second == tool_second + 1 and annotation_first == tool_first) \
                or (annotation_second == tool_second - 1 and annotation_first == tool_first) \
                or (annotation_first == tool_first + 2 and annotation_second == tool_second) \
                or (annotation_first == tool_first - 2 and annotation_second == tool_second) \
                or (annotation_second == tool_second + 2 and annotation_first == tool_first) \
                or (annotation_second == tool_second - 2 and annotation_first == tool_first)):
        
                self.result_dict[read]=annotationArray
        #print self.result_dict

                
    def finding_protein(self):
        fOuputFileStarts = open("Protein_M_Output_starts.fasta", "w")
        fOuputFileEnds = open("Protein_M_Output_ends.fasta", "w")
        
        fProteinFile = open(self.proteinFile, 'r')
        fDnaFile = open(self.dnaFile, 'r')
        linesProtein = fProteinFile.readlines()
        read = ""
        protein = ""
        dna = ""
        for line in linesProtein:
            if line.startswith(">gi"):
                read = line.strip()
            else:
                protein = line.strip()
            self.protein_dict[read] = protein
        fProteinFile.close()
        
        linesDna = fDnaFile.readlines()
        for line in linesDna:
            if line.startswith(">gi"):
                read = line.strip()
            else:
                dna = line.strip()
            self.dna_dict[read] = dna
        fDnaFile.close()
        
        for key, value in self.result_dict.items():
            if self.protein_dict.has_key(key):
                protein = self.protein_dict.get(key)
                dna = self.dna_dict.get(key)
                
                start_letter = protein[0]
                end_letter = protein[-1]
                if self.startswith_dict.has_key(start_letter):
                    count = self.startswith_dict[start_letter]
                    self.startswith_dict[start_letter] = count+1
                else:
                    self.startswith_dict[start_letter] = 1
                    
                if self.endswith_dict.has_key(end_letter):
                    count = self.endswith_dict[end_letter]
                    self.endswith_dict[end_letter] = count+1
                else:
                    self.endswith_dict[end_letter] = 1
                        
                if (start_letter == "M"):
                    fOuputFileStarts.write(key + " " + value + protein + "\n" + dna + "\n")    
                elif(end_letter == "M"):
                    fOuputFileEnds.write(key + " " + value + protein + "\n" + dna +"\n")
                else:
                    pass
                    
    def print_dicts(self):
        fOuputFileStats = open("Protein_letter_statistics.txt","w")
        
        fOuputFileStats.write("Protein startswith statistics : \n")
        for k, v in self.startswith_dict.items():
            fOuputFileStats.write(k + " " + str(v) + "\n")
        fOuputFileStats.write("\n")
        fOuputFileStats.write("Prorein endswith statistics : \n")
        for k, v in self.endswith_dict.items():
            fOuputFileStats.write(k + " " + str(v) + "\n")
 

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 2 input files \n tool results file \
            \n 2. tool protein file \n 3. tool DNA file \n", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =3 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )
    
    args = parser.parse_args()
    resultFile = args.filename[0]
    proteinFile = args.filename[1]
    dnaFile = args.filename[2]
    
    output_directory = args.output_directory
    if output_directory is None:
        path=""
    else:
       path = output_directory
       
    pf = protein_finder(path, resultFile, proteinFile, dnaFile)
    pf.making_result_dict()
    pf.finding_protein()
    pf.print_dicts()
    
if __name__ == '__main__':
    main()
