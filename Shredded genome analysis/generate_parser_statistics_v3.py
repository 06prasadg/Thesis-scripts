"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Script to generate statistics from all results

9/23/2014
making changes to not print other not imp statistics.
adding split case count
"""

from sys import argv
import argparse
import os.path

class statistics(object):
    
    exact_matched = 0
    partial_matched = 0
    in_tool_but_not_in_annotation = 0
    frame_shift = 0
    total_tool_count = 0
    total_annotation_count = 0
    in_annotation_but_not_in_tool = 0
    total_3_match_count = 0
    total_5_match_count = 0
    extra_annot_count = 0
    no_condition_count = 0
    total_split_cases_count = 0
    
    def __init__(self, result_input, tool_input, annot_input, match_3_input, match_5_input, split_cases_input, no_condition_input, extra_annot_input):
        self.result_input = result_input
        self.tool_input = tool_input
        self.annot_input = annot_input
        self.match_3_input = match_3_input
        self.match_5_input = match_5_input
        self.split_cases_input = split_cases_input
        self.extra_annot_input = extra_annot_input
        self.no_condition_input = no_condition_input
        
    ##reading result file
    def read_result_file(self):
        fResult_input = open(self.result_input,"r")
        lines = fResult_input.readlines()
        for line in lines:
            if (not line.startswith(">rand")):
                lineArr = line.split(" ")
                indicator = lineArr[2]
                if lineArr[2] == "00":
                    self.exact_matched += 1
                if lineArr[2] == "01" or lineArr[2] == "10" or lineArr[2] == "11" or lineArr[2] == "20" or lineArr[2] == "02" or lineArr[2] == "22" or lineArr[2] == "12" or lineArr[2] == "21":
                    self.partial_matched += 1
                if lineArr[2] == "TT":
                    self.in_tool_but_not_in_annotation += 1
                if lineArr[2] == "F-1" or lineArr[2] == "F+1" or lineArr[2] == "F-2" or lineArr[2] == "F+2":
                    self.frame_shift += 1
                if lineArr[2] == "AA":
                    self.in_annotation_but_not_in_tool += 1
            
    def read_tool_file(self):
        ## reading tool file       
        fTool_input = open(self.tool_input,"r")
        for line in fTool_input.readlines():
            if (not line.startswith(">rand")):
               self.total_tool_count += 1
               
    def read_annot_file(self):
        ## reading annotation input
        fAnnot_input = open(self.annot_input,"r")
        for line in fAnnot_input.readlines():
            self.total_annotation_count += 1
    
    def read_3_match_file(self):
        ## reading 3 matched file
        fMatch_3_input = open(self.match_3_input,"r")
        for line in fMatch_3_input.readlines():
            self.total_3_match_count += 1
    
    def read_5_match_file(self):
        ## reading 5 matched file
        fMatch_5_input = open(self.match_5_input,"r")
        for line in fMatch_5_input.readlines():
            self.total_5_match_count += 1
            
    def read_extras_file(self):
        ## reading extra annotation file
        fExtra_annot_input = open(self.extra_annot_input,"r")
        for line in fExtra_annot_input.readlines():
            self.extra_annot_count += 1
    
    def read_no_condition_file(self):
        ## reading No condition file
        fNo_condition_input = open(self.no_condition_input,"r")
        for line in fNo_condition_input.readlines():
            self.no_condition_count += 1
            
    def read_split_cases_file(self):
        ## reading tool file       
        fSplit_cases_input = open(self.split_cases_input,"r")
        for line in fSplit_cases_input.readlines():
            if (not line.startswith(">rand")):
               self.total_split_cases_count += 1
               
    def display_results(self):
        print "Total annotation count : ", self.total_annotation_count
        print "Total Tool count : ", self.total_tool_count
        print "Exact matched : ", self.exact_matched
        print "   Exact matched [(+/-)1/2] (3' matched) : ", self.total_3_match_count
        print "   Exact matched [(+/-)1/2] (5' matched) : ", self.total_5_match_count
        print "   Exact matched [(+/-)1/2] (both sides) : ", self.frame_shift 
        print "Partial matched : ", self.partial_matched
        print "ORFs in tool but not in annotation : ", self.no_condition_count + self.in_tool_but_not_in_annotation
        print "ORFs in Annotation but not in tool: ", self.in_annotation_but_not_in_tool + self.extra_annot_count
        print "Split cases count : ", self.total_split_cases_count 

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 6 input files \n 1. Results file \n 2. Tool input file \n 3. Annotation file \n 4. 3 Matched file \n 5. 5 Matched file \n 6. Split cases file \n 7. No condition file \n 8. Extra annotation file\n", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =8 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )

    args = parser.parse_args()
    result_input = args.filename[0]
    tool_input = args.filename[1]
    annot_input = args.filename[2]
    match_3_input = args.filename[3]
    match_5_input = args.filename[4]
    split_cases_input = args.filename[5]
    no_condition_input = args.filename[6]
    extra_annot_input = args.filename[7]
    
    statisticsObject = statistics(result_input, tool_input, annot_input, match_3_input, match_5_input, split_cases_input, no_condition_input, extra_annot_input)
    statisticsObject.read_result_file()
    statisticsObject.read_tool_file()
    statisticsObject.read_annot_file()
    statisticsObject.read_3_match_file()
    statisticsObject.read_5_match_file()
    statisticsObject.read_split_cases_file()
    statisticsObject.read_no_condition_file()
    statisticsObject.read_extras_file()
    statisticsObject.display_results()

if __name__ == '__main__':
    main()