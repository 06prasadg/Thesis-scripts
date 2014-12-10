from sys import argv
import argparse
import os.path
import re

class Scatter_data_generator(object):
    score_dict = {}
    orf_length_dict = {}
    score_match_dict = {}
    score_not_match_dict = {}
    orf_len_match_dict = {}
    orf_len_not_match_dict = {}
    
    def __init__(self, result_input, no_condition_input, match_3_input, match_5_input, mga_original_file):
        self.result_input = result_input
        self.no_condition_input = no_condition_input
        self.match_3_input = match_3_input
        self.match_5_input = match_5_input
        self.mga_original_file = mga_original_file
    
    # make a score dict of all reads and their scores to be used later    
    def parser_mga_original_file(self):
        fMga_original_file = open(self.mga_original_file,"r")
        score = ""
        lines = fMga_original_file.readlines()
        for line in lines:
            if line.startswith("# gi"):
                line = line.strip()
                readName = line.replace("# ",">")
            elif(line.startswith("gene")):
                lineArr = line.split("\t")
                first = int(lineArr[1])
                second = int(lineArr[2])
                score = lineArr[6]
                length = second - first
                length = str(length)
                # scores dictionary
                if self.score_dict.has_key(readName):
                    listItems = self.score_dict.get(readName)
                    if (isinstance(listItems,str)):
                        if not (listItems == score):
                            listItems = []
                            listItems.append(self.score_dict.get(readName))
                            listItems.append(score)
                            self.score_dict[readName] = listItems
                    else:
                        for l in listItems:
                            if not(l == score):
                               listItems.append(score)
                               self.score_dict[readName] = listItems
                else:
                    self.score_dict[readName] = score
            
                # length dictionary
                if self.orf_length_dict.has_key(readName):
                    listItems = self.orf_length_dict.get(readName)
                    if (isinstance(listItems,str)):
                        if not (listItems == length):
                            listItems = []
                            listItems.append(self.orf_length_dict.get(readName))
                            listItems.append(length)
                            self.orf_length_dict[readName] = listItems
                    else:
                        for l in listItems:
                            if not(l == length):
                               listItems.append(length)
                               self.orf_length_dict[readName] = listItems
                else:
                    self.orf_length_dict[readName] = length
     
    def parse_result_files(self):
        fResult = open(self.result_input)
        resultLines = fResult.readlines()
        fNo_condition = open(self.no_condition_input)
        noConditionLines = fNo_condition.readlines()
        f3_match = open(self.match_3_input)
        match_3_lines = f3_match.readlines()
        f5_match = open(self.match_5_input)
        match_5_lines = f5_match.readlines()
        
        score_matches = open("scores_of_exact_match.csv","w")
        #score_partial = open("Scores_partial_file.csv","w")
        score_not_match = open("scores_of_in_tool.csv","w")
        
        # opening result file
        for line in resultLines:
            if line.startswith(">gi"):
                lineArr = line.split(" ")
                indicator = lineArr[2]
                lineSplit = lineArr[0].split("/")
                key = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
                
                # exact/near match cases
                if lineArr[2] == "00" or lineArr[2] == "F-1" or lineArr[2] == "F+1" or lineArr[2] == "F-2" or lineArr[2] == "F+2":
                    if self.score_dict.has_key(key):
                        self.score_match_dict[key] = self.score_dict.get(key)
                    if self.orf_length_dict.has_key(key):
                        self.orf_len_match_dict[key] = self.orf_length_dict.get(key)
                # partial matches
                #if lineArr[2] == "01" or lineArr[2] == "10" or lineArr[2] == "11" or lineArr[2] == "20" or lineArr[2] == "02" or lineArr[2] == "22" or lineArr[2] == "12" or lineArr[2] == "21":
                #    if self.score_dict.has_key(key):
                #        listItems = self.score_dict.get(key)
                #        if (isinstance(listItems,str)):
                #            score_partial.write(key + "," + str(listItems) + "\n")
                #        else:
                #            for l in listItems:
                #                score_partial.write(key + "," + str(l) + "\n")
                # not matches
                if lineArr[2] == "TT":
                    if self.score_dict.has_key(key):
                        self.score_not_match_dict[key] = self.score_dict.get(key)             
                    if self.orf_length_dict.has_key(key):
                        self.orf_len_not_match_dict[key] = self.orf_length_dict.get(key)
                
        for line in match_3_lines:                
                    if line.startswith(">gi"):
                        lineArr = line.split(" ")
                        lineSplit = lineArr[0].split("/")
                        key = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
                        if self.score_dict.has_key(key):
                            self.score_match_dict[key] = self.score_dict.get(key)
                        if self.orf_length_dict.has_key(key):
                           self.orf_len_match_dict[key] = self.orf_length_dict.get(key)
    
        for line in match_5_lines:
                    if line.startswith(">gi"):
                        lineArr = line.split(" ")
                        lineSplit = lineArr[0].split("/")
                        key = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
                        if self.score_dict.has_key(key):
                            self.score_match_dict[key] = self.score_dict.get(key)
                        if self.orf_length_dict.has_key(key):
                            self.orf_len_match_dict[key] = self.orf_length_dict.get(key)
                                    
        for line in noConditionLines:
            if line.startswith(">gi"):
                lineArr = line.split(" ")
                lineSplit = lineArr[0].split("/")
                key = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
                if self.score_dict.has_key(key):
                    self.score_not_match_dict[key] = self.score_dict.get(key)
                if self.orf_length_dict.has_key(key):
                        self.orf_len_not_match_dict[key] = self.orf_length_dict.get(key)
        
        for key, value in self.score_match_dict.items():
            value = self.score_match_dict.get(key)
            length = self.orf_len_match_dict.get(key)
            if (isinstance(value,str)):
                score_matches.write(key + "," + length + "," + value + "\n")
            else:
                for v, l in zip(value, length):
                    score_matches.write(key + "," + l + "," + v + "\n")
                    
        for key, value in self.score_not_match_dict.items():
            value = self.score_not_match_dict.get(key)
            length = self.orf_len_not_match_dict.get(key)
            if (isinstance(value,str)):
                score_not_match.write(key + "," + length + "," + value + "\n")
            else:
                for v, l in zip(value, length):
                    score_not_match.write(key + "," + l + "," + v + "\n")
        
        score_matches.close()
        score_not_match.close()

def main():
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 5 input files \n 1. Results file \n 2. No condition file \n 3. 3 Matched file \n 4. 5 Matched file \n 5. Orginal MGA file \n", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =5 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )

    args = parser.parse_args()
    result_input = args.filename[0]
    no_condition_input = args.filename[1]
    match_3_input = args.filename[2]
    match_5_input = args.filename[3]
    mga_original_file = args.filename[4]
    sdg = Scatter_data_generator(result_input, no_condition_input, match_3_input, match_5_input, mga_original_file)
    sdg.parser_mga_original_file()
    sdg.parse_result_files()
    

if __name__ == '__main__':
    main()
