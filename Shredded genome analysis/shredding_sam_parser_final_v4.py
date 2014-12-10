"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Script to parse the annotation file predictions and put them in a standard common format
"""

import re
from sys import argv
from datetime import datetime

class shredding_parser(object):
    
    Annotation_dictionary = {}
    sam_dict = {}
    final_dict = {}
    
    def __init__(self, shreddedfile, fileAnnotation):
        self.shreddedfile = shreddedfile
        self.fileAnnotation = fileAnnotation
          
    # converting annotation file to dictionary
    def converting_annotation_to_dictionary(self):
        fAnnotation = open(self.fileAnnotation,'r')
        lines = fAnnotation.readlines()
        for line in lines:
            flag = ""
            matchObject = re.findall(r'_+[0-9]*', line)
            #value = ''.join(matchObject[2:])
            #print matchObject
            p = list(matchObject[2])
            ORF_start_coordNum = int(''.join(p[1:]))
            #print "ORF start ", ORF_start_coordNum
            q = list(matchObject[3])
            ORF_end_coordNum = int(''.join(q[1:]))
            #print "ORF end ", ORF_end_coordNum
            
            if(ORF_start_coordNum<ORF_end_coordNum):
               flag = ":smaller_first"
               self.Annotation_dictionary[ORF_start_coordNum] = str(ORF_end_coordNum)+flag  #set a flag to 1 say that start < end eg. 345:890
            else:
               flag = ":larger_first"
               self.Annotation_dictionary[ORF_end_coordNum] = str(ORF_start_coordNum)+flag
            #print ORF_start_coordNum, ORF_end_coordNum, flag
        
    def shreddedFiletoDict(self):
        fshredded = open(self.shreddedfile, 'r')
        lineLen = 0
        mFactor = ""
        lines = fshredded.readlines()
        for line in lines:
            if line.startswith("gi"):
                lineArr = line.split("\t")
                #print lineArr
                indicator = lineArr[1]
                firstNum = int(lineArr[3])
                mFactor = lineArr[5]
                seqlen = len(str(lineArr[9]))
                ##since total legnth is 300 so if say 11+300 = 311, hence seq from 11-310
                secondNum = firstNum+seqlen-1
                numArray = [firstNum,secondNum,mFactor,seqlen,indicator]
                self.sam_dict[lineArr[0]] = numArray
        #print sam_dict

    def comparison(self):
        for keySAM, valueSAM in self.sam_dict.items():
            SAM_start_coord = int(valueSAM[0])
            SAM_end_coord = int(valueSAM[1])
            cigar_factor = valueSAM[2]
            seqlen = int(valueSAM[3])
            indicator = int(valueSAM[4])
            readName = keySAM
            insertions = 0
            deletions  = 0
            
            if cigar_factor == "300M":
                pass
            
            # if cigar factor anything else other than 300M, do all this.
            else:
                numList = []
                letterList = []
                matchObjNum = re.findall(r'[0-9]*',cigar_factor)
                for m in matchObjNum:
                    if m != '':
                        numList.append(m)
                matchObjLetter = re.findall(r'[A-Z]*',cigar_factor)
                for m in matchObjLetter:
                    if m != '':
                       letterList.append(m)
                # find no of insertions
                for n,l in zip(numList,letterList):
                    if l == "I":
                       insertions = int(n)
                    if l == "D":
                       deletions = int(n)
            
            for keyAnnot, valueAnnot in self.Annotation_dictionary.items():
                ORF_start_coord = int(keyAnnot)
                ORF_end_coord = int(valueAnnot.split(":")[0])
                flag = valueAnnot.split(":")[1]
                self.find_case(SAM_start_coord, SAM_end_coord, cigar_factor, readName, seqlen, insertions, deletions, ORF_start_coord,
                          ORF_end_coord, flag, indicator)


    def find_case(self, SAM_start_coord, SAM_end_coord, cigar_factor, readName, seqlen, insertions, deletions,
                  ORF_start_coord, ORF_end_coord, flag, indicator):
        #print "Case matching"
        final_dict_key = ""
        final_dict_value = ""
        SAM_start_coord_new = 0
        SAM_end_coord_new = 0
        SAM_start_end_coord_new = ""
        SAM_old_coord = ""
        ORF_start_end_coord = ""
        
        ## case 1
        ## eg: ORF - 2189_12
        ## SAM - 60_359
        if( ((SAM_start_coord > ORF_start_coord) or (SAM_start_coord == SAM_end_coord)) and ( (SAM_end_coord < ORF_end_coord) or (SAM_end_coord == ORF_end_coord)) and (SAM_end_coord > ORF_start_coord)):
            #print seqlen
            SAM_start_coord_new = 1
            SAM_end_coord_new = seqlen
            ##making changes 9/14/2014
            if indicator == 0:
                if flag == "larger_first":
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            
            ## if indicator is 16
            else:
                if flag == "larger_first":
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            
            SAM_old_coord = "_" + str(SAM_start_coord) + "_" + str(SAM_end_coord)
            final_dict_key = readName + "_" + SAM_start_end_coord_new
            final_dict_value = ORF_start_end_coord + " " + cigar_factor + " " + SAM_old_coord
            self.final_dict[final_dict_key] = final_dict_value
        
        ## case2
        ## ORF - 2189_12
        ## SAM 3_302
        elif((SAM_start_coord < ORF_start_coord) and ((SAM_end_coord < ORF_end_coord) or (SAM_end_coord == ORF_end_coord)) and (SAM_end_coord > ORF_start_coord)):
            
            ## making changes 9/12/2014
            if indicator == 0:
                SAM_start_coord_new = ORF_start_coord - SAM_start_coord + 1
                SAM_end_coord_new = seqlen
                if flag == "larger_first":
                    # 300_10
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    # 10_300
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            
            ## if indicator is 16
            else:
                SAM_start_coord_new = SAM_end_coord - ORF_start_coord + 1
                SAM_end_coord_new = 1
                if flag == "larger_first":
                    # 1_290
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    # 290_1       
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            
            SAM_old_coord = "_" + str(SAM_start_coord) + "_" + str(SAM_end_coord)
            final_dict_key = readName + "_" + SAM_start_end_coord_new
            final_dict_value = ORF_start_end_coord + " " + cigar_factor + " " + SAM_old_coord
            self.final_dict[final_dict_key] = final_dict_value
         
        ## case 3
        ## ORF 2189_12
        ## SAM 2160_2459, 2170_2469
        elif(((SAM_start_coord > ORF_start_coord) or (SAM_start_coord == ORF_start_coord)) and  (SAM_start_coord < ORF_end_coord) and (SAM_end_coord > ORF_end_coord)):
            
            
            ##making changes 9/14/2014
            if indicator == 0:
                SAM_start_coord_new = 1
                SAM_end_coord_new = ORF_end_coord - SAM_start_coord + 1
                if flag == "larger_first":
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            else:
                ## making changes on 9/15/2014
                SAM_start_coord_new = SAM_end_coord - ORF_end_coord + 1
                SAM_end_coord_new = seqlen
                if flag == "larger_first":
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            
            SAM_old_coord = "_" + str(SAM_start_coord) + "_" + str(SAM_end_coord)
            final_dict_key = readName + "_" + SAM_start_end_coord_new
            final_dict_value = ORF_start_end_coord + " " + cigar_factor + " " + SAM_old_coord
            self.final_dict[final_dict_key] = final_dict_value
        
        ## case 4
        ## ORF 2403_2200
        ## SAM 2190_2489
        elif((SAM_start_coord < ORF_start_coord) and (SAM_end_coord > ORF_end_coord)):
            
            
             ##making changes 9/14/2014
            if indicator == 0:
                SAM_start_coord_new = ORF_start_coord - SAM_start_coord + 1 # 11
                SAM_end_coord_new = ORF_end_coord - SAM_start_coord + 1 # 214
                ## NOT SURE ON THIS CASE
                if flag == "larger_first":
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                else:
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            else:
                SAM_start_coord_new = SAM_end_coord - ORF_end_coord + 1 #87
                SAM_end_coord_new = SAM_end_coord - ORF_start_coord + 1 #290 
                if flag == "larger_first":
                    SAM_start_end_coord_new = str(SAM_start_coord_new)+"_"+str(SAM_end_coord_new)
                    ORF_start_end_coord = str(ORF_end_coord) + ":" + str(ORF_start_coord)
                ## NOT SURE ON THIS CASE
                else:
                    SAM_start_end_coord_new = str(SAM_end_coord_new)+"_"+str(SAM_start_coord_new)
                    ORF_start_end_coord = str(ORF_start_coord) + ":" + str(ORF_end_coord)
            
            SAM_old_coord = "_" + str(SAM_start_coord) + "_" + str(SAM_end_coord)
            final_dict_key = readName + "_" + SAM_start_end_coord_new
            final_dict_value = ORF_start_end_coord + " " + cigar_factor + " " + SAM_old_coord
            self.final_dict[final_dict_key] = final_dict_value
        
        else:
            pass
        
    
    def write_output(self):
        fOutput = open("Output_mapped_sam.tab" , 'w')
        for k, v in self.final_dict.items():
            fileLine = ">" + k + " " + v + "\n";
            fOutput.write(fileLine)       
        fOutput.close()
        
def main():

    shreddedfile = argv[1]
    fileAnnotation = argv[2]
    
    print "Time script started :",datetime.now().hour,":",datetime.now().minute,":",datetime.now().second
    sp = shredding_parser(shreddedfile, fileAnnotation)
    sp.converting_annotation_to_dictionary()
    sp.shreddedFiletoDict()
    sp.comparison()
    sp.write_output()
    print "Time script ended :",datetime.now().hour,":",datetime.now().minute,":",datetime.now().second 

if __name__ == '__main__':
    main()
