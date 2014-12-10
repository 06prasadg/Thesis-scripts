"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
"""

import re
import argparse
import os.path
import sys
from datetime import datetime
from operator import itemgetter

class final_parser(object):
    toolName = ""
    rFile = None
    rExtraFile = None
    ext_3 = None
    ext_5 = None
    rExtraAnnotFile = None
    rSplitCaseFile = None
    
    split_dictionary = {}
    Annotation_dictionary = {}
    tool_dictionary = {}
    total_matched_dictionary = {}
    final_dictionary = {}
    extras_dictionary = {}
    background_dictionary = {}
    full_annot = {}
    actual_annot = {}
    diff_annot = {}
    
    def __init__(self, toolFile, fileAnnotation, path):
        self.toolFile = toolFile
        self.fileAnnotation = fileAnnotation
        self.path = path
        

    def getting_tool_name(self):
        
        toolFileOutput = open(self.toolFile,'r')
        line = toolFileOutput.readline().strip()       
        lineArr = line.split(" ")
        toolArr = lineArr[1].split("=")
        self.toolName = toolArr[1]

    # converting annotation file to dictionary
    def converting_annotation_to_dictionary(self):
        fAnnotation= open(self.fileAnnotation, 'r')
        #totalAnnotation = 0
        lines = fAnnotation.readlines()
        for line in lines:
            lineArr = line.split(" ")
            lineSplit = lineArr[0].split("/")
            key = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
            matchObject = re.findall(r'_+[0-9]*', lineSplit[1])
            value = ''.join(matchObject[0:2])  
            if self.Annotation_dictionary.has_key(key):
                listItems = self.Annotation_dictionary.get(key)
                if (isinstance(listItems,str)):
                    listItems = []
                    listItems.append(self.Annotation_dictionary.get(key))
                    listItems.append(value)
                    self.Annotation_dictionary[key] = listItems
                else:
                    listItems.append(value)
                    self.Annotation_dictionary[key] = listItems
            else:
                self.Annotation_dictionary[key] = value
            #totalAnnotations = totalAnnotations + 1
        #print "total lines in annotations ",totalAnnotations


    def converting_toolOutput_to_dictionary(self):
        toolFileOutput = open(self.toolFile,'r')
        lines = toolFileOutput.readlines()
        #tool_count = 0
        for line in lines:
            if(line.startswith(">gi")):
                lineArr = line.split(" ")
                lineSplit = lineArr[0].split("/")
                key = lineSplit[0]+"/"+''.join(list(lineSplit[1])[0])
                matchObject = re.findall(r'_+[0-9]*', lineSplit[1])
                value = ''.join(matchObject[0:2])
                if self.tool_dictionary.has_key(key):
                    listItems = self.tool_dictionary.get(key)
                    if (isinstance(listItems,str)):
                        listItems = []
                        listItems.append(self.tool_dictionary.get(key))
                        listItems.append(value)
                        self.tool_dictionary[key] = listItems
                    else:
                        listItems.append(value)
                        self.tool_dictionary[key] = listItems
                else:
                    self.tool_dictionary[key] = value
               # tool_count += 1
        #print "total lines in tool ",tool_count

    def in_annotation_not_tool(self):
        #In_annotation_but_not_in_tool=0
        #total_matched_1 = 0
        for keyAnnotation, valueAnnotation in self.Annotation_dictionary.items():
            
            ## this loop will store all readnames with coords of annot which matched with tool readnames
            ## but there will be some coords which wouldnt match with any tool coords.
            ## we need to find them
            if self.tool_dictionary.has_key(keyAnnotation):
                if (isinstance(valueAnnotation,list)):
                        #In_annotation_but_not_in_tool += len(valueAnnotation)
                        for i in valueAnnotation:
                          self.full_annot[keyAnnotation+i] = 1
                else:
                        self.full_annot[keyAnnotation+valueAnnotation]=1
            
            ## prints to final dictionary the ones whose annotation read names not matched at all
            else:
                if (isinstance(valueAnnotation,list)):
                        #In_annotation_but_not_in_tool += len(valueAnnotation)
                        for i in valueAnnotation:
                          self.final_dictionary[keyAnnotation+i+" Annotation"] = "AA "
                else:
                        #In_annotation_but_not_in_tool += 1
                        self.final_dictionary[keyAnnotation+valueAnnotation+" Annotation"] = "AA "
        
    
    def in_tool_not_annotation(self):
        #In_tool_but_not_in_annotation = 0
        #total_matched = 0
        for keyTool, valueTool in self.tool_dictionary.items():
            if self.Annotation_dictionary.has_key(keyTool):
                #if (isinstance(valueTool,list)):
                #        total_matched += len(valueTool)
                #else:
                #        total_matched += 1
                # this is our main dictionary for all the comparisons
                self.total_matched_dictionary[keyTool] = valueTool
            
            ## prints to final dictionary the ones whose tool read names not matched at all
            else:
                if (isinstance(valueTool,list)):
                        #In_tool_but_not_in_annotation += len(valueTool)
                        for i in valueTool:
                          self.final_dictionary[keyTool+i+" Tool="+self.toolName] = "TT "
                else:
                        #In_tool_but_not_in_annotation += 1
                        self.final_dictionary[keyTool+valueTool+" Tool="+self.toolName] = "TT "   

    def findMultiples(self,low,high):
        longarray = []
        for i in range(low,high+1,3):
            longarray.append(i)    
        return longarray

    def findMultiplesReverse(self,high,low):
        shortarray = []
        for i in range(high,low-1,-3):
            shortarray.append(i)
        return shortarray        

    def main_comparison(self):
        for keyTool, valueTool in self.total_matched_dictionary.items():
            valueAnnotation = self.Annotation_dictionary.get(keyTool)
            # we pass the read name, and list of tool n annot co-ords
            #finding_split_cases(keyTool, valueAnnotation, valueTool)
            if (isinstance(valueTool,list)):
                #count += len(valueTool)
                for t in valueTool:
                    if(isinstance(valueAnnotation,list)):
                        for a in valueAnnotation:
                            self.find_case_matched(keyTool,a,t)
                    else:
                        self.find_case_matched(keyTool,valueAnnotation,t)
            else:
                #count += 1
                if(isinstance(valueAnnotation,list)):
                        for a in valueAnnotation:
                            self.find_case_matched(keyTool,a,valueTool)
                else:
                        self.find_case_matched(keyTool,valueAnnotation,valueTool)
                
 
    def find_case_matched(self, key_of_tool, value_of_annotation, value_of_tool):
        flag1, flag2 = (False,False)
        
        # case 1 - exact matched (00)
        if value_of_annotation == value_of_tool:
            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "00 " +  value_of_annotation
            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
            #actual_annot[key_of_tool+value_of_annotation]=1
            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
    
        ## rest of partial matched cases
        else:       
            #partial_matched += 1
            value_of_annotation_array = value_of_annotation.split("_")
            if(int(value_of_annotation_array[1]) < int(value_of_annotation_array[2])):
               #print "innnn"
               temp = value_of_annotation_array[1]
               value_of_annotation_array[1] = value_of_annotation_array[2]
               value_of_annotation_array[2] = temp
               flag1 = True
               
            value_of_tool_array = value_of_tool.split("_")
            if(int(value_of_tool_array[1]) < int(value_of_tool_array[2])):
               temp = value_of_tool_array[1]
               value_of_tool_array[1] = value_of_tool_array[2]
               value_of_tool_array[2] = temp
               flag2 = True
            self.calculations(key_of_tool, value_of_tool, value_of_annotation, value_of_tool_array, value_of_annotation_array, flag1, flag2)

    def calculations(self, key_of_tool, value_of_tool, value_of_annotation, value_of_tool_array, value_of_annotation_array, flag1, flag2):
                
        #arrays for other 4 cases
        array1 = self.findMultiplesReverse(int(value_of_annotation_array[1]),int(value_of_annotation_array[2])) #1
        array2 = self.findMultiples(int(value_of_annotation_array[2]),int(value_of_annotation_array[1])) #2
        array3 = self.findMultiples(int(value_of_annotation_array[1]),int(value_of_tool_array[1])) #3
        array4 = self.findMultiplesReverse(int(value_of_annotation_array[2]),0) #4
        
        if(not self.background_dictionary.has_key(key_of_tool+value_of_tool+" Tool="+self.toolName)):
            if (flag1 and flag2):
    
                if (value_of_tool_array[1]==value_of_annotation_array[1]):
                    backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[2]),0)   
                    forwardArray = self.findMultiples(int(value_of_annotation_array[2]),int(value_of_annotation_array[1]))
                    for b in backwardArray:
                        if b==int(value_of_tool_array[2]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "10 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            #actual_annot[key_of_tool+value_of_annotation]=1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "3'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                    for f in forwardArray:
                        if f == int(value_of_tool_array[2]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "20 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                           
                        elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "3'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                
                ##20, 10
                elif (value_of_tool_array[2]==value_of_annotation_array[2]):  
                    backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[1]),int(value_of_annotation_array[2]))
                    forwardArray = self.findMultiples(int(value_of_annotation_array[1]),int(value_of_tool_array[1]))
                    for b in backwardArray:
                        if b==int(value_of_tool_array[1]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "02 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "5'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                        
                    for f in forwardArray:        
                        if f == int(value_of_tool_array[1]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "01 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "5'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                 
                ## 22, 11, 21, 12           
                elif(self.array_results(array1,array2,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "22 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    
                elif(self.array_results(array3,array4,value_of_tool_array)):  
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "11 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)  
                        
                elif(self.array_results(array1,array4,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "12 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    
                elif(self.array_results(array3,array2,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "21 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                
                ##frameshift loops
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+1) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+1)):
                        self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F-1 " +  value_of_annotation
                        self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                        self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                               
        
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-1) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-1)):
                        self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F+1 " +  value_of_annotation
                        self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                        self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
        
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+2) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+2)):
                        self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F-2 " +  value_of_annotation
                        self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                        self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
        
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-2) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-2)):
                        self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F+2 " +  value_of_annotation
                        self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                        self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    
                else:
                    pass
            
            
            elif((not flag1) and (not flag2)):
    
                if (value_of_tool_array[1]==value_of_annotation_array[1]):
                    
                    backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[2]),0)   
                    forwardArray = self.findMultiples(int(value_of_annotation_array[2]),int(value_of_annotation_array[1]))
                    for b in backwardArray:          
                        if b==int(value_of_tool_array[2]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "01 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "5'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                            
                    for f in forwardArray:
                        if f == int(value_of_tool_array[2]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "02 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                            (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "5'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                
                ##20, 10
                elif (value_of_tool_array[2]==value_of_annotation_array[2]):
                    backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[1]),int(value_of_annotation_array[2]))
                    forwardArray = self.findMultiples(int(value_of_annotation_array[1]),int(value_of_tool_array[1]))
                    for b in backwardArray:
                        if b==int(value_of_tool_array[1]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "20 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
    
                        elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "3'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                    for f in forwardArray:        
                        if f == int(value_of_tool_array[1]):
                            self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "10 " +  value_of_annotation
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
    
                        elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                            (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.extras_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "3'_Matched " +  value_of_annotation
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                        else:
                            pass
                            
                ## 22, 11, 21, 12           
                elif(self.array_results(array1,array2,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "22 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
    
                elif(self.array_results(array3,array4,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "11 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)  
                        
                elif(self.array_results(array1,array4,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "21 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    
                elif(self.array_results(array3,array2,value_of_tool_array)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "12 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
        
                ##frameshift loops
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+1) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+1)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F-1 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
        
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-1) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-1)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F+1 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
        
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+2) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+2)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F-2 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
        
                elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-2) and
                    (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-2)):
                    self.final_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = "F+2 " +  value_of_annotation
                    self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                    self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
    
                else:
                    pass
             
            else:
                pass
    
    def making_dict(self, key_of_tool, value_of_annotation, value_of_tool):
    
        if self.actual_annot.has_key(key_of_tool+value_of_annotation):
            listItems = self.actual_annot.get(key_of_tool+value_of_annotation)
            if (isinstance(listItems,str)):
                listItems = []
                listItems.append(self.actual_annot.get(key_of_tool+value_of_annotation))
                listItems.append(value_of_tool)
                self.actual_annot[key_of_tool+value_of_annotation] = listItems
            else:
                        listItems.append(value_of_tool)
                        self.actual_annot[key_of_tool+value_of_annotation] = listItems
        else:
            self.actual_annot[key_of_tool+value_of_annotation] = value_of_tool
    
    def array_results(self, array_1,array_2,value_of_tool_array):
        for a1,a2 in zip(array_1,array_2):
            if (a1 == int(value_of_tool_array[1]) and a2 == int(value_of_tool_array[2])):
                return True
            else:
                return False
    
    def print_final_results(self):

        rFile = os.path.join(self.path,"results_file_" + str(self.toolName) + ".tab")
        ext_3 = os.path.join(self.path,"3_matched_5!_diff_" + str(self.toolName) + ".tab")
        ext_5 = os.path.join(self.path,"5_matched_3!_diff_" + str(self.toolName) + ".tab")
        results_file = open(rFile, 'w')
        extra_3_file = open(ext_3, 'w')
        extra_5_file = open(ext_5, 'w')
        
        for key, value in self.final_dictionary.items():
            results_file.write(key + " " + value + "\n")
        results_file.close()

        for key, value in self.extras_dictionary.items():
            valueArr = value.split(" ")
            if(valueArr[0] == "3'_Matched"):
                extra_3_file.write(key + " " + value + "\n")
            else:
                extra_5_file.write(key + " " + value + "\n")
     
    def finding_tool_extras(self):
        #count = 0
        rExtraFile = os.path.join(self.path,"NoCondition_satisfied_mismatch" + str(self.toolName) + ".tab")
        extra_file = open(rExtraFile,'w')  
        for key, value in self.total_matched_dictionary.items():
            valueAnnots = self.Annotation_dictionary.get(key)
            ## not background dict gives all values of read+coords of tool not matched
            if (isinstance(value,list)):
                for i in value:
                    if not(self.background_dictionary.has_key(key+i+" Tool="+self.toolName)):
                       #count +=1
                       extra_file.write(key+i+" Tool="+self.toolName+ " "+ str(valueAnnots) + "\n")
            else: 
                if not(self.background_dictionary.has_key(key+value+" Tool="+self.toolName)):
                    #count += 1
                    extra_file.write(key+value+" Tool="+self.toolName+ " "+ str(valueAnnots) + "\n")
    
        #print "Mismatched + No condition of Tool count : ",count
            
    def annot_extras(self):
        rExtraAnnotFile = os.path.join(self.path,"Extra_Annots_" + str(self.toolName) + ".tab")
        extra_annot_file = open(rExtraAnnotFile,'w')
        print "len of actual annot : ", len(self.actual_annot)
        annot_extra_count = 0
            
        #annot_not_match_count = 0
        for key in self.full_annot.keys():
            if(self.actual_annot.has_key(key)):
                pass
            else:
                self.diff_annot[key] = 1
                
        for key in self.diff_annot.keys():
            extra_annot_file.write(key + "\n")
                  
        #print "Annotation coords not matched count : ", len(diff_annot)
      
    def main_comparison_split_case(self):
        for keyTool, valueTool in self.total_matched_dictionary.items():
            valueAnnotation = self.Annotation_dictionary.get(keyTool)
            # we pass the read name, and list of tool n annot co-ords
            self.finding_split_cases(keyTool, valueAnnotation, valueTool)

    def finding_split_cases(self, keyTool,valueAnnotation, valueTool):
        flag1, flag2 = (False,False)
        temp = 0
        tool_coords_array = []
        final_value_of_tool_array = []
        
        ## Annotation co-ord is 1
        if(isinstance(valueAnnotation,str)):
            ## tool predicted > 1
            if(isinstance(valueTool,list)):
                value_of_annotation_array = valueAnnotation.split("_")
                
                #exchange nos if smaller first (keep a standard order of larger first then smaller for both tool, annoation co-ords)
                ## flags true means smaller first, larger later
                if(int(value_of_annotation_array[1]) < int(value_of_annotation_array[2])):
                    temp = int(value_of_annotation_array[1])
                    value_of_annotation_array[1] = int(value_of_annotation_array[2])
                    value_of_annotation_array[2] = temp
                    flag1 = True
                # just to convert string to int
                else:
                    value_of_annotation_array[1] = int(value_of_annotation_array[1])
                    value_of_annotation_array[2] = int(value_of_annotation_array[2])
                    
                for value_of_tool in valueTool:
                    value_of_tool_array = value_of_tool.split("_")
                    if(int(value_of_tool_array[1]) < int(value_of_tool_array[2])):
                        temp = int(value_of_tool_array[1])
                        value_of_tool_array[1] = int(value_of_tool_array[2])
                        value_of_tool_array[2] = temp
                        flag2 = True
                    # just to convert string to int
                    else:
                        value_of_tool_array[1] = int(value_of_tool_array[1])
                        value_of_tool_array[2] = int(value_of_tool_array[2])
                     
                    # convert list of strings to list of int    
                    tool_coords_array.append(map(int,value_of_tool_array[1:]))
                #print tool_coords_array, " "    
                tool_coords_array = sorted(tool_coords_array, key=itemgetter(1))
                #print tool_coords_array, " "
                minimum = tool_coords_array[0][1]
                maximum = tool_coords_array[-1][0]
                #print "min : ", minimum, " ", "max : ", maximum
                # final tool array has no of co-ords in tool for read, minimum co-ord and max co-ord
                final_value_of_tool_array = [len(tool_coords_array), int(maximum), int(minimum)]
               ## passing read name, final tool array, tool co-ords for read, annotation co-ords for read
                self.calculations_split_case(keyTool, value_of_annotation_array, final_value_of_tool_array,
                tool_coords_array, valueAnnotation, flag1, flag2, valueTool)
         
            ## tool predicted 1 so do nothing - its just normal case
            else:
                       pass
            
        ## Annotation co-ord > 1        
        else:
            if(isinstance(valueTool,list)):
                for value_of_annotation in valueAnnotation:
                    value_of_annotation_array = value_of_annotation.split("_")
                    if(int(value_of_annotation_array[1]) < int(value_of_annotation_array[2])):
                       temp = int(value_of_annotation_array[1])
                       value_of_annotation_array[1] = int(value_of_annotation_array[2])
                       value_of_annotation_array[2] = temp
                       flag1 = True
                    # just to convert string to int
                    else:
                        value_of_annotation_array[1] = int(value_of_annotation_array[1])
                        value_of_annotation_array[2] = int(value_of_annotation_array[2])
                    
                for value_of_tool in valueTool: 
                    value_of_tool_array = value_of_tool.split("_")
                    if(int(value_of_tool_array[1]) < int(value_of_tool_array[2])):
                        temp = int(value_of_tool_array[1])
                        value_of_tool_array[1] = int(value_of_tool_array[2])
                        value_of_tool_array[2] = temp
                        flag2 = True
                    # just to convert string to int
                    else:
                        value_of_tool_array[1] = int(value_of_tool_array[1])
                        value_of_tool_array[2] = int(value_of_tool_array[2])
                     
                     # convert list of strings to list of int       
                    tool_coords_array.append(map(int,value_of_tool_array[1:]))
                #print tool_coords_array, " "     
                tool_coords_array = sorted(tool_coords_array, key=itemgetter(1))
                #print tool_coords_array, " "
                minimum = tool_coords_array[0][1]
                maximum = tool_coords_array[-1][0]
                #print "min : ", minimum, " ", "max : ", maximum
                # final tool array has no of co-ords in tool for read, minimum co-ord and max co-ord
                final_value_of_tool_array = [len(tool_coords_array), int(maximum), int(minimum)]
                ## passing read name, final tool array, tool co-ords for read, annotation co-ords for read
                self.calculations_split_case(keyTool, value_of_annotation_array, final_value_of_tool_array, tool_coords_array,
                value_of_annotation, flag1, flag2, valueTool)
                
            ## tool predicted 1 so do nothing - its just normal case
            else:
                       pass
        
    def calculations_split_case(self, key_of_tool, value_of_annotation_array,
                    value_of_tool_array, tool_coords_array, value_of_annotation, flag1, flag2, valueTool):
                
        #arrays for other 4 cases
        array1 = self.findMultiplesReverse(int(value_of_annotation_array[1]),int(value_of_annotation_array[2])) #1
        array2 = self.findMultiples(int(value_of_annotation_array[2]),int(value_of_annotation_array[1])) #2
        array3 = self.findMultiples(int(value_of_annotation_array[1]),int(value_of_tool_array[1])) #3
        array4 = self.findMultiplesReverse(int(value_of_annotation_array[2]),0) #
        tool_coords_num = str(value_of_tool_array[0])
        
        if (flag1 and flag2):
    
            if(value_of_tool_array[1]==value_of_annotation_array[1] and value_of_tool_array[2]==value_of_annotation_array[2]):
               self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 00"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
               
               for t in tool_coords_array:
                   value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                   self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                   self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
               
            elif (value_of_tool_array[1]==value_of_annotation_array[1]):
                backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[2]),0)   
                forwardArray = self.findMultiples(int(value_of_annotation_array[2]),int(value_of_annotation_array[1]))
                for b in backwardArray:
                    if b==int(value_of_tool_array[2]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 10"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
                    elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 3'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
                    
                for f in forwardArray:
                    if f == int(value_of_tool_array[2]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 20"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
                    elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 3'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
            
            ##20, 10
            elif (value_of_tool_array[2]==value_of_annotation_array[2]):  
                backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[1]),int(value_of_annotation_array[2]))
                forwardArray = self.findMultiples(int(value_of_annotation_array[1]),int(value_of_tool_array[1]))
                for b in backwardArray:
                    if b==int(value_of_tool_array[1]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 02"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + str(value_of_annotation_array[1]) + ":" + str(value_of_annotation_array[2])
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 5'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + str(value_of_annotation_array[1]) + ":" + str(value_of_annotation_array[2])
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
                    
                for f in forwardArray:        
                    if f == int(value_of_tool_array[1]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 01"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 5'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
             
            ## 22, 11, 21, 12           
            elif(self.array_results(array1,array2,value_of_tool_array)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 22"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif(self.array_results(array3,array4,value_of_tool_array)):  
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 11"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    
            elif(self.array_results(array1,array4,value_of_tool_array)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 12"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[2])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif(self.array_results(array3,array2,value_of_tool_array)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 21"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
            ##frameshift loops
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+1) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+1)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F-1"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)           
        
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-1) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-1)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F+1"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+2) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+2)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F-2"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-2) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-2)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F+2"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[1])+"_"+str(t[0])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            else:
                        pass
    
        elif(not flag1 and not flag2):
    
            if(value_of_tool_array[1]==value_of_annotation_array[1] and value_of_tool_array[2]==value_of_annotation_array[2]):
               self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 00"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
               
               for t in tool_coords_array:
                   value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                   self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                   self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
               
            elif (value_of_tool_array[1]==value_of_annotation_array[1]):
                backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[2]),0)   
                forwardArray = self.findMultiples(int(value_of_annotation_array[2]),int(value_of_annotation_array[1]))
                for b in backwardArray:
                    if b==int(value_of_tool_array[2]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 01"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
                    elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 5'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
                    
                for f in forwardArray:
                    if f == int(value_of_tool_array[2]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 02"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
                    elif ( (int(value_of_annotation_array[2])+2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-2 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])+1 == int(value_of_tool_array[2])) or
                        (int(value_of_annotation_array[2])-1 == int(value_of_tool_array[2])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 5'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
            
            ##20, 10
            elif (value_of_tool_array[2]==value_of_annotation_array[2]):  
                backwardArray = self.findMultiplesReverse(int(value_of_annotation_array[1]),int(value_of_annotation_array[2]))
                forwardArray = self.findMultiples(int(value_of_annotation_array[1]),int(value_of_tool_array[1]))
                for b in backwardArray:
                    if b==int(value_of_tool_array[1]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 20"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 3'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
                    
                for f in forwardArray:        
                    if f == int(value_of_tool_array[1]):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 10"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    elif ( (int(value_of_annotation_array[1])+2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-2 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])+1 == int(value_of_tool_array[1])) or
                        (int(value_of_annotation_array[1])-1 == int(value_of_tool_array[1])) ):
                        self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 3'_Matched"] \
                        = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                        for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    else:
                        pass
             
            ## 22, 11, 21, 12           
            elif(self.array_results(array1,array2,value_of_tool_array)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 22"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif(self.array_results(array3,array4,value_of_tool_array)):  
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 11"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                    
            elif(self.array_results(array1,array4,value_of_tool_array)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 21"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif(self.array_results(array3,array2,value_of_tool_array)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " 12"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
            ##frameshift loops
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+1) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+1)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F+1"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)           
        
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-1) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-1)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F-1"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])+2) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])+2)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F+2"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            elif((int(value_of_annotation_array[1]) == int(value_of_tool_array[1])-2) and
                (int(value_of_annotation_array[2]) == int(value_of_tool_array[2])-2)):
                self.split_dictionary[key_of_tool + "_" + str(value_of_tool_array[1]) + "_" + str(value_of_tool_array[2]) + " F-2"] \
                = " " + "split_case_" + tool_coords_num + " " + str(valueTool) + " " + value_of_annotation
                for t in tool_coords_array:
                            value_of_tool = "_"+str(t[0])+"_"+str(t[1])
                            self.background_dictionary[key_of_tool+value_of_tool+" Tool="+self.toolName] = 1
                            self.making_dict(key_of_tool, value_of_annotation, value_of_tool)
                            
            else:
                        pass
        else:
            pass
        
    def print_split_cases(self):
        rSplitCaseFile = os.path.join(self.path,"split_case_file_" + str(self.toolName) + ".tab")
        split_cases_file = open(rSplitCaseFile, 'w')
        for key, value in self.split_dictionary.items():
            split_cases_file.write(key + value + "\n")
        split_cases_file.close()
    
def main():
    
    print "Time script started :",datetime.now().hour,":",datetime.now().minute,":",datetime.now().second
    parser = argparse.ArgumentParser(description = "User Manual \n -i = input_file \n This script requires 2 input file \n 1. Co-ordniates file generated from Annotation file \n 2. Output Generated in standard output format for tools any of these tools :- MetaGeneMark  OR MetaGeneAnnotator OR Orphelia or Glimmer-MG", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i" , dest = "filename", help = "Input_file" , metavar = "file" , nargs =2 , required = True)
    parser.add_argument("-o", dest = "output_directory", help = "output_directory" , metavar = "path" )
    
    args = parser.parse_args()
    fileAnnotation = args.filename[0]
    toolFile = args.filename[1]
    
    output_directory = args.output_directory
    if output_directory is None:
        path=""
    else:
        path = output_directory
    fpObject = final_parser(toolFile, fileAnnotation, path)
    fpObject.getting_tool_name()
    
    fpObject.converting_annotation_to_dictionary()
    fpObject.converting_toolOutput_to_dictionary()
    fpObject.in_annotation_not_tool()
    fpObject.in_tool_not_annotation()
    fpObject.main_comparison()
    fpObject.main_comparison_split_case()
    fpObject.finding_tool_extras()
    fpObject.annot_extras()
    fpObject.print_final_results()
    fpObject.print_split_cases()
    
    print "Time script ended :",datetime.now().hour,":",datetime.now().minute,":",datetime.now().second  

if __name__ == '__main__':
    main()
