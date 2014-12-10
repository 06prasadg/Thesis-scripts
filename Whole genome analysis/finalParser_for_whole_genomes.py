"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu

"""


import re
import argparse
import os.path
import sys

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
   
#cFile = os.path.join(path,"comparison_File_" + '.'.join(toolFile.split(".")[0:-1]) + ".tab")
#sFile = os.path.join(path,"statistics_File.tab")

fAnnotation= open(fileAnnotation, 'r')
toolFileOutput = open(toolFile,'r')

#comparison_file = open("comparison_File.tab",'w')
#statistics_file = open(sFile,'w')

Annotation_dictionary = {}
tool_dictionary = {}
notMatched_dictionary = {}
notFound_dictionary = {}
finalDictionary = {}
background_dictionary = {}
tracking_annotation_keys = {}

global totalAnnotations
global matched
global partial_matched
global frame_shift
global In_annotation_but_not_in_tool
global tool_count
global found_in_tool_but_not_in_annotation
global rFile
global results_file
global sFile
global statistics_file


def naming_result_file():
    toolFileOutput = open(toolFile,'r')
    line = toolFileOutput.readline().strip()
    #print line
    lineArr = line.split(" ")
    toolArr = lineArr[-1].split("=")
    toolName = toolArr[1]
    #print toolName
    rFile = os.path.join(path,"results_file_" + str(toolName) + ".tab")
    return rFile

def naming_statistics_file():
    toolFileOutput = open(toolFile,'r')
    line = toolFileOutput.readline().strip()
    #print line
    lineArr = line.split(" ")
    toolArr = lineArr[-1].split("=")
    toolName = toolArr[1]
    #print toolName
    sFile = os.path.join(path,"statistics_file_" + str(toolName) + ".tab")
    return sFile

# converting annotation file to dictionary
def converting_annotation_to_dictionary():
    totalAnnotations = 0
    lines = fAnnotation.readlines()
    for line in lines:
        #print line
        matchObject = re.findall(r'_+[0-9]*', line)
        #print matchObject
        value = ''.join(matchObject[2:])
        #print value
        Annotation_dictionary[value] = line.strip()
        totalAnnotations = totalAnnotations + 1
    #print Annotation_dictionary
    #print "totalAnnotations " , totalAnnotations
    return totalAnnotations

def converting_toolOutput_to_dictionary():
    lines = toolFileOutput.readlines()
    tool_count = 0
    for line in lines:
        #print line
        matchObject = re.findall(r'_+[0-9]*', line)
        #print matchObject
        value = ''.join(matchObject[1:3])
        #print value
        tool_dictionary[value] = line.strip()
        tool_count += 1
    return tool_count
    #print tool_dictionary
    #print "total: " , total
    
def main_comparison_function():
    matched = 0
    for key,value in Annotation_dictionary.items():
        #print "key: ", key
        #print "Array: ", key.split("_")
        if(tool_dictionary.has_key(key)):
            flag = 1
            matched = matched + 1
            finalDictionary[tool_dictionary.get(key)] = "00 " + key
            background_dictionary[key]=value
            tracking_annotation_keys[key]=1
        else:
            #flag = 0
            pass
    return matched

def Found_in_Annotation_but_not_in_tool():
    toolFileOutput = open(toolFile,'r')
    line = toolFileOutput.readline().strip()
    lineArr = line.split(" ")
    toolArr = lineArr[-1].split("=")
    toolName = toolArr[1]
    xFile = os.path.join(path,"XX_results_" + str(toolName) + ".tab")
    f = open(xFile,"w")
    In_annotation_but_not_in_tool = 0
    for key,value in Annotation_dictionary.items():
        #print key
        if( not(background_dictionary.has_key(key)) and not(tracking_annotation_keys.has_key(key))):
            #print value, "XX"
            f.write(value + " XX " + '\n')
            #results_file.write(value + " XX " + '\n')
            In_annotation_but_not_in_tool += 1
    return In_annotation_but_not_in_tool

def findMultiples(low,high):
    longarray = []
    for i in range(low,high+1,3):
        longarray.append(i)    
    return longarray

def findMultiplesReverse(high,low):
    shortarray = []
    for i in range(high,low-1,-3):
        shortarray.append(i)
    return shortarray
        
def to_find_nonMatching():
    
    """
Here is how I am going number/name the output:

00 = 5'3' matched
01 = 5' matched 3' long
02 = 5' matched 3' short
10 = 5' long 3' matched
20 = 5' short 3' matched
11 = 5' long 3' long
22 = 5' short 3' short
12 = 5' long 3' short
21 = 5' short 3' long
XX = found in annotation file but not in tool output
YY = found in tool output but not in annotation file

    """
    
    partial_matched = 0
    frame_shift = 0
    found_in_tool_but_not_in_annotation= 0
    
    for key,value in tool_dictionary.items():
        if not (Annotation_dictionary.has_key(key)):
            if not (notMatched_dictionary.has_key(key)):
                notMatched_dictionary[key] = value
     
    for key,value in notMatched_dictionary.items():
        keyNotMatchedArray = key.split("_")
        for key1 in Annotation_dictionary.keys():
            keyAnnotationArray = key1.split("_")
            if (keyNotMatchedArray[1]==keyAnnotationArray[1] and not (finalDictionary.has_key(value))):
                
                backwardArray = findMultiplesReverse(int(keyAnnotationArray[2]),0)
                
                forwardArray = findMultiples(int(keyAnnotationArray[2]),int(keyAnnotationArray[1]))
                
                for b in backwardArray:
                    if b==int(keyNotMatchedArray[2]):
                        #print "5' 3 Long loop"
                        finalDictionary[value]="01 " + key1
                        background_dictionary[key] = value
                        tracking_annotation_keys[key1]=1
                        partial_matched += 1
                        
                for f in forwardArray:        
                    if f == int(keyNotMatchedArray[2]):
                       #print "5' 3 Short loop"
                       finalDictionary[value]="02 " + key1
                       background_dictionary[key] = value
                       tracking_annotation_keys[key1]=1
                       partial_matched += 1
                
            elif (keyNotMatchedArray[2]==keyAnnotationArray[2] and not (finalDictionary.has_key(value))):
                
                backwardArray = findMultiplesReverse(int(keyAnnotationArray[1]),int(keyAnnotationArray[2]))
                
                forwardArray = findMultiples(int(keyAnnotationArray[1]),int(keyNotMatchedArray[1]))
                
                for b in backwardArray:
                    if b==int(keyNotMatchedArray[1]):
                        #print "5' 3 Long loop"
                        finalDictionary[value]="20 " + key1
                        background_dictionary[key] = value
                        tracking_annotation_keys[key1]=1
                        partial_matched += 1
                        
                for f in forwardArray:        
                    if f == int(keyNotMatchedArray[1]):
                       #print "5' 3 Short loop"
                       finalDictionary[value]="10 " + key1
                       background_dictionary[key] = value
                       tracking_annotation_keys[key1]=1
                       partial_matched += 1
            
            ## 
            elif (not (finalDictionary.has_key(value))):
                
                array1 = findMultiplesReverse(int(keyAnnotationArray[1]),int(keyAnnotationArray[2])) #1
                array2 = findMultiples(int(keyAnnotationArray[2]),int(keyAnnotationArray[1])) #2
                array3 = findMultiples(int(keyAnnotationArray[1]),int(keyNotMatchedArray[1])) #3
                array4 = findMultiplesReverse(int(keyAnnotationArray[2]),0) #4
                
                #other 4 cases
                for a1,a2 in zip(array1,array2):
                    if (a1 == int(keyNotMatchedArray[1]) and a2 == int(keyNotMatchedArray[2])):
                        finalDictionary[value] = "22 " + key1
                        background_dictionary[key] = value
                        tracking_annotation_keys[key1]=1
                        partial_matched += 1
                 
                for a3,a4 in zip(array3,array4):   
                    if(a3 == int(keyNotMatchedArray[1]) and a4 == int(keyNotMatchedArray[2])):
                       finalDictionary[value] = "11 " + key1
                       background_dictionary[key] = value
                       tracking_annotation_keys[key1]=1
                       partial_matched += 1
                       
                for a1,a4 in zip(array1,array4):
                    if(a1 == int(keyNotMatchedArray[1]) and a4 == int(keyNotMatchedArray[2])):
                       finalDictionary[value] = "21 " + key1
                       background_dictionary[key] = value
                       tracking_annotation_keys[key1]=1
                       partial_matched += 1
                       
                for a3,a2 in zip(array3,array2):
                    if(a3 == int(keyNotMatchedArray[1]) and a2 == int(keyNotMatchedArray[2])):
                       finalDictionary[value] = "12 " + key1
                       background_dictionary[key] = value
                       tracking_annotation_keys[key1]=1
                       partial_matched += 1
                
                ##frameshift loops
                if((int(keyAnnotationArray[1]) == int(keyNotMatchedArray[1])+1) and (int(keyAnnotationArray[2]) == int(keyNotMatchedArray[2])+1)):
                    finalDictionary[value] = " : frame shift by 1 - " + key1
                    background_dictionary[key] = value
                    tracking_annotation_keys[key1]=1
                    frame_shift += 1
                    
                elif((int(keyAnnotationArray[1]) == int(keyNotMatchedArray[1])-1) and (int(keyAnnotationArray[2]) == int(keyNotMatchedArray[2])-1)):
                    finalDictionary[value] = " : frame shift by 1 + " + key1
                    background_dictionary[key] = value
                    tracking_annotation_keys[key1]=1
                    frame_shift += 1
                    
                elif((int(keyAnnotationArray[1]) == int(keyNotMatchedArray[1])+2) and (int(keyAnnotationArray[2]) == int(keyNotMatchedArray[2])+2)):
                   finalDictionary[value] = " : frame shift by 2 - " + key1
                   background_dictionary[key] = value
                   tracking_annotation_keys[key1]=1
                   frame_shift += 1
                   
                elif((int(keyAnnotationArray[1]) == int(keyNotMatchedArray[1])-2) and (int(keyAnnotationArray[2]) == int(keyNotMatchedArray[2])-2)):
                   finalDictionary[value] = " : frame shift by 2 + " + key1
                   background_dictionary[key] = value
                   tracking_annotation_keys[key1]=1
                   frame_shift += 1
                
            else:
                pass
     
    ## predicted by tool but not in Annotation Dictionary       
    for key, value in notMatched_dictionary.items():
        if not (finalDictionary.has_key(value)):  
            finalDictionary[value] ="YY "
            background_dictionary[key] = value
            found_in_tool_but_not_in_annotation += 1
                
                        
    return partial_matched, frame_shift , found_in_tool_but_not_in_annotation

    
def generating_statistics(totalAnnotations, matched, tool_count, partial_matched, frame_shift, In_annotation_but_not_in_tool, found_in_tool_but_not_in_annotation):
    #header = "Anno_count" + "\t" + "matched_tool_ids" + "\n"
    #statistics_file.write(header)
    #line = str(total) + "\t" + str(matched) + "/" + str(total) + "\n"
    statistics_file.write("Total_annotations: " + str(totalAnnotations) + "\n")
    statistics_file.write("Exact_matched: " + str(matched) + "\n")
    statistics_file.write("Partial_matched: " + str(partial_matched) + "\n")
    statistics_file.write("Frame_shift: " + str(frame_shift) + "\n")
    statistics_file.write("In_annotation_but_not_in_tool: " + str(In_annotation_but_not_in_tool) + "\n")
    statistics_file.write("Total_tool_count: " + str(tool_count) + "\n")
    statistics_file.write("Found_in_tool_but_not_in_annotation: " + str(found_in_tool_but_not_in_annotation) + "\n")
    statistics_file.close()
    
def print_final_results():
    for key, value in finalDictionary.items():
        #print key, value
        results_file.write(key + " " + value + "\n")
    results_file.close()

totalAnnotations = converting_annotation_to_dictionary()
tool_count = converting_toolOutput_to_dictionary()
rFile = naming_result_file()
results_file = open(rFile, 'w')
sFile = naming_statistics_file()
statistics_file = open(sFile,'w')
matched = main_comparison_function()
partial_matched , frame_shift, found_in_tool_but_not_in_annotation = to_find_nonMatching()
In_annotation_but_not_in_tool = Found_in_Annotation_but_not_in_tool()

print_final_results()
generating_statistics(totalAnnotations, matched, tool_count, partial_matched, frame_shift, In_annotation_but_not_in_tool, found_in_tool_but_not_in_annotation)

