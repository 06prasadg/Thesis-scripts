"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Converts output produced by ORF callers to a format suitable to produce Venn diagrams in R
"""

import sys
import re
import csv

fMGA = open(sys.argv[1])
fMGM = open(sys.argv[2])
fOrphelia = open(sys.argv[3])
fGlimmer = open(sys.argv[4])

linesMGA = fMGA.readlines()
linesMGM = fMGM.readlines()
linesOrphelia = fOrphelia.readlines()
linesGlimmer = fGlimmer.readlines()

MGA_list = []
MGM_list = []
Orphelia_list = []
Glimmer_list = []

#list for MGA items
def MGA_function():
    MGA_total=0
    for line in linesMGA:
        if line.startswith(">gi"):
            MGA_total += 1
        lineArr1 = line.split(" ")
        if lineArr1[2] == "00":
            idArr = lineArr1[0]
            matchObject = re.findall(r'_+[0-9]*', line)
            value = ''.join(matchObject[1:3])
            MGA_list.append(value)
    return MGA_total

#list for MGM items
def MGM_function():
    MGM_total=0
    for line in linesMGM:
        if line.startswith(">gi"):
            MGM_total += 1
        lineArr1 = line.split(" ")
        if lineArr1[2] == "00":
            idArr = lineArr1[0]
            matchObject = re.findall(r'_+[0-9]*', line)
            value = ''.join(matchObject[1:3])
            MGM_list.append(value)
    return MGA_total

#list for Glimmer-MG items
def Glimmer_function():
    Glimmer_total=0
    for line in linesGlimmer:
        if line.startswith(">gi"):
            Glimmer_total += 1        
        lineArr1 = line.split(" ")
        if lineArr1[2] == "00":
            idArr = lineArr1[0]
            matchObject = re.findall(r'_+[0-9]*', line)
            value = ''.join(matchObject[1:3])
            Glimmer_list.append(value)
    return Glimmer_total
            
#list for Orphelia items
def Orphelia_function():
    Orpehelia_total=0
    for line in linesOrphelia:
        if line.startswith(">gi"):
            Orpehelia_total += 1 
        lineArr1 = line.split(" ")
        if lineArr1[2] == "00":
            idArr = lineArr1[0]
            matchObject = re.findall(r'_+[0-9]*', line)
            value = ''.join(matchObject[1:3])
            Orphelia_list.append(value)
    return Orpehelia_total

#converting list to csv file
def csv_convertor(listArg,name):
    csvfile = "Csv_output_"+name+".csv"
    with open(csvfile, "w") as output:
         writer = csv.writer(output, lineterminator='\n')
         for val in listArg:
            writer.writerow([val])
    

#finding the common elements
def find_common_elements(MGA_total,MGM_total,Glimmer_total,Orphelia_total):
    
    #pair combinations
    MGA_MGM = set(MGM_list) & set(MGA_list)
    MGA_MGM_count = len(MGA_MGM)
    
    MGA_Glimmer = set(MGA_list) & set(Glimmer_list)
    MGA_Glimmer_count = len(MGA_Glimmer)
    
    MGA_Orphelia = set(MGA_list) & set(Orphelia_list)
    MGA_Orphelia_count = len(MGA_Orphelia)
    
    MGM_Glimmer = set(MGM_list) & set(Glimmer_list)
    MGM_Glimmer_count = len(MGM_Glimmer)
    
    MGM_Orphelia = set(MGM_list) & set(Orphelia_list)
    MGM_Orphelia_count = len(MGM_Orphelia)
    
    Glimmer_Orphelia = set(Glimmer_list) & set(Orphelia_list)
    Glimmer_Orphelia_count = len(Glimmer_Orphelia)
    
    
    #triple combinations
    MGM_MGA_Glimmer = set(MGM_list) & set(MGA_list) & set(Glimmer_list)
    MGM_MGA_Glimmer_count = len(MGM_MGA_Glimmer)
    
    MGM_MGA_Orphelia = set(MGM_list) & set(MGA_list) & set(Orphelia_list)
    MGM_MGA_Orphelia_count = len(MGM_MGA_Orphelia)
    
    MGM_Glimmer_Orphelia = set(MGM_list) & set(Glimmer_list) & set(Orphelia_list)
    MGM_Glimmer_Orphelia_count = len(MGM_Glimmer_Orphelia)
    
    MGA_Glimmer_Orphelia = set(MGA_list) & set(Glimmer_list) & set(Orphelia_list)
    MGA_Glimmer_Orphelia_count = len(MGA_Glimmer_Orphelia)
    
    
    #4s combinations
    MGM_MGA_Glimmer_Orphelia = set(MGM_list) & set(MGA_list) & set(Glimmer_list) & set(Orphelia_list)
    MGM_MGA_Glimmer_Orphelia_count = len(MGM_MGA_Glimmer_Orphelia)
    
    
    #print all numbers
    print "Total MGA : ", MGA_total
    print "Total MGM : ", MGM_total
    print "Glimmer Total : ", Glimmer_total
    print "Orphelia Total : ", Orphelia_total
    #print rest of them all combinations
    
MGA_total = MGA_function()
MGM_total = MGM_function()
Glimmer_total = Glimmer_function()
Orphelia_total = Orphelia_function()
find_common_elements(MGA_total,MGM_total,Glimmer_total,Orphelia_total)

##call to csv generator file function
csv_convertor(MGA_list,"MGA")
csv_convertor(MGM_list,"MGM")
csv_convertor(Glimmer_list,"Glimmer")
csv_convertor(Orphelia_list,"Orphelia")

        
    
