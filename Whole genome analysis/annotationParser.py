"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
This script converts whole genome annotation predictions to a common standard format
"""

from sys import argv
file1 = argv[1]

f=open(file1,'r')
f1=open("coordinates.fasta",'w')

newLine= ''
protein_id = ''
num1 = ''
num2 = ''

for i, line in enumerate(f):
    if i==0:
        header=line.strip()
    if i > 0 :
        array1= line.split('\t')
        #print array1
        if len(array1) == 3 and array1[2] == "CDS\n":
            num1 = array1[0]
            num2 = array1[1]
            #print num1 , num2
        if len(array1) > 3 and array1[3] == "protein_id":
            protein_id = array1[4].strip()
            #print protein_id
            newLine = header + protein_id+ "_" + num1 + "_" + num2 + "\n"
            print newLine
            #f1.write(newLine)
     
