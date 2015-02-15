"""
Written by Prasad Gajare, Center for Bioinformatics and Computational Biology, Delaware Biotechnology Institute, University of Delaware.
Please report bugs at prasadg@udel.edu
Script to find ORFans from the two input files
"""


from sys import argv


btab_file = argv[1]
fasta_file = argv[2]
#orfan_btab_file = argv[3]

btab = open(btab_file,'r')
fasta = open(fasta_file,'r')
#orfan = open(orfan_btab_file, 'w')

btab_readlines = btab.readlines()
fasta_readlines = fasta.readlines()

btab_dict = {}
fasta_dict = {}

def ORFan_finder():
	# for btab file
	for line in btab_readlines:
		split_line = line.split('\t')
		btab_dict[">" + split_line[0]] = 1
        # for fasta file
	value1 = ''
	value2 = ''
	for line in fasta_readlines:
		if line.startswith('>'):
			value1 = ''
			value2 = ''
			split_line = line.split(' ')
			key = split_line[0]
			value1 = line.strip()
		else:
			value2 = value2 + line
		fasta_dict[key] = value1 + value2
	
	for k in fasta_dict.keys():
		if k in btab_dict:
			pass
		else:
			print fasta_dict[k]

ORFan_finder()

