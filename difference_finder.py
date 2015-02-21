from sys import argv

file_2014 = argv[1]
file_2013 = argv[2]
diff_file = argv[3]

f_2014 = open(file_2014,'r')
f_2013 = open(file_2013,'r')
d_file = open(diff_file,'w')


f_2014_readlines = f_2014.readlines()
f_2013_readlines = f_2013.readlines()

dict_2013 = {}
dict_2014 = {}

def diff_finder():
    for line in f_2013_readlines:
        split_line1 = line.split('\t')
        key1 = split_line1[0]
        dict_2013[key1] = 1
        
    for line in f_2014_readlines:
        split_line = line.split('\t')
        key2 = split_line[0]
        dict_2014[key2] = 1
            
    for key1 in dict_2013.keys():
        if dict_2014.has_key(key1):
            pass
        else:
            d_file.write(key1 + '\n') 
            
def main():
    diff_finder()
    
if __name__ == '__main__':
	main()
