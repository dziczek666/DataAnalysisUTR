import sys

# read program parameters from the command line
GTF = sys.argv[1] 
input_class = sys.argv[2] 
out = open(sys.argv[3], 'w')  


for line in open(GTF):
  if line.split('\t')[8] == input_class:
    out.write(line)
    