import sys

topfile = sys.argv[1]
outname = sys.argv[2]

top = open(topfile,'r')
out = open(outname,'w')

line=top.readline()
while line:
	if not '%FLAG SCEE_SCALE_FACTOR' in line:
		out.write(line)
	else:
		while not '%FLAG SOLTY' in line:
			line = top.readline() 
		out.write(line)
	line = top.readline()
out.close()
top.close()
print "DONE"
