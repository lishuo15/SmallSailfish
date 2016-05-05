import sys

def readfile(filename) :
  with open(filename,'r') as f:
    s = f.read().split('>')
  ret = []
  for one in s :
    k = one.splitlines()
    ret.append(''.join(k[1:]))
  return  ret

if (len(sys.argv) > 2) : 
  outf = open(sys.argv[2], 'w')
else :
  outf = open('transcripts.txt', 'w')
for x in readfile(sys.argv[1]) :
  if len(x) > 3 :
    outf.write(x)
    outf.write('\n')


