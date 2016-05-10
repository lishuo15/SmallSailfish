import sys

def readfile(filename) :
  with open(filename,'r') as f:
    s = f.read().split('>')
  ret = []
  for one in s :
    if len(one) > 3 :
      k = one.splitlines()
      print k
      ret.append( (k[0], ''.join(k[1:]) ) )
  return  ret

if (len(sys.argv) > 2) : 
  outf = open(sys.argv[2], 'w')
else :
  outf = open('transcripts.txt', 'w')
for x,y in readfile(sys.argv[1]) :
  if len(x+y) > 3 :
    outf.write(x + '\n' + y + '\n')
outf.close()
