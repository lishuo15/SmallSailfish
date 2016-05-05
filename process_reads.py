import sys

def  readfile(filename) :
  with open(filename,'r') as f :
    s = f.read().split('@')
  ret  =  []
  for x in s :
    if len(x) < 10 : continue
    k = x.splitlines()
    ret.append(k[1])
  return  ret

outf = open('reads.txt', 'w')

s = readfile(sys.argv[1])
outf.write('I'.join(s))
outf.close()
