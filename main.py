K = 25
transcripts = [] 
read = []

""" The input should have three files, 
parameters.txt : K (the K in k-mer)
transcripts.txt : every non-empty line is a transcript
read.txt : the read.  
in read and transcripts, only ACGT characters are considered, other characters are ignored
"""
K = int (open('parameters.txt', 'r').read())
#raw = filter(lambda x : x in ['A','C','G','T', '\n'], open('transcripts.txt','r').read())
raw = open('transcripts.txt','r').read()
raw = filter(lambda x : len(x) > 0, raw.split('\n'))

transcripts = raw[1::2]
names = raw[::2]

import re, math, pickle, sys
def occurrences(text, sub):
  return len(re.findall('(?={0})'.format(re.escape(sub)), text))

NEW_INDEX = False
if len(sys.argv) > 1 : 
  NEW_INDEX = True

if NEW_INDEX :
  print "creating index"
  all_kmers = []
  for line in transcripts  :
    kmers = set([line[i:i+K] for i in range(len(line) - K + 1)])
    all_kmers += kmers
  sortedKmers = sorted(list(all_kmers))
  print "Sorted : ", len(sortedKmers)
  apperance = [[0 for y in range(len(transcripts))] for x in range(len(sortedKmers))]
  for i, kmer in enumerate(sortedKmers) :
    for j, trans in enumerate(transcripts) :
      apperance[i][j] += occurrences(trans, kmer)
  EquivClass = [-1] * len(sortedKmers)
  NumOfEqCls = 0
  for i in range(len(sortedKmers)) :
    if EquivClass[i] >= 0: continue
    EquivClass[i] = NumOfEqCls
    for j in range(i + 1, len(sortedKmers)) :
      if all((x == y for x, y in zip(apperance[i],  apperance[j]))) :
        EquivClass[j] = NumOfEqCls
    NumOfEqCls += 1
  print "Number of equivlent class :", NumOfEqCls

  pickle.dump(sortedKmers, open('index1.sav', 'w'))
  pickle.dump(EquivClass, open('index2.sav', 'w'))
  pickle.dump(apperance, open('index3.sav', 'w'))

else :
  sortedKmers = pickle.load(open('index1.sav', 'r'))
  EquivClass = pickle.load(open('index2.sav', 'r'))
  NumOfEqCls = max(EquivClass) + 1
  apperance = pickle.load(open('index3.sav', 'r'))

hashtable = dict()
for i in range(len(sortedKmers)) :
  hashtable[sortedKmers[i]] = i

# process read
#read = filter(lambda x : x in ['A', 'C', 'G', 'T'], open('read.txt').read())
read = open('reads.txt').read()
L = [0.0] * NumOfEqCls
mapped = 0
for i in range(len(read) - K + 1) :
  kmer = read[i:i+K]
  if kmer in hashtable :
    j = hashtable[kmer]
    #print NumOfEqCls, j, EquivClass[j]
    L[EquivClass[j]] += 1
    mapped += 1

"""
alpha[j, i] := kmer j allocate to transcript i
"""

# initialize
alpha = [[0.0] * len(transcripts)] * NumOfEqCls
u1 = [0.0] * len(transcripts)
u = [0.0] * len(transcripts)
kmer_occurence = [sum(apperance[i][j] for j in range(len(transcripts))) 
                    for i in range(len(sortedKmers))]
apperance_total = [0.0] * NumOfEqCls
for l, kmer in enumerate(sortedKmers) :
  apperance_total[EquivClass[l]] += kmer_occurence[l]
for l, kmer in enumerate(sortedKmers) :
  j = EquivClass[l]
  for i in range(len(transcripts)) :
    alpha[j][i] += L[j] * apperance[l][i] / apperance_total[j]
  
# apperance : kmer trans
set_of_trans_has_kmer = [0] * NumOfEqCls
set_of_kmer_in_trans = [0] * len(transcripts)
for j in range(NumOfEqCls) :
  set_of_trans_has_kmer[j] = filter(lambda t : apperance[j][t] > 0, range(len(transcripts)))
for i in range(len(transcripts)) :
  setofkmers = filter(lambda j : apperance[j][i] > 0, range(len(sortedKmers)))
  set_of_kmer_in_trans[i] = list(set( EquivClass[j] for j in setofkmers))

l_prime = [len(t) - K + 1 for t in transcripts]
print "initialization done"

def step_get_alpha(u1) :
  global set_of_trans_has_kmer
  for j in range(NumOfEqCls) :
    w = sum(u1[t] for t in set_of_trans_has_kmer[j])
    for i in set_of_trans_has_kmer[j] :
      alpha[j][i] = u1[i] * L[j] / w
  return alpha

def step_get_mu1(alpha) :
  global set_of_kmer_in_trans
  y = 0.0
  u_star = [0.0] * len(transcripts)
  for t, trans in enumerate(transcripts) :
    c_t = sum(alpha[j][t] for j in set_of_kmer_in_trans[t])
    u_star[t] = c_t / l_prime[t]
    y += u_star[t]
  u11 = [ ut / y for ut in u_star ]
  return u11


def EMstep(mu1) :
  alpha = step_get_alpha(mu1)
  return step_get_mu1(alpha)

def norm(vec) :
  s = sum(map(lambda x : x * x, vec))
  return math.sqrt(s)
  
def SQUAREMstep(u0) :
  u1 = EMstep(u0)
  u2 = EMstep(u1)
  r = [a1-a0 for a1, a0 in zip(u1,u0)]
  v = [a2-a1-rr for a2, a1,  rr in  zip(u2,u1,r)] 
  gamma = -norm(r) / norm(v)
  #gamma = modify(gamma) ;; currently no such step
  temp = [u_0 - 2.0 * gamma * r_ + gamma*gamma*v_ for u_0, r_, v_ in zip(u0, r, v)]
  u3 = [max(0.0, tmp) for tmp in temp]
  return EMstep(u3)

def run_EM(u1) :
  for i in range(100) : # limited iterations
    #u1 = SQUAREMstep(u1)
    print u1
    u1 = EMstep(u1)
  return u1


u1 = step_get_mu1(alpha)
u1 = run_EM(u1)
for t in range(len(transcripts)) :
  print names[t], u1[t], u1[t] * 1000000, u1[t] * (10 ** 9) / mapped

