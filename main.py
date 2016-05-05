K = 25
transcripts = [] 
read = []

def read_input() :
  """ The input should have three files, 
  parameters.txt : K (the K in k-mer)
  transcripts.txt : every non-empty line is a transcript
  read.txt : the read.  
  in read and transcripts, only ACGT characters are considered, other characters are ignored
  """
  K = int (open('parameters.txt', 'r').read())
  raw = filter(lambda x : x in ['A','C','G','T', '\n'], open('transcripts.txt','r').read())
  transcripts = filter(lambda x : len(x) > 0, raw.split('\n'))
  read = filter(lambda x : x in ['A', 'C', 'G', 'T'], open('read.txt'))

  all_kmers = set()
  for line in transcripts  :
    kmers = set([line[i:i+K] for i in range(len(line) - K + 1)])
    all_kmers += kmers
  sortedKmers = sorted(list(all_kmers))
  hashtable = dict()
  for i in range(len(sortedKmers)) :
    hashtable[sortedKmers[i]] = i
  apperence = [[0] * len(transcripts) for x in range(len(sortedKmers))]
  for i, kmer in enumerate(sortedKmers) :
    for j, trans in enumerate(transcripts) :
      if kmer in trans :
        apperance[i][j] += 1
  EquivClass = [-1] * len(sortedKmers)
  NumOfEqCls = 0
  for i in range(len(sortedKmers)) :
    if EquivClass[i] > 0 continue
    NumOfEqCls += 1
    EquivClass[i] = NumOfEqCls
    for j in range(i + 1, len(sortedKmers)) :
      if all((x == y for x, y in zip(apperence[i],  apperence[j]))) :
        EquivClass[j] = NumOfEqCls
  EquivKmer = filter(lambda x : len(x > 0), 
      map(lambda i : sortedKmers[i] if EquivClass[i] == i else "", range(len(sortedKmers))))

def count_kmers() :
  L = [0] * NumOfEqCls
  for idx, kmer in enumerate(sortedKmers) :
    L[EquivClass[idx]] += read.count(kmer)

def init_mu() :
  # fix later
  pass
  
def EMstep(mu1) :
  alpha = step_get_alpha(mu1)
  return step_get_mu1(alpha)

def norm(vec) :
  s = sum(map(lambda x : x * x, vec))
  return sqrt(s)
  
def SQUAREMstep(u0) :
  u1 = EMstep(u0)
  u2 = EMstep(u1)
  r = [a1-a0 for a1, a0 in zip(u1,u0)]
  v = [a2-a1-rr for a2, a1,  rr in  zip(u2,u1,r)] 
  gamma = -norm(r) / norm(v)
  #gamma = modify(gamma) ;; currently no such step
  u3 = [0.0] * len(u1)
  return EMstep(u3)

def run_EM(u1) :
  for i in range(10000) : # limited iterations
    u1 = SQUAREMstep(u1)
  return u1

def output() :
  global mu, mu1
  print mu
  print mu1

if __name__ == 'main' :
  read_input()
  build_index()
  count_kmers()
  run_EM()
  output()


# TO DO : vector operations
# Initialization
# confirm step 6 and 7 in SQUAREM
# Test on data
