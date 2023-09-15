#python bs

import os
import numpy as np
import msprime, pyslim
import tskit
import sys
import pandas as pd
import math
from numpy.random import default_rng

rng = default_rng()

#pull in variables from .sbatch file
#print(len(sys.argv))
#print(str(sys.argv))
WORKINGDIR = sys.argv[1]
TREEFILE = sys.argv[2]
prefix = WORKINGDIR + TREEFILE
f = open(prefix+"_popsize", "r")
N = f.read()
N = int(N)
print(f"WORKINGDIR is ", WORKINGDIR)
print(f"TREEFILE is ", TREEFILE)
print(f"prefix is ", prefix)
print(f"N is ", N)

ts = tskit.load(prefix)
ts = ts.simplify()
recap_ts = pyslim.recapitate(ts, recombination_rate=1e-8, ancestral_Ne=2*N)
ts = msprime.sim_mutations(recap_ts, rate=1e-8, random_seed=12345)

afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
N = ts.num_individuals

freq = []
for i in range(len(afs)):
   a = (afs[i]*((i/(N))**2+(1-(i/(N)))**2))
   freq.append(a)

s = np.sum(freq)/np.sum(afs)

with open(prefix+"-true_s.txt","w") as indfile:
  indfile.write(s)

exit()