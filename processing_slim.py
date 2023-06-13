#python bs

import os
import numpy as np
import msprime, pyslim
import tskit
import sys
import pandas as pd
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
recap_ts = pyslim.recapitate(recombination_rate=1e-8, ancestral_Ne=2*N)
ts = msprime.sim_mutations(recap_ts, rate=1e-7, random_seed=12345)

#specific places
alive = pyslim.individuals_alive_at(ts, 0)

groups = {
    'alive' : np.random.choice(alive, size=100, replace=False),
}

group_order = ['alive']
sampled_nodes = [[] for _ in groups]
for j, k in enumerate(group_order):
   for ind in groups[k]:
      sampled_nodes[j].extend(ts.individual(ind).nodes)

ind_nodes = []
ind_group = []
ind_ids = []
for j, group in enumerate(group_order):
   for ind in groups[group]:
      ind_ids.append(ind)
      ind_nodes.append(ts.individual(ind).nodes)
      ind_group.append(group_order[j])

nind = len(ind_ids)
pairs = [(i, j) for i in range(nind) for j in range(nind)]
ind_div = ts.divergence(ind_nodes, indexes=pairs)

x = []
for i in ind_ids:
   ind = ts.individual(i)
   label = f"tsk_{ind.id}"
   x.append(label)

b = np.reshape(ind_div, (100,100))
panda_df = pd.DataFrame(data = b, columns = x)
panda_df.to_csv(prefix+"-pi.csv", sep=' ', index=True)

ind_Fst = ts.Fst(ind_nodes, indexes=pairs)

c = np.reshape(ind_Fst, (100,100))
panda2_df = pd.DataFrame(data = c, columns = x)
panda2_df.to_csv(prefix+"-Fst.csv", sep=" ", index=True)

indivlist = []
indivnames = []
with open(prefix+"-pi_locs.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for group in group_order:
     for i in groups[group]:
        indivlist.append(i)
        ind = ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

exit()