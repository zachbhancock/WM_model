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

def effective_dispersal(ts, max_edge_length=math.inf):
    """
    Returns the mean standardized squared distance (in x and y) between all
    parent-child relationships in the recorded tree sequence separated
    by a time interval less than or equal to `max_edge_length`. This is an
    estimator for the variance of a Gaussian dispersal kernel.
    """
    assert type(ts) == tskit.trees.TreeSequence
    nx = ny = 0
    mx = my = 0
    for edge in ts.edges():
        parent_id = edge.parent
        child_id = edge.child
        edge_length = ts.node(parent_id).time - ts.node(child_id).time
        if edge_length <= max_edge_length:
            parent_individual_id = ts.node(parent_id).individual
            child_individual_id = ts.node(child_id).individual
            if parent_individual_id and child_individual_id:
                parent_loc = ts.individual(parent_individual_id).location
                child_loc = ts.individual(child_individual_id).location
                if len(parent_loc) >= 2 and len(child_loc) >= 2:
                    dx = (parent_loc[0] - child_loc[0]) / math.sqrt(edge_length)
                    dy = (parent_loc[1] - child_loc[1]) / math.sqrt(edge_length)
                    if nx:
                        nx += 1
                        mx += (dx*dx - mx) / nx
                    else:
                        nx = 1
                        mx = dx*dx
                    if ny:
                        ny += 1
                        my += (dy*dy - my) / ny
                    else:
                        ny = 1
                        my = dy*dy
    return [mx, my]



eff_dispersal = effective_dispersal(ts)

def mean(eff_dispersal):
  return sum(eff_dispersal) / len(eff_dispersal)

eff_dispersal = mean(eff_dispersal)
eff_dispersal = str(eff_dispersal)
with open(prefix+"-eff_disp.txt","w") as indfile:
  indfile.write(eff_dispersal)


generation_ago = pyslim.individuals_alive_at(ts, 1)
eff_density = len(generation_ago) / 625
eff_density = str(eff_density)
with open(prefix+"-eff_dens.txt","w") as indfile:
  indfile.write(eff_density)


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