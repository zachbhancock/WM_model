#python bs

import os
import numpy as np
import msprime, pyslim
import tskit
import sys
import pandas as pd
import math
import glob
from numpy.random import default_rng

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

def min(a, b):
     
    if a <= b:
        return a
    else:
        return b

def effective_dispersal_torus(ts, max_edge_length=math.inf):
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
                    dx_wrapped = (25 - dx)/math.sqrt(edge_length)
                    dy_wrapped = (25 - dy)/math.sqrt(edge_length)
                    dx = min(dx_wrapped, dx)
                    dy = min(dy_wrapped, dy)
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



def eff_N(file_path):
   f = open(file_path+"_popsize", "r")
   N = f.read()
   N = int(N)
   ts = tskit.load(file_path)
   ts = ts.simplify()
   recap_ts = pyslim.recapitate(ts, recombination_rate=1e-8, ancestral_Ne=2*N)
   ts = msprime.sim_mutations(recap_ts, rate=1e-7, random_seed=12345)
   ancestors = pyslim.individuals_alive_at(ts, 1)
   ancestor_N = len(ancestors)
   area = 2500
   eff_dispersal = effective_dispersal_torus(ts)
   def mean(eff_dispersal):
     return sum(eff_dispersal) / len(eff_dispersal)
   
   eff_dispersal = mean(eff_dispersal)
   eff_WN = 4*math.pi*((2*ancestor_N)/area)*eff_dispersal
   eff_WN = str(eff_WN)
   with open(file_path+"-eff_N.txt","w") as indfile:
     indfile.write(eff_WN)

list_of_prefixes = glob.glob('wmModel_*_slimIter_*_sigma_*_K_*0')

for file_path in list_of_prefixes:
    # Make sure it's a file and not a directory
    if os.path.isfile(file_path):
        # Run your script on the file
        eff_N(file_path)
  

exit()