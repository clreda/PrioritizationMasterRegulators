# coding:utf-8

###############################################################
## Influence maximization algorithm run on a Boolean network ##
## for the detection of master regulator genes               ##
###############################################################

from glob import glob
import pandas as pd
import numpy as np
import subprocess as sb
import sys
import os
from random import seed
from scipy.stats.mstats import gmean
import maboss

from params import *

for a in im_params:
    globals()[a] = im_params[a]

seed(seed_number)
np.random.seed(seed_number)

from utils_grn import load_grn, load_ensemble_grn, set_initial_states, get_states_proba
from utils_state import compare_states

####################
## Spread process ##
####################

#' @param network string: filename of the network in .bnet (needs to be pickable)
#' @param S string list: subset of node names
#' @param gene_outputs string list: list of node names to check
#' @param states Pandas DataFrame or None: if None, pick state_window states at random and geom mean their results (all possible states = too computationnally intensive)
#' @return the mean change in the final attractor states reachable from the modified initial state where genes in S are flipped
## Note: probably not submodular...
def spread(network_name, S, gene_outputs, seednb=0, states=None, verbose=True):
    seed(seednb)
    np.random.seed(seednb)
    master_simulation = load_ensemble_grn("/".join(network_name.split("/")[:-1]), pert_=None, sign_=None, ChromHMM_sig_=None)
    state_len = states.shape[1]
    dists = [None]*state_len
    print_states = None
    gene_outputs_ = [g for g in gene_outputs if (g not in S)]
    for si in range(state_len):
        initial_state = states[[si]]
        mutant_simulation_init = master_simulation.copy()
        ## The considered state is set as initial state
        mutant_simulation_init = set_initial_states(mutant_simulation_init, [initial_state])
        ## Set the gene outputs (genes with at least one regulator different from itself and not perturbed)
        mutant_simulation_init.set_outputs(gene_outputs_)
        mutant_result_init = mutant_simulation_init.run()
        ## Retrieve reachable attractor states from this initial state WITHOUT the current set of perturbations
        attrs_init, probas_init = get_states_proba(mutant_result_init, gene_outputs_)
        ## Directions of perturbations
        mutations_off, mutations_on = [[g for g in S if ((g in states[[si]].dropna().index) and (states.loc[g][si]==v))] for v in [1,0]]
        ## Apply the perturbations
        mutant_simulation = maboss.copy_and_mutate(mutant_simulation_init, mutations_off, "OFF")
        mutant_simulation = maboss.copy_and_mutate(mutant_simulation, mutations_on, "ON")
        mutant_simulation = set_initial_states(mutant_simulation, [initial_state])
        mutant_simulation.set_outputs(gene_outputs_)
        mutant_result = mutant_simulation.run()
        ## Retrieve reachable attractor states from this initial state WITH the current set of perturbations
        attrs, probas = get_states_proba(mutant_result, gene_outputs_)
        sims_, nb_gene = compare_states(attrs, attrs_init, gene_outputs_)
        assert nb_gene == len(gene_outputs_)
        dists[si] = 1-np.max(sims_) #max: minimum of change in (*different*) attractors induced by the subset S
    spd = gmean([d+1 for d in dists])-1 # aggregate the values across the set of initial states
    if (verbose):
        print("Spread = %.4f\tS = %s" % (spd, S))
    return spd

#######################################
## INFLUENCE MAXIMIZATION ALGORITHM  ##
#######################################

# Kempe et al. (2003) greedy algorithm for influence maximization
# Finds iteratively the maximum spreader and adds it to the list until the list is of size k
#' @param network Bonesis network
#' @param states Pandas DataFrame list
#' @param k integer
#' @return S string list: nodes in the spreader set, spreads float list: value of spread associated with each node in spreader set
import random
from joblib import Parallel, delayed

def greedy(network_name, k, genes=None, states=None, ncpus=1, sublist=None, state_len=100):
    S = []
    S_unfold = []
    if (str(genes)=="None"):
        network = load_grn(network_name)
        ## remove nodes without external regulators from the process
        genes = [x.split(" <- ")[0] for x in str(network).split("\n")[:-1]]
        gene_outputs = [x.split(" <- ")[0] for x in str(network).split("\n")[:-1] if (x.split(" <- ")[1] not in [x.split(" <- ")[0], "0", "1"])]
    else:
        gene_outputs = genes
    if (str(sublist)=="None"):
        sublist = genes
    ## Choose binary states at random
    if (states is None):
        states = pd.DataFrame([np.random.choice([0,1], p=[0.5,0.5], size=len(genes)).tolist() for _ in range(state_len)], columns=genes, index=range(state_len)).T
    spreads = {}
    start_k = 0
    spread_values = {}
    if (os.path.exists(root_folder+name+"/application_regulators.csv")):
        res = pd.read_csv(root_folder+name+"/application_regulators.csv", index_col=0)
        if (res.shape[0]>0):
            spreads = res.to_dict()
            S = eval(list(res.columns)[-1])
            S_unfold = [cc for c in eval(list(res.columns)[-1]) for cc in c]
            start_k = res.shape[1]
    for current_k in range(start_k, k):
        print("Iteration #%d" % (current_k+1))
        largest_spread = 0 if (len(S_unfold)==0) else spread(network_name, S_unfold, genes, states=states)
        gene_list = [g for g in sublist if (g not in S_unfold)]
        if (ncpus==1):
            sprds = [spread(network_name, S_unfold+[g], genes, states=states) for g in gene_list]
        else:
            seeds = [np.random.randint(1000000) for _ in range(len(gene_list))]
            sprds = Parallel(n_jobs=ncpus, backend='loky')(delayed(spread)(network_name, S_unfold+[gene_list[gid]], gene_outputs, states=states, seednb=sn) for gid, sn in enumerate(seeds))
        spread_values.update({(current_k+1):{g: sprds[ig] for ig, g in enumerate(gene_list)}})
        ## can't find spreader that strictly increases the spread of the previous spreader
        if (np.max(sprds)<=largest_spread):
            break
        ## Keep ex-aequo
        ex_aequo = np.argwhere(np.array(sprds)==np.max(sprds)).flatten().tolist()
        max_spread, nodes = np.max(sprds), [gene_list[x] for x in ex_aequo]
        largest_spread = max_spread
        if (len(nodes)>1):
            nodes = tuple(nodes)
            S_unfold += nodes
            S += [nodes]
        else:
            nodes = nodes[0]
            S_unfold += [nodes]
            S += [[nodes]]
        spreads.setdefault(str(S), {"SpreadValue": max_spread})
        pd.DataFrame(spreads).to_csv(root_folder+name+"/application_regulators.csv")
        pd.DataFrame(spread_values).to_csv(root_folder+name+'/spread_values.csv')
    spreads_df = pd.DataFrame(spreads)
    gene_list = spreads_df.sort_values(by=spreads_df.columns[0], ascending=False).astype(np.float64)
    gene_list.columns = ["Spread value"]
    vals = pd.DataFrame(spread_values)
    vals.sort_values(by=vals.columns[0], ascending=False).to_csv(root_folder+name+"/spread_values.rnk", sep="\t", header=None)
    print(pd.DataFrame(spreads))
    return S, spreads

#############################
## Application to epilepsy ##
#############################

if (__name__=="__main__"):
    ## Solution to test for master regulator detection
    fname_ls = glob(solution_folder+"solution.bnet")
    assert len(fname_ls)>0
    grn_fname = fname_ls[0]
    grn_fname_reformat = "influence_maximization/"+grn_fname.split("/")[-1]
    sb.call("mkdir -p influence_maximization/", shell=True)
    sb.call("sed 's/ <- /, /g' "+grn_fname+" > "+grn_fname_reformat, shell=True)
    ## Get M30 genes
    with open(grn_fname, "r") as f:
         genes = [x.split(" <- ")[0] for x in f.read().split("\n")[:-1]]
    ## epileptic hippocampi (private communication of normalized count matrix from the raw count data in ArrayExpress)
    if (not os.path.exists(path_to_initial_states)):
        raise ValueError("File does not exists.")
    ## mat: normalized gene expression matrix (rows=genes, columns=samples)
    mat = pd.read_csv(path_to_initial_states, index_col=0)
    mat = mat.loc[~mat.index.duplicated()]
    ## Aggregate and binarize initial state
    from utils_state import binarize_experiments
    bin_mat = binarize_experiments(mat, thres=0.5, method="binary")
    states = bin_mat.astype(int)
    states.columns = range(states.shape[1])
    ## If you wanted to run this on a strict subset of genes
    sublist = None
    state_len = states.shape[1]
    S, spreads = eval(method)(grn_fname_reformat, k=k, sublist=sublist, ncpus=njobs, states=states)
    sb.call("rm -rf influence_maximization/", shell=True)
