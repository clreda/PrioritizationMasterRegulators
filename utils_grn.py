# coding: utf-8

from params import njobs, file_folder, name, maboss_params

import bonesis
bonesis.settings["parallel"] = njobs
bonesis.settings["solutions"] = "subset-minimal"
import mpbn

import numpy as np
import pandas as pd
from itertools import product
from functools import reduce
import sys
import subprocess as sb
from zipfile import ZipFile
import os
from glob import glob
from tqdm import tqdm

import maboss

from utils_data import quantile_normalize

## [Babichev et al., 2019]
## https://cran.r-project.org/web/packages/desirability/vignettes/desirability.pdf
def desirability(f, A=0,B=1,s=1):
    return np.exp(-np.exp(-((B-A)+s/(B-A)*f)))
def gmean(iterable):
    a = np.log(iterable)
    return np.exp(np.mean(a))
## symmetrize influences
def symmetrize_influences(influences__):
    genes = list(set(list(influences__.index)+list(influences__.columns)))
    influences_di = influences__.to_dict()
    for t in [tt for tt in genes if (tt not in influences__.columns)]:
        influences_di.setdefault(t, {g : 0 for g in genes})
    influences_di = pd.DataFrame(influences_di).T.to_dict()
    for g in [gg for gg in genes if (gg not in influences__.index)]:
        influences_di.setdefault(g, {t : 0 for t in genes})
    influences = pd.DataFrame(influences_di).T.loc[genes][genes]
    return influences

#' @param influences data table genes x genes with 1's (resp. -1's) for positive (resp. negative) interactions
#' @param perc_influences data table genes x genes with percentage of observation of interaction (comprised between 0 and 100)
#' @return Harrington's desirability index value
def general_topological_parameter(influences_, perc_influences_=None, weights=None):
    assert str(weights)!="None"
    influences = symmetrize_influences(influences_).fillna(0)
    if (str(perc_influences_)!="None"):
        perc_influences = symmetrize_influences(perc_influences_)
        perc_influences = perc_influences.loc[influences.index][influences.columns]
    else:
        perc_influences = perc_influences_
    genes = list(set(list(influences_.index)+list(influences_.columns)))
    influences_di = influences_.to_dict()
    for t in [tt for tt in genes if (tt not in influences_.columns)]:
        influences_di.setdefault(t, {g : 0 for g in genes})
    influences_di = pd.DataFrame(influences_di).T.to_dict()
    for g in [gg for g in genes if (g not in influences_.index)]:
        influences_di.setdefault(g, {t : 0 for t in genes})
    influences = pd.DataFrame(influences_di).T.loc[genes][genes]
    n = len(list(influences.index))
    network_outdegree = np.sum(np.abs(influences.values),axis=1)
    k_mean = np.mean(network_outdegree)
    k_max = np.max(network_outdegree)
    ## To maximize (in [0,1])
    DS = 2*np.sum(np.sum(np.triu(np.abs(influences.values)),axis=1),axis=0)/float(n*(n-1))
    ## To maximize (in [0,1])
    connect = np.sum(np.triu((np.abs(influences)+np.abs(influences).T).values),axis=1)
    max_connect = 0.5*np.multiply(network_outdegree,network_outdegree-1)
    max_connect[max_connect<=0] = 0.5
    numCL = np.sum(np.divide(connect,max_connect),axis=0)
    CL = numCL/n
    ## To maximize (in [0,+inf])
    Centr = n/(n-2)*(k_max/(n-1)-DS)
    ## To maxmize (in [0,+inf])
    GT = np.sqrt(np.var(network_outdegree))/max(k_mean,1)
    criteria = ["DS","CL","Centr","GT"]
    values = [DS, CL, Centr, GT]
    Ds = [desirability(values[if_], s=weights[f]) for if_, f in enumerate(criteria)]
    return gmean(Ds)

def finetuning_edge_selection(influences, thres_range=range(1,10,1), verbose=True):
    best_criterion = -float("inf")
    best_thres = None
    thres_range = list([x*0.1 for x in thres_range])
    for thres in thres_range:
        influences_ = np.zeros((influences.shape[0],influences.shape[1]))
        for u,v in np.argwhere(np.abs(influences.values)>=thres).tolist():
            influences_[u,v] = abs(influences.iloc[u,v])
        influences_ = pd.DataFrame(influences_,index=influences.index,columns=influences.columns)
        criterion = general_topological_parameter(influences_)
        if (verbose):
            print("Thres=%.2f\tD=%.5f" % (thres,criterion))
        if (criterion>=best_criterion):
            best_criterion = criterion
            best_thres = thres
    influences_ = np.zeros((influences.shape[0],influences.shape[1]))
    for u,v in np.argwhere(np.abs(influences.values)>=best_thres).tolist():
        influences_[u,v] = (-1)**int(influences.iloc[u,v]<0)
    influences = pd.DataFrame(influences_,index=influences.index,columns=influences.columns)
    if (verbose):
        print("Best value=%.5f for thres=%.2f" % (best_criterion,best_thres))
    return influences

## Convert influences DataFrame to GRN using the default regulatory function (inspired by CasQ)
## if target <- activators, inhibitors, then grf(target) = and(a, a in activators) AND and(not i, i in inhibitors)
#' @param influences Pandas DataFrame
#' @param fname string
def influences2bnet(influences, fname, sep="<-", write=True, default=0):
    genes = []
    grfs = {}
    for nx,ny in np.argwhere(influences.values!=0).tolist():
        x,y,v = influences.index[nx], influences.index[ny], influences.values[nx,ny]
        genes = list(set(genes+[x,y]))
        regs_y = grfs.get(y, [])
        val = x if (v>0) else "(!"+x+")"
        if (val not in regs_y):
            regs_y.append(val)
        grfs.update(dict([[y,regs_y]]))
    bnet = []
    for g in grfs:
        bnet.append(g+" "+sep+" "+"&".join(grfs[g]))
    for g in genes:
        if (g not in grfs):
            bnet.append(g+" "+sep+" "+str(default))
    bnet = "\n".join(bnet)
    if (write):
        with open(fname, "w+") as f:
            f.write(bnet)
    return bnet

def load_ensemble_grn(path, pert_=None, sign_=None, ChromHMM_sig_=None, params=maboss_params):
    assert os.path.exists(path)
    simulation = maboss.Ensemble(path=path)
    mutation_di = {}
    if (str(ChromHMM_sig_)!="None"):
        sig = ChromHMM_sig_.dropna().astype(int)
        sig[sig == 0] = "OFF"
        sig[sig == 1] = "ON"
        sig = sig.loc[[g for g in sig.index if (g in simulation.nodes)]]
        if (len(sig)>0):
            mutation_di.update(sig.to_dict()[list(sig.keys())[0]])
    if ((pert_!="None") and (str(sign_)!="None") and (pert_ in simulation.nodes)):
        mutation_di.update(dict([[pert_,sign_]]))
    if (len(mutation_di)>0):
        simulation.mutations = mutation_di
    simulation.param.update(params)
    simulation.set_outputs([g for g in simulation.nodes if (g not in simulation.mutations)]) 
    return simulation

def set_initial_states(simulation, initial_states):
    initial = pd.concat(initial_states, axis=1).dropna()
    initial = initial.loc[[g for g in initial.index if ((g not in simulation.mutations) and (g in simulation.nodes))]]
    prob = initial.mean(axis=1)
    for g in prob.index:
        p = float(prob.loc[g])
        simulation.set_istate(g, (1.-p, p))
    for g in [gene for gene in simulation.nodes if (gene not in initial.index)]:
        simulation.set_istate(g, (0.5,0.5))
    return simulation

def get_states_proba(result, genes):
    probas = result._get_raw_probas()
    states = result._get_raw_states()
    probs = list(probas[-1])
    N = len(probs)
    res_mat = np.matrix([[int(g in s) for s in states[-1]] for g in genes])
    res_states = pd.DataFrame(res_mat, index=genes, columns=["S%d" % i for i in range(len(states[-1]))]).copy()
    return res_states, np.array(probs)

## weighted influences
def save_ensemble_grn(fname, influences, thres_range=None, num=100, encode=False):
    if (str(thres_range)=="None"):
        vals = np.sort(np.unique(influences.loc[influences.values!=0].values)).tolist()
        n1, n2 = max((len(vals)//2)-(num//2), 0), min((len(vals)//2)+(num//2), len(vals))
        thres_range = vals[n1:n2]
        if(len(thres_range)!=num):
            thres_range = np.linspace(vals[0],vals[1],num,endpoint=True).tolist()
        assert len(thres_range)==num
    grn_list = [None]*len(thres_range)
    for it, thres in enumerate(thres_range):
        influences_ = np.multiply(np.array(np.abs(influences.fillna(0).values)>=thres, dtype=int),(-1)**np.array(influences.fillna(0).values<=-thres, dtype=int))
        influences_ = pd.DataFrame(influences_,index=influences.index,columns=influences.columns)
        nedges = np.sum(np.abs(influences_.values))
        if (nedges > 0):
            grn_list[it] = influences2bnet(influences_, "", sep=",", write=False)
    grn_list = [grn for grn in grn_list if (str(grn)!="None")]
    with ZipFile(fname, "w") as bundle:
        for i, bn in enumerate(grn_list):
            with bundle.open(f"bn{i}.bnet", "w") as fp:
                fp.write(str.encode(bn) if (encode) else bn)
    return None

## DFS on undirected network
def get_weakly_connected(edges, genes):
    N = len(genes)
    components = []
    adjacency = np.zeros((N,N))
    ## build undirected adjacency matrix
    for g1,g2 in edges:
        ig1,ig2 = genes.index(g1),genes.index(g2)
        adjacency[ig1,ig2] = adjacency[ig2,ig1] = 1
    ## DFS
    to_visit = [N-1]
    visited = [False]*N
    while (not all(visited)):
        component = []
        while (True):
            node = to_visit.pop()
            visited[node] = True
            if (genes[node] not in component):
                component.append(genes[node])
            children = np.argwhere(adjacency[node,:]==1).flatten().tolist()
            to_visit = [child for child in children if (not visited[child])]+to_visit
            if (len(to_visit)==0):
                break
        components.append(component)
        to_visit = [np.argmin(visited)]
    return components

#' @param path_to_PPI path to STRING PPI
#' @return list of genes, undirected interactions
def get_genes_interactions_from_PPI(path_to_PPI, connected=False, score_STRING=0, score="score"):
    ppi = pd.read_csv(path_to_PPI, sep="\t", header=0)
    ppi.index = ["-".join(list(sorted(x))) for x in zip(ppi["preferredName_A"], ppi["preferredName_B"])]
    ppi = ppi.sort_values(by=score, ascending=False)
    ppi = ppi[~ppi.index.duplicated(keep="first")]
    def ppi2edges(ppi_input):
        edges_output = list(zip(list(ppi_input["preferredName_A"]),list(ppi_input["preferredName_B"])))
        return edges_output
    edges = ppi2edges(ppi)
    genes = list(set([g for e in edges for g in e]))
    if (connected):
        components = get_weakly_connected(edges, genes)
        main_component_id = np.argmax([len(c) for c in components])
        ## Genes which are isolated in the full STRING PPI
        isolated_genes = []
        if (len(components)>1):
            for c_id,c in components:
                if (c_id != main_component_id):
                    isolated_genes += c
        for g in isolated_genes:
            ppi = ppi.loc[(ppi["preferredName_A"]!=g)&(ppi["preferredName_B"]!=g)]
        Ngenes = len(components[main_component_id])
        t = min(score_STRING/1000., ppi[score].max())
        ppi_ = ppi.loc[ppi[score]>=t]
        edges_ = ppi2edges(ppi_)
        genes_ = list(set([g for e in edges_ for g in e]))
        components = get_weakly_connected(edges_, genes_)
        main_component_id = np.argmax([len(c) for c in components])
        isolated_genes = [g for g in genes if (g not in genes_)]
        if (len(components)>1):
            for c_id,c in enumerate(components):
                if (c_id != main_component_id):
                    isolated_genes += c
        edges = [e for e in edges_]
        ## Add edges in decreasing score order until the PPI is connected and all genes are present
        if (len(isolated_genes)>0):
            ppi_ = ppi.loc[ppi[score]<t]
            if (ppi_.shape[0]>0):
                ids_isolated = [idx for idx in ppi_.index if ((ppi_.loc[idx]["preferredName_A"] in isolated_genes) or (ppi_.loc[idx]["preferredName_B"] in isolated_genes))]
                ppi_ = ppi_.loc[ids_isolated]
                ppi__scores = list(set(ppi_[score]))
                i, done = 0, False
                edges__ = ppi2edges(ppi_)
                while (not done and i < len(ppi__scores)):
                    edges__ = ppi2edges(ppi_.loc[ppi_[score]==ppi__scores[i]])
                    edges += edges__
                    components = get_weakly_connected(edges, genes)
                    done = (len(components)==1) and (len(components[0])==Ngenes)
                    i += 1
        assert len(components)==1 and len(set([g for c in components for g in c]))==Ngenes
        genes = list(set([g for e in edges for g in e]))
    gene_a, gene_b = [e for e,_ in edges], [e for _,e in edges]
    return genes, [gene_a, gene_b]

#' @param influences pandas DataFrame genes x genes
#' @return influences trimmed out isolated nodes
def trim_isolated_out(influences):
    non_isolated = [u for nu, u in enumerate(influences.index) if (np.sum(np.abs(influences.values[nu,:]))>0 or np.sum(np.abs(influences.values[:,nu]))>0)]
    influences = influences.loc[non_isolated][non_isolated]
    return influences

#' @param influences pandas DataFrame genes x genes
#' @param fname filename of png file
#' @param optional should interactions be drawn as optional (dashed lines)
#' @return None
# change engine by checking out man dot for graphviz
def influences2graph(influences, fname, optional=False, compile2png=True, engine="sfdp"):#"dot"):
    dotfile = fname+".dot"
    filename = fname+".png"
    graph = ["digraph {"]
    for idx in influences.columns:
        graph += [idx+" [shape=circle];"]
    for x,y in np.argwhere(influences.values != 0).tolist():
        nx, ny = influences.index[x], influences.columns[y]
        sign = influences.values[x,y]
        signs = [-1,1] if (sign == 2) else [sign]
        for sign in signs:
            edge = nx+" -> "+ny
            edge += " ["
            edge += "label="+("\"+\"" if (sign > 0) else "\"-\"")
            edge += ", color="+("red" if (sign < 0) else "green")
            edge += ",arrowhead=tee" if (sign < 0) else ""
            edge += ",style=dashed" if (optional) else ""
            edge += "];"
            graph += [edge]
    graph += ["}"]
    with open(dotfile, "w+") as f:
        f.write("\n".join(graph))
    if (compile2png):
        sb.call("sfdp -Tpng "+dotfile+" > "+filename+" && rm -f "+dotfile, shell=True)
    return None

#' @param gene_list Python list of genes
#' @param pregraph influences matrix or None (all interactions in that matrix are preserved)
#' @param expr_mat_file file name of experiment files
#' @param interactions known interactions (list of edge input genes, list of edge output genes)
#' @param full_graph consider the full graph (N^2 edges for N nodes)
#' @param tau threshold for correlation noise
#' @param direct_graph try to direct edges using the CBDN method
#' @param cor_method string character for method implemented in scipy
#' @param verbose print stuff
#' @param beta power for the adjacency graph: the higher it is, the more stronger edges are selected
#' @return influences DataFrame for the signed adjacency matrix associated with the GRN
def build_influences(gene_list, pregraph=None, expr_mat_file=None, interactions=None, full_graph=False, tau=0., add_other_genes=False, direct_graph=False, cor_method="pearson", verbose=1, beta=1):
    N1,N = [len(gene_list)]*2
    assert len(interactions[0])==len(interactions[1])
    if (full_graph or (str(interactions)=="None")):
        pe = np.ones((N,N))
    else:
        gene_a, gene_b = interactions
        if (pregraph is not None):
            gene_aa = gene_a+[pregraph.index[i] for i, j in np.argwhere(pregraph.values!=0).tolist()]
            gene_bb = gene_b+[pregraph.columns[j] for i,j in np.argwhere(pregraph.values!=0).tolist()]
            gene_interactions = list(set(["--".join([gene_aa[i],gene_bb[i]]) for i in range(len(gene_aa))]))
            gene_aa = [g.split("--")[0] for g in gene_interactions]
            gene_bb = [g.split("--")[1] for g in gene_interactions]
        else:
            gene_aa,gene_bb = gene_a,gene_b
        gene_list = list(set(gene_aa+gene_bb))
        N = len(gene_list)
        pe = np.zeros((N,N))
        gene_a_id, gene_b_id = [gene_list.index(x) for x in gene_aa], [gene_list.index(x) for x in gene_bb]
        pe[gene_a_id, gene_b_id] = 2
        pe[gene_b_id, gene_a_id] = 2
        ## x2 for direction, x2 for sign
        #assert np.sum(pe)==4*len(interactions[0]) #valid if pregraph is None
    if (verbose):
        print("#edges in PPI: %d+%d=%d (x2 for direction, x2 for sign, additional unique edges, genes %d/%d)" % (len(interactions[0])*4, int(np.sum(pe))-4*len(interactions[0]), int(np.sum(pe)), N, N1))
    interactions = [gene_aa,gene_bb]
    # Correlation matrix
    mat = pd.read_csv(expr_mat_file, sep=",", header=0, index_col=0).iloc[:-3,:].apply(pd.to_numeric)
    mat = quantile_normalize(mat)
    mat = mat.loc[[g for g in mat.index if (g in gene_list)]]
    coexpr = np.power(mat.T.corr(method=cor_method), beta)
    N_ = len(coexpr.index)
    if (direct_graph):
        # CBDN: partial correlation: pcor(i,j) is the partial correlation between genes i and all genes k != j,i w.r.t. gene j
        coexpr_ = coexpr.values
        def pcor(i, j):
            num = (coexpr_[i,:]-np.multiply(coexpr_[i,:], coexpr_[j,:]))
            denom = np.power((1-np.power(coexpr_[i,:],2))*(1-np.power(coexpr_[j,:],2)), 0.5)
            denom[i] = denom[j] = 1
            pc = np.divide(num, denom)
            pc[i] = pc[j] = 0
            return pc
        def influence(j,i):
            return 1/(N_-1)*np.linalg.norm(coexpr_[i,:]-pcor(i,j), 1)
        arr_ids = np.matrix([[e1,e2] for e1,e2 in product(range(N_), range(N_))])
        direct_mat = np.array(np.vectorize(influence)(arr_ids[:,0], arr_ids[:,1])).reshape((N_,N_))
        # no two-gene cycles
        direct_mat[(direct_mat-direct_mat.T)<0] = 0
        direct_mat = np.array(direct_mat>0,dtype=int)
    else:
        direct_mat = np.ones((N_,N_))
    def build_influ():
        corr_mat = coexpr.values
        ids = [gene_list.index(g) for g in coexpr.index]
        infl_mat = np.multiply(np.multiply(corr_mat, np.asarray(pe[ids, :][:,ids]>0, dtype=bool)), direct_mat) 
        if (pregraph is not None):
            for i,j in np.argwhere(pregraph.values!=0).tolist():
                icol,jcol = [list(coexpr.index).index(v) for v in [pregraph.index[i], pregraph.columns[j]]]
                v = pregraph.values[i,j]
                if (v != 2):
                    infl_mat[icol][jcol] = pregraph.values[i,j]
                    infl_mat[jcol][icol] = pregraph.values[j,i]
                else:
                    infl_mat[icol][jcol] = corr_mat[icol][jcol]
                    infl_mat[jcol][icol] = corr_mat[jcol][icol]
        df = pd.DataFrame(infl_mat, index=coexpr.index, columns=coexpr.index)
        notisolated_genes = list(set([g for g in df.index if (df.loc[g].abs().sum()>0)]+[g for g in df.columns if (df[g].abs().sum()>0)]))
        ## remove isolated nodes
        df = df.loc[notisolated_genes][notisolated_genes]
        return df
    def ppi2edges(ppi_input):
        edges_output = list(zip(list(ppi_input["regulator"]),list(ppi_input["target"])))
        return edges_output
    infl_mat = build_influ()
    genes = list(set(list(infl_mat.index)+list(infl_mat.columns)))
    vals = [[list(infl_mat.index)[i], list(infl_mat.columns)[j], infl_mat.values[i,j]] for i, j in np.argwhere(infl_mat.values!=0).tolist()]
    ppi = pd.DataFrame(vals, index=range(len(vals)), columns=["regulator","target","score"]).sort_values(by="score", ascending=False)
    Ngenes = len(genes)
    edges = ppi2edges(ppi)
    components = get_weakly_connected(edges, genes)
    print("#edges = %d, #genes = %d/%d (non isolated), #components = %d" % (len(edges),Ngenes,N1,len(components)))
    ## Keep edges with value > tau threshold
    ppi_ = ppi.loc[ppi["score"].abs()>=tau]
    edges_ = ppi2edges(ppi_)
    genes_ = list(set([g for e in edges_ for g in e]))
    components = get_weakly_connected(edges_, genes_)
    print("#edges (corr>=tau=%.2f) = %d, #genes = %d, #components = %d" % (tau, len(edges_), len(genes_), len(components)))
    main_component_id = np.argmax([len(c) for c in components])
    ## add genes not in subset + not from the main (biggest) component
    isolated_genes = [g for g in genes if (g not in genes_)]+[g for cid, c in enumerate(components) if (cid!=main_component_id) for g in c]
    print("tau=%s: #isolated genes = %d" % (tau,len(isolated_genes)))
    assert len(components[main_component_id])+len(isolated_genes)==Ngenes
    edges = [e for e in edges_]
    ## Add edges in decreasing score order until the PPI is connected and all genes are present
    if (len(isolated_genes)>0):
        ppi_ = ppi.loc[ppi["score"].abs()<tau]
        if (ppi_.shape[0]>0):
            ids_isolated = [idx for idx in ppi_.index if ((ppi_.loc[idx]["regulator"] in isolated_genes) or (ppi_.loc[idx]["target"] in isolated_genes))]
            ppi_ = ppi_.loc[ids_isolated]
            ppi__scores = list(sorted(list(set(list(ppi_["score"].abs()))), reverse=True))
            i, done = 0, False
            edges__ = ppi2edges(ppi_)
            while (not done and i < len(ppi__scores)):
                edges__ = ppi2edges(ppi_.loc[ppi_["score"].abs()==ppi__scores[i]])
                edges += edges__
                genes_ = list(set([g for e in edges for g in e]))
                components = get_weakly_connected(edges, genes_)
                assert len(genes_)==sum([len(c) for c in components])
                print(" - #edges (corr>=%.5f) = %d, #genes = %d, #components = %d" % (ppi__scores[i], len(edges), len(genes_), len(components)))
                done = (len(components)==1) and (len(components[0])==Ngenes)
                i += 1
    assert len(components)==1 and len(components[0])==Ngenes
    vals = [[g1,g2,float(infl_mat.loc[g1][g2])] for g1,g2 in edges]
    infl_mat = pd.DataFrame(vals, index=range(len(vals)), columns=["regulator","target","score"]).pivot(index="regulator",columns="target",values="score").fillna(0)
    infl_mat[infl_mat<0] = -1
    infl_mat[infl_mat>0] = 1
    assert np.sum(np.abs(infl_mat.values))==len(edges)
    print("#edges %d\t#genes %d" % (np.sum(np.abs(infl_mat.values)), Ngenes))
    other_genes = []
    if (add_other_genes):
        ## adding genes not in coexpr
        other_genes = [g for g in gene_list if (g not in genes)]
        infl_mat = np.hstack((2*np.ones((N_,N-N_)), infl_mat))
        infl_mat = np.vstack((2*np.ones((N-N_,N)), infl_mat))
    influences = infl_mat
    return influences

#' @param influences DataFrame from build_influences
#' @param maxclause integer: small means risking not finding any solution, large means inference is time-consuming
#' @param exact boolean: takes into account all edges in inference
#' @param verbose prints stuff
#' @return InfluenceGraph from bonesis
#' @param exact
def create_grn(influences, maxclause=None, exact=False, verbose=1):
    if (str(maxclause)=="None"):
        maxclause = get_maxdegree(influences)
    influences_list = []
    ## x -> y
    for nx,ny in np.argwhere(influences.values!=0).tolist():
        x,y,v = list(influences.index)[nx], list(influences.columns)[ny], influences.values[nx,ny]
        assert v in [2,1,-1]
        if (v==2):
            signs = [1,-1]
        else:
            signs = [v]
        for sign in signs:
            influences_list.append((x,y,dict(sign=sign)))
    grn = bonesis.InfluenceGraph(graph=influences_list, maxclause=maxclause, exact=exact)
    return grn

## In-degree
def get_maxdegree(influences):
    maxindegree=int(np.max(np.sum(influences.values, axis=1)))
    return maxindegree

#' @param grn InfluenceGraph from bonesis
#' @param signatures DataFrame from binarize_experiments
#' @param ChromHMM_sig cell-specific mask for profiles
#' @param exps_ids subset of experiments
#' @param attractor_states list of attractor states
#' @return BoNesis object from bonesis
def build_observations(grn, signatures, disease_cell_line, ChromHMM_sig=None, exps_ids=[], attractor_states=[]):
    exps = [x for x in signatures.columns if ("initial" not in x)]
    if (len(exps_ids)==0):
        exps_ids = range(len(exps))
    print("\t".join([exps[i] for i in exps_ids]))
    if (str(ChromHMM_sig)!="None"):
        cell_lines = list(ChromHMM_sig.columns)
    else:
        cell_lines = []
    ## Instantiate experiments
    ## zero state
    data_exps = {"zero": {g: 0 for g in list(set(list(signatures.index)+list(grn.nodes)))}}
    for exp_nb in exps_ids:
        cell = exps[exp_nb].split("_")[-1]
        cols = ["Exp%d_"%(exp_nb+1)+x for x in ["init", "final"]]
        data_df = signatures[["initial_"+cell, exps[exp_nb-1]]]#.dropna()
        data_df = data_df.T
        ## Compatible with perturbation experiment
        pert, sign = exps[exp_nb].split("_")[0], 0 if (exps[exp_nb].split("_")[1]=="KD") else 1
        data_df[pert] = sign 
        data_df = data_df.T
        data_df.columns = cols
        for col in cols:
            data_exps.update(data_df[[col]].dropna().astype(int).to_dict())
    ## Instantiate attractor states (epilepsy)
    for state_id, state in enumerate(attractor_states):
        data_df = pd.DataFrame(state.values, index=state.index, columns=state.columns)
        data_df.columns = ["attractor_%d" % (state_id+1)]
        for col in data_df.columns:
            data_exps.update(data_df[[col]].dropna().astype(int).to_dict())
    print("Experiments")
    print_exps = pd.DataFrame.from_dict(data_exps, orient="index").fillna(-1).astype(int)
    print_exps[print_exps==-1] = ""
    print(print_exps)
    BO = bonesis.BoNesis(grn, data_exps)
    ## Instantiate external perturbagens
    ## reachability & fixed point constraints
    for exp_nb in exps_ids:
        exp = exps[exp_nb]
        pert, sign = exp.split("_")[0], 0 if (exp.split("_")[1]=="KD") else 1
        cell = exp.split("_")[-1]
        mutant_e = {}
        if (cell in cell_lines and str(ChromHMM_sig)!="None"):
            di = ChromHMM_sig[[cell]].dropna().astype(int).to_dict()[cell]
        else:
            di = {}
        mutant_e.update(di)
        mutant_e.update(dict([[pert, sign]]))
        with BO.mutant(mutant_e) as m:
            final_m = m.fixed(m.obs("Exp%d_final" % (exp_nb+1)))
            ## Trajectory : init -> trap space
            ~m.obs("Exp%d_init" % (exp_nb+1)) >= final_m #~m.obs("Exp%d_final" % (exp_nb+1))
            ## Universal reachable fixed points
            ~m.obs("Exp%d_init" % (exp_nb+1)) >> "fixpoints" ^ {m.obs("Exp%d_final" % (exp_nb+1))}
        ## No trivial trajectory from zero state
        di = data_exps["Exp%d_init" % (exp_nb+1)]
        di.update(mutant_e)
        if (all([di.get(g,1)==0 for g in grn.nodes])):
            continue
        ## any attractor state from this init should match with the considered fixed point in the mutant GRN
        ~BO.obs("zero") / ~BO.obs("Exp%d_final" % (exp_nb+1))
    ## Fixed point
    for state_id, state in enumerate(attractor_states):
        if (str(ChromHMM_sig)!="None"):
            cell_mask = ChromHMM_sig[[disease_cell_line]].dropna().astype(int).to_dict()[disease_cell_line]
            with BO.mutant(cell_mask) as m:
                # fixed point
                #m.fixed(~m.obs("attractor_%d" % (state_id+1)))
                # trap space
                fm = m.fixed(m.obs("attractor_%d" % (state_id+1)))
        else:
            # fixed point
            #BO.fixed(~BO.obs("attractor_%d" % (state_id+1)))
            # trap space
            fmm = BO.fixed(BO.obs("attractor_%d" % (state_id+1)))
    BO.boolean_networks().standalone(output_filename=file_folder+name+".asp")
    return BO

def solution2influences(sol):
    N = sol.shape[0]
    nodes = list(sol.index)
    influences = np.zeros((N,N))
    grfs = get_grfs_from_solution(sol.T)
    for k in grfs:
        ku = nodes.index(k)
        for r in grfs[k]:
            if (r in ["0","1"]):
                continue
            if (r in nodes):
                ru = nodes.index(r)
                influences[ru,ku] = grfs[k][r]
            else:
                influences = np.concatenate((influences, np.zeros((1,N))), axis=0)
                influences = np.concatenate((influences, np.zeros((N+1,1))), axis=1)
                influences[N,ku] = grfs[k][r]
                N += 1
                nodes.append(r)
    return pd.DataFrame(influences, index=nodes, columns=nodes)

def zip2df(fname):
    with ZipFile(fname, "r") as zip:
        zip.extractall()
    grn_fnames = glob("*.bnet")
    grns = []
    for grn_fname in grn_fnames:
        grn = load_grn(grn_fname)
        grns.append(grn)
    solutions = pd.DataFrame(grns)
    sb.call("rm -f *.bnet", shell=True)
    return solutions

def load_grn(fname):
    return mpbn.MPBooleanNetwork(fname)

#' @param R DataFrame columns = solutions, rows = GRFs
#' @return id of the solution with the fewest number of "trivial" GRFs (i.e., no regulators) and #of trivial grfs in this solution
def get_trivial_solution(R, dtype):
    assert dtype in ["in", "out", "all"]
    degrees = get_degrees(R, dtype=dtype)
    nb_zerodeg = np.sum(np.array(degrees.values==0, dtype=int), axis=0)
    sol_id = np.argmax(nb_zerodeg)
    return R.columns[sol_id], np.max(nb_zerodeg)

#' @param R DataFrame columns = solutions, rows = GRFs
#' @return degrees of all genes for each solution
def get_degrees(R, dtype="all", weighted=False):
    assert dtype in ["all", "in", "out"]
    degrees = np.zeros(R.shape)
    if (weighted):
        ## STRING combined score
        STRING_score = pd.read_csv(file_folder+"network_score.tsv", sep="\t")
        STRING_score = STRING_score[["preferredName_A","preferredName_B","score"]]
        STRING_score = STRING_score.drop_duplicates(keep="first")
        STRING_score.index = ["--".join(list(sorted(list(STRING_score.loc[x][["preferredName_A","preferredName_B"]])))) for x in STRING_score.index]
        STRING_score = STRING_score["score"]
    for ic, col in enumerate(R.columns):
        grfs = get_grfs_from_solution(R[col])
        for k in grfs:
            if (weighted):
                weights_k = [float(STRING_score.loc["--".join(list(sorted([r,k])))]) for r in grfs[k]]
            if (dtype in ["all","in"]):
                degrees[list(R.index).index(k), ic] += sum(weights_k) if (weighted) else len(grfs[k])
            if (dtype in ["all", "out"]):
                #for ir, r in enumerate(grfs[k]):
                #    degrees[list(R.index).index(r), ic] += weights_k[ir] if (weighted) else 1
                degrees[[list(R.index).index(r) for r in grfs[k]], ic] += np.array([1]*len(grfs[k]) if (weighted) else weights_k)
    return pd.DataFrame(degrees, index=R.index, columns=R.columns)

#' @param R DataFrame columns = solutions, rows = GRFs
#' @return id of the solution with the fewest number of edges
def get_minimal_edges(R, connected=False, maximal=False):
    minimal_id = None
    nb_edges = float("inf") if (not maximal) else -float("inf")
    for ic, col in enumerate(R.columns):
        infl = solution2influences(R[col])
        nedges = int(np.sum(np.abs(infl.values)))
        if (nedges>400):
            print(nedges)
        if (connected):
            influences = solution2influences(R[col])
            e1,e2 = np.nonzero(influences.values!=0)
            edges = list(zip(e1.tolist(),e2.tolist()))
            components = get_weakly_connected(edges, list(R.index))
            if (len(components)>1):
                continue
        if (((nedges < nb_edges) and not maximal) or ((nedges > nb_edges) and maximal)):
            nb_edges = nedges
            minimal_id = ic
    return R.columns[minimal_id], nb_edges

#' @param solution DataFrame with GRFs
#' @return grfs dictionary version of solution: keys are genes, values are dictionaries of regulators with sign (-1, 1)
def get_grfs_from_solution(solution):
    sol = solution.to_dict()
    grfs = {}
    for sk, k in enumerate(sol):
        grf = "".join("".join(" ".join(" ".join(str(sol[k]).split("|")).split("&")).split("(")).split(")"))
        regulators = {}
        if (str(grf) not in ["0", "1"]):
            regs = list(set(grf.split(" ")))
            for r in regs:
                if (r[0] == "!"):
                    regulators.setdefault(r[1:], -1)
                else:
                    regulators.setdefault(r, 1)
        grfs.setdefault(k, regulators)
    return grfs

#' @param solution DataFrame with GRFs
#' @param fname file name for solution
#' @param verbose print to stdout
def save_grn(solution, fname, sep="<-", verbose=True, max_show=4, write=True):
    sol = solution.to_dict()
    print_sol = []
    for sk, k in enumerate(sol):
        print_sol.append("%s %s %s" % (str(k), sep, str(sol[k])))
    if (verbose):
        print("\n".join(print_sol[:max_show]+["..." if (len(print_sol)>max_show) else ""]))
    print_sol = "\n".join(print_sol)
    if (write):
        with open(fname, "w+") as f:
            f.write(print_sol)
    return None

def save_solutions(bnetworks, fname, limit):
    with ZipFile(fname, "w") as bundle:
        n = 0
        for i, bn in enumerate(bnetworks):
            with bundle.open(f"bn{i}.bnet", "w") as fp:
                fp.write(bn.source().encode())
            n += 1
            if n == limit:
                break
    return n

#' @param BO BoNesis object from bonesis
#' @param fname file name 
#' @param limit maximum number of solutions to look for
#' @param verbose prints stuff
#' @return list of # solutions per iteration
def infer_network(BO, fname="solutions", limit=50, verbose=1,use_diverse=True,niterations=1):
    if (use_diverse):
        infer = lambda bo : bo.diverse_boolean_networks
    else:
        infer = lambda bo : bo.boolean_networks
    param_di = {"skip_supersets": True}
    if (str(limit)!="None"):
        param_di.setdefault("limit",limit)
    bnetworks = infer(BO)(**param_di)
    nsolutions = []
    for niter in tqdm(range(niterations)):
        if (not os.path.exists(fname+"_"+str(niter+1)+".zip")):
            nsolutions.append(save_solutions(bnetworks, fname+"_"+str(niter+1)+".zip", limit))
    return nsolutions

def predict(model, state=None):
    if (str(state)=="None"):
        attractors = model.attractors()
    else:
        attractors = model.attractors(reachable_from=state)
    return attractors

def simplify(solution,maxvar=18):
    from itertools import product
    from quine_mccluskey.qm import QuineMcCluskey as QM
    gene_di = {ig:g for ig, g in enumerate(solution.index)}
    formula_s = []
    for g in solution.index:
        sol = str(solution.loc[g])
        sol = " ".join(sol.split("("))
        sol = " ".join(sol.split(")"))
        sol = " ".join(sol.split("!"))
        sol = " ".join(sol.split("|"))
        sol = " ".join(sol.split("&"))
        gene_list = list(set([x for x in sol.split(" ") if (len(x)>0)]))
        if (all([x in ["0","1"] for x in gene_list])):
            if ("1" in gene_list and "0" not in gene_list):
                terms_s = "1"
            else:
                terms_s = "0"
        elif (len(gene_list)>maxvar):
            terms_s = str(solution.loc[g])
        else:
            gene_di = {ig: g for ig, g in enumerate(gene_list)}
            n = len(gene_di)
            values = list(product([0,1], repeat=n))
            grf_g = solution.loc[g]
            grf_g = " or ".join(grf_g.split("|"))
            grf_g = " not ".join(grf_g.split("!"))
            grf_g = " and ".join(grf_g.split("&"))
            def eval_g(x, grf):
                evaluated_grf = "".join([s for s in grf])
                for ix, xx in enumerate(x):
                    evaluated_grf = str(xx>0).join(evaluated_grf.split(gene_di[ix]))
                return eval(evaluated_grf)
            terms = ["".join([str(vv) for vv in v]) for v in values if (eval_g(v, grf_g))]
            terms_qm = QM().simplify_los(terms)
            if (str(terms_qm)=="None"):
                terms_qm = []
            terms_s = "("+(")|(".join(["&".join([("" if (t=="1") else "!")+gene_di[it] for it, t in enumerate(term) if (t != "-")]) for term in terms_qm]))+")"
            if (len(terms_s)==0):
                terms_s = "1" if (len(terms)>0) else "0"
        formula_s.append(terms_s)
    simplified_solution = pd.Series(formula_s, index=solution.index)
    return simplified_solution

##https://doi.org/10.5281/zenodo.3938904 "Tumour - Synthesis with BoNesis"
def has_cyclic(bn):
    mbn = mpbn.MPBooleanNetwork(bn)
    for a in mbn.attractors():
        if "*" in a.values():
            return True
    return False
