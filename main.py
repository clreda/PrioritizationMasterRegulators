# coding:utf-8

import sys
import requests
import json
import pandas as pd
import subprocess as sb
import os
import numpy as np
import pickle
from glob import glob

from itertools import chain, combinations
#https://docs.python.org/3/library/itertools.html#itertools-recipes
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

from params import *
from utils_exp import get_experimental_constraints, plot_signatures, plot_distributions, plot_discrete_distributions, profiles2signatures
from utils_grn import build_influences, create_grn, build_observations, infer_network, influences2graph, get_genes_interactions_from_PPI, symmetrize_influences

sys.path.insert(1, "utils/")
from LINCS_utils import *
from STRING_utils import *
from DISGENET_utils import *

np.random.seed(seed_number)
from random import seed
seed(seed_number)

for folder in [file_folder, plot_folder, solution_folder]:
    sb.call("mkdir -p "+folder, shell=True)

with open(solution_folder+"params.json", "w+") as f:
    f.write(str(args))

fnames = []
## remove previous solutions
#fnames += glob(root_folder+"*solutions-*")
## remove previous interaction constraints
#fnames += glob(root_folder+"*influences.*")
## remove previous profiles
#fnames += glob(root_folder+"*profiles.csv")
## remove previous signatures
#fnames += glob(root_folder+"*signatures.csv")
for fname in fnames:
    sb.call("rm "+fname, shell=True)

#<><><><><>
# Genes 
#<><><><><>

if (str(path_to_genes)=="None"):
    ## See DisGeNet API documentation for the definitions of the parameters 
    gene_df = get_genes_proteins_from_DISGENET([disease_cid], min_score=disgenet_args["min_score"], min_ei=disgenet_args["min_ei"], min_dsi=disgenet_args["min_dsi"], min_dpi=disgenet_args["min_dpi"])
    model_genes = list(set(list(gene_df["Gene Name"])[0].split("; ")))
    with open(file_folder+"model_gene_set.txt", "w") as f:
        f.write(("\n".join(model_genes))+"\n")
else:
    ## Import preselected set of genes of interest
    with open(path_to_genes, "r") as f:
        model_genes = f.read().split("\n")[:-1]

print("#genes in model %d" % len(model_genes))

fname_network = file_folder+"network_score="+str(score_STRING)+".tsv"
if (not os.path.exists(fname_network.split("=")[0]+".tsv")):
    network = get_network_from_STRING(model_genes, score=0, taxon_id=taxon_id, fname=fname_network.split("=")[0]+".tsv", quiet=True)

## Gridsearch for score_STRING
if (not os.path.exists(file_folder+"/score_STRING_results.csv")):
    build_idxs = lambda ls1, ls2 : ["-".join(list(sorted(list(e)))) for e in list(zip(list(ls1), list(ls2)))]
    score_range = list(range(1000,95,-5))
    score_results = pd.DataFrame([], index=score_range, columns=["#edges"])
    network_fname = fname_network.split("=")[0]+".tsv"
    ground_truth_fname = file_folder+"network_score.tsv"
    assert os.path.exists(ground_truth_fname)
    ground_truth = pd.read_csv(ground_truth_fname, sep="\t", index_col=0)[["preferredName_A", "preferredName_B"]]
    ground_truth.columns = ["gene1", "gene2"]
    binding = build_idxs(list(ground_truth["gene1"]), list(ground_truth["gene2"]))
    overlap = [0]*len(score_range)
    for iss,s in enumerate(score_range):
        print("Score = %d (%d/%d)" % (s, iss+1, len(score_range)))
        genes, interactions = get_genes_interactions_from_PPI(network_fname, connected=True, score_STRING=s)
        nedges = len(interactions[0])
        score_results.loc[s] = nedges
        intersect = build_idxs(interactions[0],interactions[1])
        overlap[iss] = len(list(set(intersect).intersection(set(binding))))*100/len(set(binding))
    network = pd.read_csv(network_fname, index_col=0, sep="\t")
    network.index = build_idxs(list(network["preferredName_A"]), list(network["preferredName_B"]))
    network = network.sort_values(by="score", ascending=False)
    network = network[~network.index.duplicated(keep="first")] #highest score is kept for every undirected pair
    score_results["#edges score>=thres"] = [network.loc[network["score"]>=s/1000.].shape[0] for s in score_range]
    score_results["overlap physical (binding)"] = overlap
    score_results.to_csv(file_folder+"score_STRING_results.csv")
    print(score_results)
    score_STRING = score_range[np.argmin(list(score_results["#edges"]))]
    print("Select score_STRING = "+str(score_STRING))

if (not os.path.exists(fname_network)):
    model_genes_, interactions = get_genes_interactions_from_PPI(fname_network.split("=")[0]+".tsv", connected=True, score_STRING=score_STRING)
    network_df = pd.DataFrame(interactions, index=["nodeA","nodeB"]).T
    network_df.to_csv(fname_network, sep="\t",header=None,index=None)
network_df = pd.read_csv(fname_network,header=None,sep="\t")
interactions = [e[0] for e in network_df.values.tolist()], [e[1] for e in network_df.values.tolist()]
model_genes_ = list(set(interactions[0]+interactions[1]))

model_genes = model_genes_
Ngenes = len(model_genes)

print("score_STRING %d\t#genes (non isolated in PPI) %d\t#edges in PPI %d" % (score_STRING, Ngenes, len(interactions[0])))

if (plot_it):
    influences = np.zeros((Ngenes, Ngenes))
    for i in range(len(interactions[0])):
        r, t = model_genes.index(interactions[0][i]), model_genes.index(interactions[1][i])
        influences[r,t] = influences[t,r] = 2
    influences = pd.DataFrame(influences, index=model_genes, columns=model_genes)
    influences2graph(influences, plot_folder+"ppi_string="+str(score_STRING), optional=True)

#<><><><><>
# Experimental constraints
#<><><><><>

# Get experiments from LINCS (Level 3)

## Convert gene names into EntrezGene CID
entrezgene_fname=file_folder+"entrezgenes_ids.csv"
if (not os.path.exists(entrezgene_fname)):
    from utils_data import request_biodbnet
    probes = request_biodbnet(model_genes, from_="Gene Symbol and Synonyms", to_="Gene ID", taxonId=taxon_id, chunksize=100)
    other_ids = list(probes[probes["Gene ID"]=="-"].index)
    probes = probes[probes["Gene ID"]!="-"]
    df_list = [probes]
    if (len(other_ids)>0):
        other_probes = request_biodbnet(other_ids, from_="Ensembl Gene ID", to_="Gene ID", taxonId=taxon_id, chunksize=100)
        other_ids = list(other_probes[other_probes["Gene ID"]=="-"].index)
        other_probes = other_probes[other_probes["Gene ID"]!="-"]
        df_list.append(other_probes)
    if (len(other_ids)>0):
        other_other_probes = request_biodbnet(other_ids, from_="HGNC ID", to_="Gene ID", taxonId=taxon_id, chunksize=100)
        other_ids = list(other_other_probes[other_other_probes["Gene ID"]=="-"].index)
        other_other_probes = other_other_probes[other_other_probes["Gene ID"]!="-"]
        df_list.append(other_other_probes)
    if (len(other_ids)>0):
        output_format = "tsv"
        method = "get_string_ids"
        params = {
            "identifiers" : "\r".join(other_ids), # your protein list
            "species" : taxon_id, # species NCBI identifier
            "limit" : 1, # only one (best) identifier per input protein
            "echo_query" : 1, # see your input identifiers in the output
            "caller_identity" : app_name # your app name
        }
        request_url = "/".join([string_api_url, output_format, method])
        results = requests.post(request_url, data=params).text
        sleep(1)
        from io import StringIO
        results = pd.read_csv(StringIO(results), sep="\t")
        other_ids = []
        for idx in list(results["queryItem"]):
            if (idx == "ENSP00000451560"):
                other_ids.append("TPPP2")
            else:
                other_ids.append(idx)
        results = request_biodbnet(other_ids, from_="Gene Symbol and Synonyms", to_="Gene ID", taxonId=taxon_id, chunksize=100)
        df_list.append(results)
    if (len(df_list)>0):
        probes = pd.concat(tuple(df_list), axis=0)
    probes.to_csv(entrezgene_fname)
probes = pd.read_csv(entrezgene_fname,index_col=0)

if (len(list(probes[probes["Gene ID"]=="-"].index))>0):
    print("Not found genes:")
    print(list(probes[probes["Gene ID"]=="-"].index))
probes = probes[probes["Gene ID"]!="-"]
print("#genes available %d (can get their EntrezID: if lower than the previous number %d, check that all input genes are gene symbols or Ensembl IDs or HGNC IDs or STRING IDs, and check you gave the correct taxon id)" % (probes.shape[0], len(model_genes)))
model_genes = list(probes.index)

## convert EntrezGene ids back to gene symbols
if (not os.path.exists(file_folder+"entrezid2symbols.csv")):
    pert_inames = [None]*len(model_genes)
    entrez_ids = [None]*len(model_genes)
    for ig, g in enumerate(model_genes):
        endpoint="genes"
        method="filter"
        all_entrezid = str(probes.loc[g]["Gene ID"]).split("; ")
        for entrezid in all_entrezid:
            params = {"where":{"entrez_id": str(entrezid)},"fields":["gene_symbol"]}
            request_url = build_url(lincs_api_url, endpoint, method, params=params, key=user_key)
            data = post_request(request_url, quiet=True, stime=0)
            if (len(data)==0):
                continue
            else:
                pert_inames[ig] = data[0]["gene_symbol"]
                entrez_ids[ig] = entrezid
                if (pert_inames[ig]==g):
                    break
        print("\t".join([g, pert_inames[ig], str(ig+1), str(len(model_genes))]))
    pert_iname_ids = [ip for ip,p in enumerate(pert_inames) if (str(p)!="None")]
    pd.DataFrame([[pert_inames[i] for i in pert_iname_ids], [entrez_ids[i] for i in pert_iname_ids]], columns=[model_genes[i] for i in pert_iname_ids], index=["Gene Symbol", "Entrez ID"]).T.to_csv(file_folder+"entrezid2symbols.csv")
pert_df = pd.read_csv(file_folder+"entrezid2symbols.csv", index_col=0)
pert_inames = list(pert_df["Gene Symbol"])
entrez_ids = list(pert_df["Entrez ID"])
model_genes = list(pert_df.index)

## Get all cell lines
fname_lineages = file_folder+"lincs_lines_lineages.pck"
if (not os.path.exists(fname_lineages)):
    if (len(cell_lines)==0):
        ## Select cell lines where M30 genes are perturbed
        endpoint = "sigs"
        method = "filter"
        cell_lines = []
        for ig, g in enumerate(pert_inames):
            print("Gene %s %d/%d" % (g, ig+1, len(model_genes)))
            params = {"where": {"pert_iname": g}, "fields": ["cell_id"]}
            request_url = build_url(lincs_api_url, endpoint, method, params=params, key=user_key)
            response = requests.get(request_url)
            assert response.status_code == 200
            data = json.loads(response.text)
            cell_lines_gene = list(set([d["cell_id"] for d in data]))
            print("#unique cells %d" % len(cell_lines_gene))
            cell_lines = list(set(cell_lines_gene+cell_lines))
    lincs_lines_lineages = {}
    for ic, cell in enumerate(cell_lines):
        print("Cell %d/%d" % (ic+1, len(cell_lines)))
        endpoint = "cells"
        method = "filter"
        params = {"where": {"cell_iname": cell}, "fields": ["cell_lineage"]}
        request_url = build_url(lincs_api_url, endpoint, method, params=params, key=user_key)
        response = requests.get(request_url)
        if (response.status_code != 200):
            continue
        data = json.loads(response.text)
        lincs_lines_lineages.setdefault(cell, [x for d in data for x in d["cell_lineage"]])
    with open(fname_lineages, "wb") as f:
        pickle.dump(lincs_lines_lineages, f)
with open(fname_lineages, "rb") as f:
        lincs_lines_lineages = pickle.load(f)

cell_lines = [x for x in lincs_lines_lineages]
matchings_ChromHMM_lincs = {c: c for c in cell_lines}

lineages_cell = { cell : {'LINCS': lincs_lines_lineages[cell], 'ChromHMM': matchings_ChromHMM_lincs[cell]} for cell in cell_lines}
pd.DataFrame(lineages_cell).to_csv(file_folder+"single-gene_"+"-".join(cell_lines[:4])+"_lineages.csv")

print("#genes = %d (no change expected)" % len(model_genes))

## Experimental profiles
profiles_fname = file_folder+"single-gene_"+"-".join(cell_lines[:4])+"_profiles_thresiscale"+str(thres_iscale)+".csv"
if (not os.path.exists(profiles_fname)):
    profiles = get_experimental_constraints(cell_lines, pert_types, entrez_ids, pert_inames, selection, thres_iscale=thres_iscale, nsigs=2) #genes
    profiles.to_csv(profiles_fname)
profiles = pd.read_csv(profiles_fname, header=0, index_col=0)

# Discard genes from subset which expression is not reported in LINCS
gene_tf_list = [list(pert_df.index)[list(pert_df["Entrez ID"]).index(int(g))] for g in list(profiles.iloc[:-3,:].index)]
a, b = interactions
new_a, new_b = [], []
for i in range(len(a)):
    if ((a[i] in gene_tf_list) and (b[i] in gene_tf_list)):
        new_a.append(a[i])
        new_b.append(b[i])
interactions = [new_a, new_b]
print("#genes = %d (present in LINCS L1000)" %len(gene_tf_list))
print("#interactions = %d (with both genes in LINCS L1000)" % len(interactions[0]))

#<><><><><>
# Interaction constraints
#<><><><><>

if (plot_it):
    Ngenes = len(gene_tf_list)
    influences = np.zeros((Ngenes, Ngenes))
    for i in range(len(interactions[0])):
        r, t = gene_tf_list.index(interactions[0][i]), gene_tf_list.index(interactions[1][i])
        influences[r,t] = influences[t,r] = 2
    influences = pd.DataFrame(influences, index=gene_tf_list, columns=gene_tf_list)
    influences2graph(influences, plot_folder+"ppi_preinfl_string="+str(score_STRING), optional=True)

influences_fname = file_folder+"single-gene_"+"-".join(cell_lines[:4])+"_influences_tau="+str(tau)+"_beta="+str(beta)+".csv"
if (not os.path.exists(influences_fname)):
    gene_tf_list_ = [str(g) for g in list(profiles.iloc[:-3,:].index)]
    interactions_ = [[str(pert_df.loc[g]["Entrez ID"]) for g in ls] for ls in interactions]
    grn_influences = build_influences(gene_tf_list_, pregraph=None, expr_mat_file=profiles_fname, interactions=interactions_, full_graph=full_graph, add_other_genes=add_other_genes, direct_graph=direct_graph, cor_method=cor_method, beta=beta, tau=tau) #genes
    print("#interactions (preprocessing, trimming out isolated nodes): %d" % np.sum(np.abs(grn_influences.values)))
    grn_influences = symmetrize_influences(grn_influences).fillna(0)
    ## use gene symbols
    grn_influences.index = [list(pert_df.index)[list(pert_df["Entrez ID"]).index(int(g))] for g in grn_influences.index]
    grn_influences.columns = [list(pert_df.index)[list(pert_df["Entrez ID"]).index(int(g))] for g in grn_influences.columns]
    influences = grn_influences
    influences.to_csv(influences_fname)
influences = pd.read_csv(influences_fname, index_col=0, header=0)
genes = list(influences.index) #
nedges = np.sum(np.abs(influences.values))
print("#filtered interactions = %d, non isolated LINCS L1000 genes %d" % (nedges, len(genes)))
assert nedges > 0

# Print graph
if (plot_it):
    influences2graph(influences, plot_folder+influences_fname.split(".csv")[0].split("/")[-1], optional=True, compile2png=True, engine="sfdp")

#<><><><><>
# Signatures
#<><><><><>

# Binarize experiments
sigs_fname = file_folder+"single-gene_"+"-".join(cell_lines[:4])+"_signatures_binthres="+str(bin_thres)+".csv"
if (not os.path.exists(sigs_fname)):
    signatures = profiles2signatures(profiles, save_fname=file_folder+"data_profiles", thres=bin_thres, selection=selection, backgroundfile=True, bin_method=bin_method if (not use_CD) else bin_method+"_CD")
    signatures.to_csv(sigs_fname)
signatures = pd.read_csv(sigs_fname, index_col=0, header=0)
signatures = signatures.dropna(how="all")
signatures.index = [list(pert_df.index)[list(pert_df["Entrez ID"]).index(idx)] for idx in signatures.index]
print("#experiments %d with %d cell lines on %d genes" % (len([c for c in signatures.columns if ("initial" not in c)]), len([c for c in signatures.columns if ("initial" in c)]), signatures.shape[0]))

#<><><><><>
# Remove genes which do not appear in signatures
#<><><><><>

gene_list = list(set([g for g in genes if (g in signatures.index)]))
signatures = signatures.loc[[g for g in gene_list]] 

from copy import deepcopy
signatures_copy = deepcopy(signatures)
signatures_copy[signatures_copy==0] = -1
signatures_copy = signatures_copy.fillna(0)
print("Frobenius norm signature matrix: %f" % np.linalg.norm(signatures_copy.values))

# Visualization
if (plot_it):
    from utils_state import binarize_experiments
    plot_distributions(profiles, fname=plot_folder+"gene_expression_distribution.png")
    vals = profiles.iloc[:-3,:].apply(pd.to_numeric)
    norm_vals = binarize_experiments(vals, thres=bin_thres, return_bin=False)
    norm_vals = norm_vals[[c for c in norm_vals.columns if (c in profiles.columns)]]
    profiles_bin = pd.concat((norm_vals, profiles.iloc[-3:,:]), axis=0)
    plot_distributions(profiles_bin, fname=plot_folder+"gene_expression_normalized_distribution.png", thres=bin_thres)
    plot_discrete_distributions(signatures, fname=plot_folder+"signature_expression_distribution.png")
    perturbed_genesymbols = list(set([i.split("_")[0] for i in profiles.loc["signame"]]))
    perturbed_entrezid = [str(pert_df.loc[g]["Entrez ID"]) for g in perturbed_genesymbols]
    plot_signatures(signatures, fname=plot_folder+"signatures_allgenes.png", max_show=len(gene_list))
    plot_signatures(signatures, fname=plot_folder+"signatures_50mostvariablegenes.png", max_show=50)

#<><><><><>
# Boolean Network
#<><><><><>

# Create, add exp. constraints, infer boolean network
solution_fname=root_folder+"/"+name+"/solutions-"+str(limit)+"_binthres="+str(bin_thres)+"_"+"-".join(cell_lines[:4])
fname_ls = glob(solution_fname+"_*.zip")

# If there is no known solution, perform the inference
if (len(fname_ls) != niterations):
    grn = create_grn(influences, exact=exact)
    gene_list = grn.nodes
    signatures = signatures[[c for c in signatures.columns if (("initial" in c) or any([(g in c) for g in gene_list]))]]
    signatures = signatures[list(sorted(list(signatures.columns)))]
    print("#experiments %d"%(len([c for c in signatures.columns if ("initial" not in c)])))
    total = len([x for x in signatures.columns if ("initial" not in x)])
    subsets = list(powerset(list(range(total))))[1:]
    subsets.reverse() # start with largest subset, removing the empty set
    for subset in subsets:
        subset = use_subset_exp if (use_subset_exp is not None) else subset
        print(subset)
        BO = build_observations(grn, signatures, disease_cell_line, exps_ids=subset, attractor_states=[], ChromHMM_sig=None)
        nsolutions = infer_network(BO, fname=solution_fname, limit=limit, verbose=1, use_diverse=use_diverse, niterations=niterations)
        if (sum(nsolutions)==0):
            print("No solution found. Try decreasing value bin_thres="+str(bin_thres)+" in [0,0.5].")
            sb.call("rm -f "+solution_fname+"*.zip", shell=True)
        else:
            with open(root_folder+"/"+name+"/subset-"+str(limit)+"_binthres="+str(bin_thres)+"_"+"-".join(cell_lines[:4])+".txt", "w") as f:
                f.write(str(subset))
        exit()
