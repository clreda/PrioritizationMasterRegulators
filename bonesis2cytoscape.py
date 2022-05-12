#coding: utf-8

import sys
import subprocess as sb
import numpy as np
import pandas as pd
import os
from glob import glob

from params import *

solution_fname = solution_folder+"solution"
sigs = glob(file_folder+"single-gene_"+"-".join(disease_cell_line)+"_signatures_binthres="+str(args["bin_thres"])+".csv")[0]
perturbed_genes = [[c.split("_")[0], "+" if (c.split("_")[1]=="OE") else "-"] for c in pd.read_csv(sigs, index_col=0).columns if ("initial" not in c)]
unique_perturbed_genes = list(set([g for g, s in perturbed_genes]))
perturbed_genes = dict([[g, "".join(list(sorted(list(set([x for gg, x in perturbed_genes if (gg==g)])))))] for g in unique_perturbed_genes])

results_folder=root_folder+name+"/cytoscape/"
sb.call("mkdir -p "+results_folder, shell=True)

sys.path.insert(1,"utils/")

from utils_grn import get_grfs_from_solution
from style_file import style_file

#############################################################
## CONVERT solution INTO CYTOSCAPE FILE+STYLE FILE IN .XML ##
#############################################################
## solution file is created after running model_selection.py
if (not os.path.exists(results_folder+solution_fname.split("/")[-1]+".sif")):
    ## Create SIF file
    target = []
    sol = pd.read_csv(solution_fname+".bnet", sep=" <- ", index_col=0, header=None, engine="python")
    gene_names = list(sol.index)
    grfs = get_grfs_from_solution(sol[sol.columns[0]])
    print(pd.DataFrame(grfs).fillna(0).astype(int))
    target = [" ".join([r, "->"+("+" if (grfs[g][r]>0) else "-"), g]) for g in grfs for r in grfs[g]]
    with open(results_folder+solution_fname.split("/")[-1]+".sif", "w+") as f:
        f.write("\n".join(target))

    ## Build associated Cytoscape style file
    target = []
    target += style_file[0]
    target.append("<visualProperty default=\""+solution_fname.split("/")[-1]+"\" name=\"SOLUTION\"/>")
    target += style_file[1]
    target += ["<visualProperty default=\"#FFCC99\" name=\"NODE_FILL_COLOR\">", "<discreteMapping attributeName=\"name\" attributeType=\"string\">"]
    colors = {"+-": "#f5eb2a", "+": "#3bd98a", "-": "#f77c39", "-+": "#f5eb2a"}
    target += ["<discreteMappingEntry attributeValue=\""+g+"\" value=\""+colors[perturbed_genes[g]].upper()+"\"/>" for g in gene_names if (g in perturbed_genes)]
    target += ["</discreteMapping>", "</visualProperty>"]

    target += style_file[2]

    with open(results_folder+solution_fname.split("/")[-1]+".xml", "w+") as f:
        f.write("\n".join(target))

#################################################
## CONVERT COLLAPSED MODEL INTO CYTOSCAPE FILE ##
#################################################
## collapsed model file is created after running model_selection.py
if (not os.path.exists(results_folder+"collapsed_model.tsv")):
    collapsed_model = pd.read_csv(root_folder+name+"/collapsed_model.csv", index_col=0)
    gene_names = list(collapsed_model.index)
    ## Create SIF file
    edges = [[collapsed_model.index[i],collapsed_model.columns[j],"->"+("-" if (collapsed_model.values[i,j]<0) else "+"), str(int(np.abs(collapsed_model.values[i,j])))] for i,j in np.argwhere(collapsed_model.values!=0).tolist()]
    edges_df = pd.DataFrame(edges, index=range(len(edges)), columns=["regulator","regulated","interaction","count"])
    edges_df.to_csv(results_folder+"collapsed_model.tsv", sep="\t")

#################################################
## ANNOTATE EDGES ACCORDING TO STRING DATABASE ##
#################################################
## Once the .sif file has been opened in Cytoscape and the edge table has been saved (with the default filename)

edge_fname = results_folder+"solution.sif default edge.csv"
string_fname = file_folder+"network_score.tsv"

if (not os.path.exists(edge_fname)):
    edges = pd.read_csv(edge_fname)
    edges_list = ["--".join(list(sorted([x.split(" ")[i] for i in [0,2]]))) for x in list(edges["name"])]
    edges.index = edges_list
    string = pd.read_csv(string_fname, sep="\t")[['preferredName_A','preferredName_B', 'score', 'nscore', 'fscore', 'pscore', 'ascore', 'escore', 'dscore', 'tscore']]
    string.index = ["--".join(list(sorted(list(string.loc[idx][["preferredName_A",'preferredName_B']])))) for idx in string.index]
    string = string[string.columns[2:]]
    string = string.loc[edges_list]
    string = string.loc[~string.index.duplicated()]
    assert len(list(set(edges_list))) == string.shape[0] ## remove direction sign
    ## https://string-db.org/help//api/#retrieving-similarity-scores-of-the-protein-set
    edges["interaction type"] = ["Interaction" if (string.loc[e]["escore"]>0) else ("Coexpression" if (string.loc[e]["ascore"]>0) else "Comention") for e in edges_list]
    edges.to_csv(edge_fname.split(".csv")[0]+"_annotated.csv")
