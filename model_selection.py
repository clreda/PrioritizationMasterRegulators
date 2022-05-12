#coding: utf-8

import sys
import requests
import json
import pandas as pd
import subprocess as sb
import os
import numpy as np
import pickle
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
from zipfile import ZipFile
from copy import deepcopy
from functools import reduce

from params import *
from utils_grn import get_degrees, get_minimal_edges, save_grn, influences2graph, solution2influences, zip2df, get_weakly_connected, general_topological_parameter

# Create, add exp. constraints, infer boolean network
root_folder += name+"/"
solution_fname=root_folder+"solutions-"
assert len(glob(solution_fname+"*_binthres="+str(bin_thres)+"_*.zip"))>0
cmd = "cd "+root_folder+" && file $(ls | grep '.zip') | grep 'Zip archive data' | cut -f1 -d':'"
fname_ls = sb.check_output(cmd,shell=True).decode("utf-8").split("\n")[:-1]
fname_ls = list(sorted([root_folder+x for x in fname_ls], key=lambda x: int(x.split(".zip")[0].split("_")[-1])))

if (len(fname_ls)>0):
    degrees_ls = []
    edges_ls = []
    sol_ls = []
    nsol = 1
    try:
        for fi, fname in enumerate(list(sorted(fname_ls))):
            solutions = zip2df(fname)
            nb_sol = solutions.shape[0]
            cols = list(range(nsol,nsol+nb_sol))
            nsol = np.max(cols)+1
            sol_ls.append(pd.DataFrame(solutions.T.values, index=solutions.columns, columns=cols))
            degrees = get_degrees(solutions.T)
            degrees_out = get_degrees(solutions.T, dtype="out")
            edges = pd.DataFrame(degrees.sum(axis=0), index=degrees.columns, columns=["It. #%d" % int(fname.split("_")[-1].split(".zip")[0])])
            edges_ls.append(edges)
            degrees = pd.DataFrame(degrees.mean(axis=1), index=degrees.index, columns=["It. #%d" % int(fname.split("_")[-1].split(".zip")[0])])
            degrees.columns = [int(fname.split(".zip")[0].split("_")[-1])]
            degrees_ls.append(degrees)
    except:
        print("'"+fname+"' not loaded.")
    nsol -= 1
    ## Plot average node degree distributions across iterations
    degrees_df = degrees_ls[0].join(degrees_ls[1:], how="outer")
    degrees_df = degrees_df[range(1,len(degrees_ls)+1)]
    if (len(degrees_ls)==50):
        plt.figure(figsize=(7,3)) #for the 50 solutions
    else:
        plt.figure(figsize=(4,3))
    degrees_df.boxplot(rot=90)
    plt.xlabel("Solution ID")
    plt.ylabel("Node degree")
    plt.savefig(plot_folder+"avg_degree_boxplot.png", bbox_inches="tight")
    plt.close()
    ## Plot number of edges across iterations
    edges_df = edges_ls[0].join(edges_ls[1:], how="outer")
    print(edges_df)
    edges_df.index = [""]
    plt.figure(figsize=(1,3))
    edges_df.T.boxplot(rot=90)
    plt.ylabel("Edge number")
    plt.savefig(plot_folder+"edge_number_boxplot.png", bbox_inches="tight")
    plt.close()
    ## All solutions
    sols = sol_ls[0].join(sol_ls[1:], how="outer")
    ## Number of unique solutions
    print("#unique solutions %d / %d" % (sols.T.drop_duplicates().shape[0], nsol))
    print(sols)
    ## Frequency of interactions
    perc_influences_df = sum([np.abs(solution2influences(sols[c])) for c in sols.columns])/nsol
    thres_perc=1
    rows_genes=perc_influences_df.sum(axis=1).loc[perc_influences_df.sum(axis=1)>thres_perc].index
    cols_genes=perc_influences_df.sum(axis=0).loc[perc_influences_df.sum(axis=0)>thres_perc].index
    influences_df_display = perc_influences_df.loc[rows_genes][cols_genes]
    influences_df_display = influences_df_display.loc[influences_df_display.sum(axis=1)>0]
    influences_df_display = influences_df_display.T.loc[influences_df_display.T.sum(axis=1)>0].T
    print(influences_df_display.shape)
    influences_df = deepcopy(perc_influences_df)
    influences_df["Regulator"] = influences_df.index
    redundant_values = pd.melt(influences_df, id_vars='Regulator', value_vars=[c for c in influences_df.columns if (c != "Regulator")]).sort_values(by="value", ascending=False)
    solutions = [solution2influences(sols[c]) for c in sols.columns]
    for thres_redundant in [0.75,1]:
        redundant_values_ = redundant_values.loc[redundant_values["value"]>=thres_redundant]
        ## in practice, should be always the same sign
        signs = ["".join(list(sorted(list(set(["+" if (s.loc[redundant_values_.loc[c]["Regulator"]][redundant_values_.loc[c]["variable"]]>0) else "-" for s in solutions]))))) for c in redundant_values_.index]
        redundant_values_["sign"] = signs
        print(redundant_values_[["Regulator","variable","sign"]])
        print(redundant_values_.shape)
    g = sns.clustermap(influences_df_display, square=True, figsize=(9,9), method="average", metric="correlation", row_cluster=False, col_cluster=False, cbar_pos=(0,0,0.02,0.8), cbar_kws={"ticks":[0,0.2,0.4,0.6,0.8,1]})
    ax = g.ax_heatmap
    ax.set_xticks(range(influences_df_display.shape[1]))
    ax.set_xticklabels(influences_df_display.columns, rotation=90, fontsize=10)
    ax.set_yticks(range(influences_df_display.shape[0]))
    ax.set_yticklabels(influences_df_display.index, fontsize=10)
    plt.savefig(plot_folder+"genewise_heatmap.png", bbox_inches="tight")
    plt.close()
    ## Histogram of different Boolean functions
    plt.figure(figsize=(5,25))
    nb_functions = pd.DataFrame([], columns=sols.index, index=["genes"])
    for g in sols.index:
        nb_functions[g] = len(list(set(map(str, list(sols.loc[g])))))
    nb_functions = nb_functions.T.sort_values(by="genes", ascending=False)
    im = nb_functions.boxplot()
    print(nb_functions)
    print(pd.DataFrame([[f(nb_functions.values) for f in [np.min, lambda x : np.percentile(x,25), lambda x : np.percentile(x,50), np.mean, lambda x : np.percentile(x,75), np.max]]], index=["#RFs"], columns=["min", "25perc", "median", "mean", "75perc", "max"]))
    t_unique = 7
    print("Genes with strictly more than %d unique GRFs" % t_unique)
    print(nb_functions.loc[nb_functions["genes"]>t_unique])
    for ig, g in enumerate(list(nb_functions.index)[:10]):
        im.text(1+(0.01 if (ig%2 == 0) else -0.06), int(nb_functions.loc[g]["genes"]), g, fontsize=20)
    im.set_ylabel("Number of unique reg. functions across solutions", fontsize=20)
    im.set_xlabel("", fontsize=20)
    plt.savefig(plot_folder+"unique_grfs_boxplot.png", bbox_inches="tight")
    plt.close()
    ## Number of unique solutions
    print("#unique solutions %d / %d" % (sols.T.drop_duplicates().shape[0], nsol))
    ## Model selection criterion
    ## * Minimal model (number of edeges) *
    minimal_id, minimal = get_minimal_edges(sols, connected=False)
    print("Minimal edge (not connected) solution (#edges %d)" % (minimal))
    influences_minimal = solution2influences(sols[minimal_id])
    influences2graph(influences_minimal, plot_folder+"inferred_minimal_solution", optional=False)
    ## * Maximal model (number of edges) *
    maximal_id, maximal = get_minimal_edges(sols, connected=False, maximal=True)
    print("Maximal edge (not connected) solution (#edges %d)" % (maximal))
    influences_maximal = solution2influences(sols[maximal_id])
    influences2graph(influences_maximal, plot_folder+"inferred_maximal_solution", optional=False)
    ## * Most connected model (minimal number of weakly connected components) *
    connectivity = float("inf") 
    id_sol = None
    for ic, ccol in enumerate(sols.columns):
        influences_ = solution2influences(sols[ccol])
        edges = [[influences_.index[u], influences_.columns[v]] for u,v in np.argwhere(influences_.values!=0).tolist()]
        genes = list(set(influences_.index).union(influences_.columns))
        components = get_weakly_connected(edges, genes)
        connection = len(components[np.argmax([len(c) for c in components])])
        if (connectivity > connection):
            connectivity = connection
            id_sol = ccol
    influences2graph(solution2influences(sols[id_sol]), plot_folder+"inferred_most_connected_solution", optional=False)
    ## * Maximizer of general topological criterion (GTP)
    Ds_ = [general_topological_parameter(solution2influences(sols[c]), perc_influences_df, weights) for c in sols.columns]
    plt.figure(figsize=(3,5))
    plt.boxplot(Ds_)
    print(pd.DataFrame([[f(Ds_) for f in [np.min, lambda x : np.percentile(x,25), lambda x : np.percentile(x,50), np.mean, lambda x : np.percentile(x,75), np.max]]], index=["value"], columns=["min", "25perc", "median", "mean", "75perc", "max"]))
    plt.ylabel("GTP")
    plt.xlabel("")
    plt.savefig(root_folder+"PLOTS/general_topological_parameter_solutions.png", bbox_inches="tight")
    plt.close()
    max_criterion_solution = sols[sols.columns[np.argmax(Ds_)]]
    max_criterion_influences = solution2influences(max_criterion_solution)
    influences2graph(max_criterion_influences, plot_folder+"inferred_max_criterion_solution", optional=False)
    ################################
    ## The final inferred network ##
    ################################
    save_grn(max_criterion_solution, root_folder+"solution.bnet")
    ## * Collapsed model (sign: sign of the interaction, absolute value: number of times the interactions appears) *
    collapsed_model = reduce(lambda x,y : x+y.loc[x.index][x.index], [solution2influences(sols[c]) for c in sols.columns])
    plt.figure(figsize=(3,5))
    plt.boxplot(collapsed_model.abs().values.flatten()*100/sols.shape[1])
    print(pd.DataFrame([[f(collapsed_model.abs().values.flatten()) for f in [np.min, lambda x : np.percentile(x,25), lambda x : np.percentile(x,50), np.mean, lambda x : np.percentile(x,75), np.max]]], index=["Edge count"], columns=["min", "25perc", "median", "mean", "75perc", "max"]))
    plt.ylabel("Frequency")
    plt.xlabel("")
    plt.savefig(root_folder+"PLOTS/frequency_edges.png", bbox_inches="tight")
    plt.close()
    collapsed_model.to_csv(root_folder+"collapsed_model.csv")
