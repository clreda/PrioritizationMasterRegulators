# coding: utf-8

import numpy as np
import pandas as pd
import subprocess as sb
import sys
import matplotlib.pyplot as plt
import requests
import json
import os

sys.path.insert(1, './utils/')
sys.path.insert(1, './utils/credentials/')

from LINCS_utils import *
from params import root_folder, name, file_folder
from utils_state import binarize_experiments

#' @param profiles
#' @return signatures
def profiles2signatures(profiles_, save_fname="profiles+LINCS", backgroundfile=False, selection="distil_ss", thres=0.5, non_constant_genes=None, bin_method="binary"):
    assert thres >= 0 and thres <= 0.5
    cell_lines = list(set([x.split("_")[1] for x in list(profiles_.loc["sigid"])]))
    signatures_ = []
    for icell, cell in enumerate(cell_lines):
        cell_save_fname = save_fname+"_"+cell+".csv"
        profiles = profiles_[[profiles_.columns[ix] for ix, x in enumerate(profiles_.loc["sigid"]) if (x.split("_")[1]==cell)]]
        initial_cols = profiles.loc["annotation"].values=="1.0"
        final_profiles = profiles[profiles.columns[~initial_cols]]
        conditions = list(set(final_profiles.loc["signame"]))
        initial_profiles = profiles[profiles.columns[initial_cols]].iloc[:-3,:].apply(pd.to_numeric)
        initial_profiles.columns = ["Ctrl_rep%d" % (i+1) for i in range(initial_profiles.shape[1])]
        final_profiles.columns = [x+"_%d" % (ix+1) for ix, x in enumerate(list(final_profiles.loc["signame"]))]
        final_profiles = final_profiles.iloc[:-3,:].apply(pd.to_numeric)
        if (not os.path.exists(cell_save_fname)):
            print("Cell %s (%d/%d)" % (cell, icell+1, len(cell_lines)))
            if (backgroundfile):
                endpoint = "sigs"
                method = "filter"
                params = {
                        "where": {"cell_id": cell, "pert_type": "trt_sh"}, 
                        "fields": ["distil_cc_q75", selection, "pct_self_rank_q25", "distil_id", "brew_prefix"]
                }
                request_url = build_url(lincs_api_url, endpoint, method, params=params, key=user_key)
                response = requests.get(request_url)
                assert response.status_code == 200
                data = json.loads(response.text)
                data = [dt for dt in data if (("distil_id" in dt) and ("distil_cc_q75" in dt) and ("pct_self_rank_q25" in dt))]
                assert len(data)>0
                print("#profiles %d" % len(data))
                ## good signatures
                data = [dt for dt in data if ((len(dt["distil_id"])>1) and (dt["distil_cc_q75"]>=0.2) and (dt["pct_self_rank_q25"]<=0.05))]
                assert len(data)>0
                print("# (good) profiles %d" % len(data))
                mselection = np.min([dt[selection] for dt in data])
                selection_min = 4
                max_selection = np.argsort(np.array([dt[selection] if (dt[selection]>=selection_min) else mselection for dt in data]))
                max_selection = max_selection[-min(50,len(max_selection)):]
                vals_selection = [dt[selection] for dt in [data[i] for i in max_selection]]
                print("# (best) profiles (capped at 50 or min>=%d) %d (%s max=%.3f, min=%.3f)" % (selection_min, len(max_selection), selection, np.max(vals_selection), np.min(vals_selection)))
                bkgrnd = create_restricted_drug_signatures([did for dt in [data[i] for i in max_selection] for did in dt["distil_id"]], cell, [int(g) for g in list(profiles.index)[:-3]], which_lvl=[3], strict=False, correct_sigids=True)
                bkgrnd.index = [str(g) for g in bkgrnd.index]
            ## aggregate replicates by median values
            initial_profile = initial_profiles.median(axis=1)
            final_profile = final_profiles.T
            final_profile.index = ["_".join(idx.split("_")[:-1]) for idx in final_profile.index]
            data = final_profile.groupby(level=0).median().T
            data.columns = conditions
            data["initial"] = initial_profile.loc[data.index]
            if (backgroundfile):
                data = data.join(bkgrnd, how="inner")
            data.to_csv(cell_save_fname)
        data = pd.read_csv(cell_save_fname, index_col=0)
        signatures = binarize_experiments(data, thres=thres, method=bin_method.split("_CD")[0], strict=not ('CD' in bin_method))
        signatures = signatures[conditions+["initial"]]
        if ("_CD" in bin_method):
            for c in conditions:
                from LINCS_utils import binarize_via_CD
                import random
                np.random.seed(0)
                random.seed(0)
                df = pd.concat((final_profiles[[col for col in final_profiles.columns if (c in col)]], initial_profiles), axis=1)
                samples = [int("Ctrl_" in col)+1 for col in df.columns]
                sigs = binarize_via_CD(df, samples=samples, binarize=1,nperm=10000)
                signatures[c] = list(sigs["aggregated"])
        signatures.columns = [s+"_"+cell for s in signatures.columns]
        signatures_.append(signatures)
    signatures_ = signatures_[0].join(signatures_[1:], how="outer")
    if (str(non_constant_genes)!="None"):
        profiles__ = profiles_.loc[[g for g in non_constant_genes if (g in profiles_.index)]]
        signames = list(profiles_.loc["signame"])
        sigids = list(signatures_.columns)
        vals = signatures_.values
        genes = list(signatures_.index)
        subgenes = list(profiles__.index)
        argmin_gene = np.argmin(profiles__.values, axis=1)
        argmax_gene = np.argmax(profiles__.values, axis=1)
        for ig, m in enumerate(argmin_gene):
            idx = [i for i, x in enumerate(sigids) if (signames[m] in x)][0]
            ig = genes.index(subgenes[ig])
            vals[ig, idx] = 0
        for ig, m in enumerate(argmax_gene):
            idx = [i for i,x in enumerate(sigids) if (signames[m] in x)][0]
            ig = genes.index(subgenes[ig])
            vals[ig, idx] = 1
        signatures_ = pd.DataFrame(vals, index=signatures_.index, columns=signatures_.columns)
    return signatures_

#' @param cell_lines Python list of string characters
#' @param pert_types Python list of string characters
#' @param thres_iscale Python number or None
#' @param gene_list Python list of string characters
#' @param thres_iscale Python nonnegative number
#' @param fname filename for saving signatures
#' @param thres_iscale Python number or None
#' @param nsigs Number of replicates for each condition
#' @param verbose Boolean to print things
#' @returns fname filename where signatures were saved as csv
def get_experimental_constraints(cell_lines, pert_types, gene_list, pert_inames, selection, thres_iscale=None, nsigs=2, verbose=1):
    assert str(thres_iscale) == "None" or thres_iscale >= 0
    assert len(cell_lines) > 0
    assert len(pert_types) > 0
    assert len(gene_list) > 0 and len(gene_list)==len(pert_inames)
    signatures = []
    perturbed_genes = []
    for line in cell_lines:
        ids = {}
        for pert_type in [p for p in pert_types if (p != "ctl_untrt")]:
            endpoint = "sigs"
            method = "filter"
            ## Select perturbation type, cell line
            where = {"pert_type": pert_type, "cell_id": line}
            ## Return perturbed gene name, perturbagen dose, perturbation type, cell line, exposure time, 
            ## distil id (unique to profile), brew prefix (unique to set of replicates)
            params = {
                    "where": where,
                    "fields": ["pert_iname", "pert_idose", "pert_type", "cell_id", "pert_itime", "distil_id", "brew_prefix"]
            }
            ## Create URL
            request_url = build_url(lincs_api_url, endpoint, method, params=params, key=user_key)
            ## Get results from request
            data_pert_ = post_request(request_url, quiet=(not verbose))
            #data_pert_ = []
            #for gene in pert_inames:
            #    p = {'where': {"pert_type": pert_type, "cell_id": line, "pert_iname": gene}, "fields": params["fields"]}
            #    res_pert_ = post_request(build_url(lincs_api_url, endpoint, method, p, key=user_key), quiet=True)
            #    if (len(res_pert_)>0):
            #        print("Result (%s,%s,%s): %d" % (gene,line,pert_type,len(res_pert_)))
            #    data_pert_ += res_pert_
            ## Trim out genes not in the considered gene list without enough replicates
            data_pert_ = [dt for dt in data_pert_ if ((dt["pert_iname"] in pert_inames) and (len(dt["distil_id"]) >= nsigs))]
            if (len(data_pert_)>0):
                entrez_id = str(gene_list[pert_inames.index(data_pert_[0]["pert_iname"])])
            for data in data_pert_:
                if (verbose):
                    print("* #Experiments so far: "+str(len(signatures)))
                signame = "_".join([str(data["pert_iname"]), "OE" if ("_oe" in pert_type) else "KD", data["pert_idose"] if (data["pert_idose"]!="-666") else "NA"])
                treatment = str(data["pert_iname"])
                ## avoid duplicates
                if (treatment in perturbed_genes):
                    print("*** duplicated treatment %s, cell %s, pert_type %s" % (treatment, str(data["cell_id"]), str(data["pert_type"])))
                    continue
                print("** Treatment %s (entrez_id %d)" % (treatment, int(entrez_id)))
                perturbed_genes.append(treatment)
                ## returns control & treated profiles from LINCS L1000
                sigs = get_treated_control_dataset(treatment, pert_type, line, {}, gene_list, entrez_id=entrez_id, add_untreated=False, 
                        which_lvl=[3], nsigs=nsigs, same_plate=True, selection=selection, 
                        quiet=(not verbose), iscale_thres=thres_iscale, trim_w_interference_scale=(str(thres_iscale)!="None"))
                if (str(sigs)=="None" or len(sigs)==0):
                    continue
                sigs.loc["signame"] = [signame]*sigs.shape[1]
                sigs.loc["sigid"] = list(sigs.columns)
                nexp = len(signatures)+1
                sigs.columns = ["Exp"+str(nexp)+":"+str(i)+"-rep"+str(ai+1) for ai, i in enumerate(list(sigs.loc["annotation"]))]
                signatures.append(sigs)
    signatures = signatures[0].join(signatures[1:], how="outer")
    signatures.to_csv(file_folder+"exp_profiles_cells="+"-".join(cell_lines)[:4]+"_perttypes="+"-".join(pert_types)+".csv")
    if (str(thres_iscale)!="None"):
        signatures = signatures[[c for ic, c in enumerate(signatures.columns) if (float(list(signatures.loc["interference_scale"])[ic])>thres_iscale)]]
        signatures = signatures.loc[[i for i in signatures.index if (i != "interference_scale")]]
    return signatures

#' @param signatures DataFrame from binarize_experiments
#' @param width, height dimension of image
#' @param max_show maximum number of genes shown (maximum variance across signatures)
#' @param fname name of image
#' @return None
def plot_signatures(signatures, width=10, height=10, max_show=50, fname="signatures.png"):
    from copy import deepcopy
    from matplotlib import colors,rc
    rc("ytick", labelsize=5)
    rc("xtick", labelsize=10)
    sigs = deepcopy(signatures)
    sigs[sigs == 0] = -1
    sigs[pd.isnull(sigs)] = 0 
    max_genes = np.argsort(np.var(sigs.values, axis=1))[-max_show:]
    max_genes = sigs.index[max_genes]
    selected_genes = list(set([y for y in [s.split("_")[0] for s in signatures.columns] if (y != "initial")]))
    for g in selected_genes:
        if (g not in sigs.index):
            sigs.loc[g] = 0
    max_genes = selected_genes+[g for g in max_genes if (g not in selected_genes)]
    sigs_ = sigs.loc[max_genes]
    fig, ax = plt.subplots(figsize=(width,height), nrows=1, ncols=1)
    cmap = colors.ListedColormap(['red', 'black', 'green'])
    bounds=[-1, -0.5, 0.5, 1]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    pos = ax.imshow(sigs_, cmap=cmap, origin="lower", interpolation='nearest', norm=norm) #remove score line
    plt.colorbar(pos, ax=ax, ticks=[-1, 0, 1], shrink=0.12, aspect=4)
    ax.set_xticks(range(sigs_.shape[1]))
    ax.set_xticklabels(sigs.columns, rotation=90)
    ax.set_yticks(range(sigs_.shape[0]))
    ax.set_yticklabels(sigs_.index)
    plt.savefig(fname, bbox_inches="tight")
    plt.close()
    return None

#' @param profiles DataFrame from get_experimental_constraints
#' @param fname name of image
#' @return None
def plot_distributions(profiles, fname="gene_expression_distribution.png", thres=None):
    bp = profiles.iloc[:-3,:].T.apply(pd.to_numeric).boxplot(rot=90, figsize=(25,15))
    if (str(thres)!="None"):
        K = profiles.shape[0]-3
        for t in [thres, 1-thres]:
            plt.plot(list(range(K)), [t]*K, "r--")
    plt.savefig(fname, bbox_inches="tight")
    plt.close()
    return None

def plot_discrete_distributions(signatures, fname="signature_expression_distribution.png"):
    N = int(np.sqrt(signatures.shape[1]))+1
    fig, axes = plt.subplots(nrows=N, ncols=N, figsize=(25,15))
    for iax in range(signatures.shape[1]):
        i,j = iax//N, iax%N
        sig = signatures.iloc[:, iax].fillna(0.5)
        axes[i,j].set_xticks([0,0.5,1])
        axes[i,j].set_xticklabels(('0','NaN','1'))
        axes[i,j].hist(sig.values.T)
        axes[i,j].set_ylabel("#genes")
        axes[i,j].set_title(signatures.columns[iax])
    plt.savefig(fname, bbox_inches="tight")
    plt.close()
    return None
