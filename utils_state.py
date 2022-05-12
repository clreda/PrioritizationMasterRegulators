# coding: utf-8

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

from utils_data import quantile_normalize
from params import njobs

#' @param data Pandas DataFrame
#' @returns binarized dataframe
def binarize_experiments(data, thres=0.5, method="binary", return_bin=True, strict=False):
    assert method in ["binary","probin"]
    if (method == "probin"):
        ## /!\ needs to install profilebinR (see install_env_synthesis.sh)
        from profile_binr import ProfileBin
        probin = ProfileBin(data.T)
        probin.fit(njobs)
        data = probin.binarize().T
        return data
    else:
        assert thres <= 0.5 and thres >= 0
        signatures = quantile_normalize(data)
        min_df = pd.concat([signatures.min(axis=1)] * signatures.shape[1], axis=1, ignore_index=True)
        min_df.columns = signatures.columns
        signatures = (signatures-min_df)
        max_df = pd.concat([signatures.max(axis=1)] * signatures.shape[1], axis=1, ignore_index=True)
        max_df.columns = signatures.columns
        signatures = signatures/max_df
        if (not return_bin):
            return signatures
        if (strict):
            signatures[signatures < thres] = 0
            signatures[signatures > 1-thres] = 1
            signatures[(signatures>=thres)&(signatures<=1-thres)] = np.nan
        else:
            signatures[signatures<=thres] = 0
            signatures[signatures>=1-thres] = 1
            signatures[(signatures>thres)&(signatures<1-thres)] = np.nan
        return signatures

#' @param b1 Pandas DataFrame: binary state (0, 1, NaN)
#' @param b2 Pandas DataFrame: binary state (0, 1, NaN)
#' @param genes string list: 
#' @param metric string: metric implemented in scikit-learn
#' @returns similarity between b1 and b2, number of genes on which it was computed
def compare_states(x, y, genes=None, metric="cityblock", nans="replace"):
    assert metric in ["cosine", "cityblock", "jaccard", "rogerstanimoto", "euclidean", "correlation"]
    assert nans in ["replace", "ignore"]
    normalized_metrics = ["jaccard", "rogerstanimoto", "cosine", "correlation"]
    xx = pd.DataFrame(x.values, index=x.index, columns=["X%d" %d for d in range(x.shape[1])])
    yy = pd.DataFrame(y.values, index=y.index, columns=["Y%d" %d for d in range(y.shape[1])])
    z = xx.join(yy, how="outer")
    if (nans == "ignore"):
        z = z.dropna()
        if (z.shape[0]==0):
            return 0., 0
    if (str(genes)!="None"):
        gene_list = [g for g in genes if (g in z.index)]
        z = z.loc[gene_list]
        if (len(gene_list)==0):
            raise ValueError
    N = z.shape[0]
    x_, y_ = [z[u.columns] for u in [xx,yy]]
    renorm = 1 if (metric in normalized_metrics) else (1/np.sqrt(N) if (metric == "euclidean") else 1/N)
    if (nans == "replace"):
        ## Compute separately distances between 1's and 0's
        X_pos, Y_pos, X_neg, Y_neg = [np.mod(u.fillna(v).T.values+1, 2).astype(bool if (metric in ["jaccard", "rogerstanimoto"]) else float) for v in [0,1] for u in [x_,y_]]
        dists_pos, dists_neg = [pairwise_distances(X, Y, metric=metric)*renorm for X,Y in [[X_pos,Y_pos],[X_neg,Y_neg]]]
        dists_pos[np.isnan(dists_pos)] = 0
        dists_neg[np.isnan(dists_neg)] = 0
        dists = np.power(np.multiply(dists_pos,dists_neg), 0.5) # geometric mean
    else:
        X, Y = [u.T.values.astype(bool if (metric in ["jaccard", "rogerstanimoto"]) else float) for u in [x_, y_]]
        dists = pairwise_distances(X, Y, metric=metric)*renorm
        dists[np.isnan(dists)] = 0
    ## numerical approximations
    dists[np.isclose(dists, 0)] = 0
    dists[np.isclose(dists, 1)] = 1
    sims = 1-dists
    return sims, N
