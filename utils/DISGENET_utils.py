# coding: utf-8

#source: https://www.disgenet.org/static/disgenet_rest/example_scripts/disgenet_api_request.py

import requests
import sys
import json
sys.path.insert(1, './credentials/')

import pandas as pd

from credentials_DISGENET import *

api_host = "https://www.disgenet.org/api/"

#' @param disease_list string list: list of Concept IDs
#' @param limit integer
#' @param source string: accepted by DisGeNet API
#' @param min_score float: minimum global score
#' @param min_ei float: minimimum Evidence Index
#' @param min_dsi float: minimum Disease Specificity Index
#' @param min_dpi float: minimum Disease Pleiotropy Index
#' @param sz integer: size of chunks
#' @returns res_df Pandas DataFrame: rows/diseases x columns/[Protein name, Gene symbol]
def get_genes_proteins_from_DISGENET(disease_list, limit=3000, source="CURATED", min_score=0, min_ei=0, min_dsi=0.25, min_dpi=0, sz=100):
    assert user_key
    assert sz <= 100 and sz > 1
    res_df = []
    for i in range(0, len(disease_list), sz):
        print("%d/%d" % (i+1, len(disease_list)))
        request_url = api_host+"gda/disease/"+",".join(disease_list[i:(i+sz)])
        params = {
                "source": source, 
                "format": "json", 
                "limit": limit, 
                "min_score": min_score, #GDA Score
                "max_score": 1, 
                "min_ei": min_ei, #Evidence Level
                "max_ei": 1,
                "min_dsi": min_dsi, #Disease Specificity Index: min 0.25
                "max_dsi": 1,
                "min_dpi": min_dpi, #Disease Pleiotropy Index
                "max_dpi": 1,
        }
        headers = {"Authorization": "Bearer %s" % user_key}
        request_url += "?"+"&".join([p+"="+str(params[p]) for p in params])
        r = requests.get(request_url, headers=headers)
        if (r.status_code not in [200, 404]):
            print(r.text)
            print(request_url)
            raise ValueError("Request failed.")
        if (r.status_code == 404):
            return None
        res = pd.DataFrame(json.loads(r.text))
        res.index = res["diseaseid"]
        vals_genes = res[["gene_symbol"]].groupby(level=0).apply(lambda x : "; ".join(list(sorted(set(list(x.values.flatten()))))))
        vals_proteins = res[["uniprotid"]].groupby(level=0).apply(lambda x : "; ".join(list(map(str,set(x.values.flatten())))))
        res = pd.concat([vals_proteins, vals_genes], axis=1)
        res_df.append(res)
    res_df = pd.concat(res_df)
    res_df.columns = ["Protein", "Gene Name"]
    return res_df

#' @param gene_list string list
#' @param disease string: Concept ID
#' @param limit integer: # evidence
#' @param source string: recognized by DisGeNet API
#' @param min_score float: minimum evidence score
#' @param sz integer: size of chunks
#' @returns res_df Pandas DataFrame: rows/genes x ["symbol", "description", "references", "year", "evidence score"]
def get_genes_evidences_from_DISGENET(gene_list, disease, limit=3000, source="CURATED", min_score=0, sz=100):
    assert user_key
    assert sz <= 100 and sz > 1
    res_df = []
    for i in range(0, len(gene_list), sz):
        print("%d/%d" % (i+1, len(gene_list)))
        request_url = api_host+"gda/evidences/gene/"+",".join(gene_list[i:(i+sz)])
        params = {
                "disease": disease,
                "source": source, 
                "format": "json", 
                "limit": limit, 
                "min_score": min_score, #GDA Score
                "max_score": 1, 
        }
        headers = {"Authorization": "Bearer %s" % user_key}
        request_url += "?"+"&".join([p+"="+str(params[p]) for p in params])
        r = requests.get(request_url, headers=headers)
        if (r.status_code not in [200, 404]):
            print(r.text)
            print(request_url)
            raise ValueError("Request failed.")
        if (r.status_code == 404):
            return None
        res = pd.DataFrame(json.loads(r.text))
        ls = list(res["results"])
        res = pd.DataFrame({d: lss for d,lss in enumerate(list(res["results"]))}).T
        res = res.loc[(res["score"]>=min_score)&(res["disease_id"]==disease)]
        if (len(res)>0):
            res_df.append(res)
    if (len(res_df)>0):
        res_df = pd.concat(res_df)
        res_df = res_df[["gene_symbol","sentence","associationtype","pmid","year","score"]]
    else:
        res_df = None
    return res_df

## Test
if __name__ == "__main__":
    epilepsy = "C0014544"
    gene_df = get_genes_proteins_from_DISGENET([epilepsy], min_score=0)
    print(len(list(gene_df["Gene Name"])[0].split("; ")))
    gene_df = get_genes_proteins_from_DISGENET([epilepsy], min_score=0.5)
    if (str(gene_df)!='None'):
        print(len(list(gene_df["Gene Name"])[0].split("; ")))
    else:
        print(0)
    gene_ls = list(gene_df["Gene Name"])[0].split("; ")
    evidences = get_genes_evidences_from_DISGENET(gene_ls, epilepsy, min_score=0)
    print(evidences.shape[0])
    evidences = get_genes_evidences_from_DISGENET(gene_ls, epilepsy, min_score=0.65)
    print(evidences.shape[0])
