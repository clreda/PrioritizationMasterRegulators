#coding: utf-8

import sys
import subprocess as sb

sys.path.insert(1, './credentials/')
sys.path.insert(1, "../")

from params import data_folder, taxon_id
from credentials_LINCS import *

import os
import pandas as pd
import numpy as np
import cmapPy.pandasGEXpress.parse_gctx as parse_gctx

from functools import reduce

path_to_lincs = data_folder+"lincs/"
lincs_specific_ctl_genes = ["ACTB", "CHMP2A", "EEF1A1", "EMC7", "GAPDH", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "TUBA1A", "VCP"]
ctl_gene_fname = "utils/LINCS_data/control_genes.csv"
sb.call("mkdir -p "+"/".join(ctl_gene_fname.split("/")[:-1]), shell=True)
if (not os.path.exists(ctl_gene_fname)):
    from utils_data import request_biodbnet
    lincs_specific_ctl_genes_ = [np.min(list(map(int, s.split("; ")))) for s in request_biodbnet(lincs_specific_ctl_genes, from_="Gene Symbol and Synonyms", to_="Gene ID", taxonId=taxon_id, chunksize=100)["Gene ID"]]
    pd.DataFrame(lincs_specific_ctl_genes_, columns=["Entrez_ID"], index=lincs_specific_ctl_genes).to_csv(ctl_gene_fname)
lincs_specific_ctl_genes = list(map(int,pd.read_csv(ctl_gene_fname, index_col=0)["Entrez_ID"]))

####################
## API CALLS      ##
####################

lincs_api_url = "https://api.clue.io/api"

#' @param api_url Python character string
#' @param endpoint Python character string
#' @param method Python character string
#' @param Python dictionary
#' @param key Python character string
#' @return url Python character string: URL of request
def build_url(api_url, endpoint, method=None, params={}, key=None):
	assert key
	request_url = "/".join([lincs_api_url, endpoint])
	convert_dict = lambda lst : '"'.join(str(lst).split("\'"))
	if (len(params) > 0):
		params_concat = "&".join([k+"="+convert_dict(params[k]) for k in list(params.keys())])
	else:
		params_concat = ""
	if (method == "count"):
		if (len(params) > 0):
			params_concat = "&".join([k+"="+convert_dict(params[k]) for k in list(params.keys())])
		else:
			params_concat = ""
		request_url += "/"+method+"?"+params_concat
	elif (method == "filter"):
		if (len(params) > 0):
			params_concat = "{"+(",".join(['"'+k+'":'+convert_dict(params[k]) for k in list(params.keys())]))+"}"
		else:
			params_concat = ""
		request_url += "?"+method+"="+params_concat
	elif (method == "distinct"):
		assert params["field"]
		if ("where" in list(params.keys())):
			request_url += "/"+method+"?where="+convert_dict(params["where"])+"&field="+params["field"]
		else:
			request_url += "/"+method+"?field="+params["field"]
	else:
		raise ValueError
	return request_url+"&user_key="+user_key

#' @param url Python character string
#' @param quiet OPTIONAL Python bool
#' @return data Python character string (JSON): result of request
def post_request(url, quiet=False, stime=1):
	from time import sleep
	import subprocess as sb
	import json
	if (not quiet):
		print("> POST "+url.split("&")[0])
	response = sb.check_output("wget -O - -q \'" + url + "\'", shell=True)
	data = json.loads(response)
	if ("/distinct?field=" in url):
		data = list(set([y for x in data for y in x]))
	sleep(stime)
	return data

##<><><><><><><><><><><><><><><><><><><><><><><><><><><><>##
##   Binarization of Level 3 signatures                   ##
##<><><><><><><><><><><><><><><><><><><><><><><><><><><><>##

## Implementation of CD in http://www.maayanlab.net/CD/
#' @param df Pandas DataFrame of non-binary signatures
#' @param samples Python list of integers (1: control or 2: treated)
#' @param nperm number of iterations for p-value computation
#' @return binary_sig binary (UP/DOWN 1/0 regulated NA otherwise) signature for treated
def binarize_via_CD(df, samples=[], binarize=1, nperm=10000, thres_pval=0.01):
    sb.check_output("pip install git+https://github.com/Maayanlab/geode.git", shell=True)
    from geode import chdir
    assert len(df.columns) >= 4
    assert len(df.columns) == len(samples)
    assert len(np.argwhere(np.array(samples) == 1).tolist()) > 1
    assert len(np.argwhere(np.array(samples) == 2).tolist()) > 1
    df = df.dropna()
    ## returns a list of tuples (signed magnitude CD, gene name, p-value) sorted by the absolute value in descending 
    ## order characteristic directions of genes
    ## automatically computes significant genes (calculate_sig=1 computes p-values,  sig_only=1 returns only significant genes)
    ## for memory reasons
    if (len(df.index) > 25000):
        from random import sample
        print("WARNING: The number of genes has been trimmed. Consider filtering your genes beforehand!")
        df = df.iloc[sample(range(len(df.index)), 25000),:]
    chdir_res = chdir(df.values, samples, list(df.index), calculate_sig=binarize, nnull=nperm, sig_only=binarize)
    significant_genes = list(map(lambda x : x[:2], chdir_res))
    genes = list(df.index)
    significant_gene_list = list(map(lambda y : y[1], significant_genes))
    if (binarize>0):
        signature = [int(significant_genes[significant_gene_list.index(x)][0] > 0) if (x in significant_gene_list) else np.nan for x in genes]
    else:
        signature = [significant_genes[significant_gene_list.index(x)][0] if (x in significant_gene_list) else np.nan for x in genes]
    sig = pd.DataFrame(signature, index=genes, columns=["aggregated"])
    return sig

####################
## DOWNLOAD FILES ##
####################

#' @param path Python character string
#' @param file_name Python character string
#' @param base_url Python character string
#' @param file_sha Python character string
#' @param check_SHA Python bool
#' @return None
def download_file(path, file_name, base_url, file_sha, check_SHA=True):
    import subprocess as sb
    if (not os.path.exists(path+file_name)):
        if (not os.path.exists(path+file_name+".gz")):
            cmd = "wget -O "+path+file_name+".gz "+base_url+file_name+".gz"
            print(cmd)
            sb.call(cmd, shell=True)
            if (check_SHA):
                ## Checks file integrity
                sha_table = pd.read_csv(path+file_sha, sep="  ", names=["sha", "file"], engine='python')
                true_sha_id = list(sha_table[sha_table["file"] == file_name+".gz"]["sha"])[0]
                print(path+file_name)
                sha_id = sb.check_output("sha512sum "+path+file_name+".gz", shell=True).decode("utf-8").split("  ")[0]
                assert sha_id == true_sha_id
        sb.call("gzip -df "+path+file_name+".gz", shell=True)
        print(file_name+" successfully downloaded")
    else:
        ## Checks file integrity
        if (check_SHA and os.path.exists(path+file_name+".gz")):
            sha_table = pd.read_csv(path+file_sha, sep="  ", names=["sha", "file"], engine='python')
            true_sha_id = list(sha_table[sha_table["file"] == file_name+".gz"]["sha"])[0]
            sha_id = sb.check_output("sha512sum "+path+file_name+".gz", shell=True).decode("utf-8").split("  ")[0]
            if (not (sha_id == true_sha_id)):
                cmd = "wget -c -O "+path+file_name+".gz "+base_url+file_name+".gz"
                print(cmd)
                sb.call(cmd, shell=True)
                if (check_SHA):
                    ## Checks file integrity
                    sha_table = pd.read_csv(path+file_sha, sep="  ", names=["sha", "file"], engine='python')
                    true_sha_id = list(sha_table[sha_table["file"] == file_name+".gz"]["sha"])[0]
                    sha_id = sb.check_output("sha512sum "+path+file_name+".gz", shell=True).decode("utf-8").split("  ")[0]
                    assert sha_id == true_sha_id
                sb.call("gzip -df "+path+file_name+".gz", shell=True)
                print(file_name+" successfully downloaded")

#' @param path Python character string
#' @param which_lvl Python integer
#' @return gene_files, sig_files, lvl3_files, lvl5_files Python lists of character strings
def download_lincs_files(path, which_lvl):
	assert all([x in [3,5] for x in which_lvl])
	## Lastest versions of LINCS
	lincs_gse = {
		"phase1": {
			"acc": "GSE92742",
			"file_lvl5": "Level5_COMPZ.MODZ_n473647x12328.gctx",
			"file_lvl3": "Level3_INF_mlr12k_n1319138x12328.gctx",
			"gene_file": "gene_info.txt",
			"sig_file": "sig_info.txt"
		},
		"phase2": {
			"acc": "GSE70138",
			"file_lvl5": "Level5_COMPZ_n118050x12328_2017-03-06.gctx",
			"file_lvl3": "Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx",
			"gene_file": "gene_info_2017-03-06.txt",
			"sig_file": "sig_info_2017-03-06.txt"
		}
	}
	keys = ["gene_file", "sig_file", "file_lvl3", "file_lvl5"]
	file_di = {}
	for k in keys:
		file_di.setdefault(k, [])
	for key in list(lincs_gse.keys()):
		acc = lincs_gse[key]["acc"]
		base_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/"+acc[:5]+"nnn/"+acc+"/suppl/"
		for k in keys:
			file_di[k] = file_di[k]+[acc+"_Broad_LINCS_"+lincs_gse[key][k]]
		file_sha = acc+"_SHA512SUMS.txt"
		download_file(path, file_sha, base_url, None, check_SHA=False)
		## 5GB! ~20 minutes with a good Internet connection (level 5)
		## 48.8GB and 12.6GB! (level 3)
		## if the download stops, you can resume it using option "-c" of wget
		for lvl in which_lvl:
			download_file(path, file_di["file_lvl"+str(lvl)][-1], base_url, file_sha)
		for k in keys[:2]:
			download_file(path, file_di[k][-1], base_url, file_sha)	
	return [file_di[k] for k in keys]

####################
## SIGNATURES     ##
####################

## Create drug signature dataset from LINCS L1000 dataset
#' @param sig_ids Python list of character strings
#' @param cell_line Python character string
#' @param gene_list Python list of character string
#' @param which_lvl Python integer
#' @param path Python character string
#' @param strict return all signatures in sig_ids or return None
#' @return sigs Pandas DataFrame of signatures
def create_restricted_drug_signatures(sig_ids, cell_line, gene_list, which_lvl=[3], path=path_to_lincs, strict=True, correct_sigids=False):
    assert all([x in [3,5] for x in which_lvl])
    gene_files, sig_files, lvl3_files, lvl5_files = download_lincs_files(path, which_lvl=which_lvl)
    sigs = None
    cid = ("sig_id" if (5 in which_lvl) else "distil_id")
    for idf, sf in enumerate(sig_files):
        df = pd.read_csv(path+sf, sep="\t", engine='python')
        if (3 in which_lvl):
            sig_selection = []
            distil_proc = lambda sig : "_".join(list(map(lambda x : x.split(".")[0], sig.split("_"))))
            for sig in sig_ids:
                cids = list(filter(lambda x : "_".join(sig.split("_")[:2]).split(".")[0] in x and sig.split("_")[-1] in x, list(df[cid])))
                if (len(cids) > 0):
                    cids = [y for x in cids for y in x.split("|")]
                    cid__ = list(filter(lambda x : distil_proc("_".join(sig.split("_")[:3])) in x and sig.split("_")[-1] in x, cids))
                    cid_ = list(filter(lambda x : sig.split("_")[3] in x, cid__))
                    if (len(cid_) == 0 or all([c in sig_selection for c in cid_])):
                        cid_ = list(filter(lambda x : x not in sig_selection and x not in map(distil_proc, sig_ids), cid__))
                        if (len(cid_) == 0):
                            cid_ = list(filter(lambda x : sig.split("_")[1].split(".")[0] in x and sig.split("_")[-1].split(":")[0] in x and x not in sig_selection, cids))
                if (len(cid_)==0):
                    return None
                sig_selection.append(cid_[0])
            sig_selection = list(set(sig_selection))
        else:
            sig_selection = list(filter(lambda sig : any([sig == cids for cids in list(df[cid])]), sig_ids))
        sig_selection = sig_ids if (not correct_sigids) else sig_selection
        if (len(sig_selection) > 0):
            df = pd.read_csv(path+gene_files[idf], sep="\t", engine='python')
            if (gene_list is None):
                gene_selection = list(df.index)
            else:
                gene_selection = [s for s in df.index if (int(df.loc[s]["pr_gene_id"]) in gene_list)]
            try:
                #https://github.com/cmap/cmapPy/blob/master/cmapPy/pandasGEXpress/parse_gctx.py
                dataset = parse_gctx.parse(path+(lvl5_files if (which_lvl == 5) else lvl3_files)[idf], cid=sig_selection, ridx=gene_selection).data_df
            except Exception as e:
                # means "not in metadata for the considered GCTX file"
                print(e)
                continue
            dataset.index = list(map(int,list(df.loc[gene_selection]["pr_gene_id"])))
            assert dataset.shape[0]==len(gene_selection)
            assert dataset.shape[1]==len(sig_selection)
            if (str(sigs) == "None"):
                sigs = dataset
            else:
                sigs = sigs.join(dataset, how="outer")
    if (str(sigs) == "None"):
        print("\nWARNING! None of the following signature/profile ids was found:"),
        print("-> {"+"  ".join(sig_ids)+"}")
        if (strict):
            raise ValueError
        else:
            return None
    assert (len(sig_ids)==sigs.shape[1]) or (not strict)
    if (str(sigs)!="None" and ((len(sigs.columns) < len(sig_ids)) and strict)):
        print("\nWARNING! Some of the following signature/profile ids were not found:"),
        print("-> {"+"  ".join(list(filter(lambda x : distil_proc(x) not in sigs.columns, sig_ids)))+"}")
        return None
    return sigs

## Select "best" signature according to filters
#' @param params Python dictionary
#' @param filters Python dictionary
#' @param selection Python character string
#' @param nsigs Python integer
#' @param same_plate Python bool
#' @param quiet Python bool
#' @return data Python list of dictionaries
def select_best_sig(params, filters, selection="distil_ss", nsigs=2, same_plate=True, quiet=False, iunit=None):
	assert len(selection) == 0 or (selection in list(filters.keys()) or selection in params.get("fields", []))
	id_fields = ["brew_prefix"]#,"distil_ss"]
	if (any([not f in params.get("fields", []) for f in id_fields]) or (len(params.get("fields", [])) == 0 or (len(filters) > 0 and any([not (f in params.get("fields", [])) for f in list(filters.keys())])))):
		fields = params.get("fields", [])
		fields += list(filters.keys())
		params["fields"] = list(set(fields+id_fields))
	endpoint = "sigs"
	method = "filter"
	request_url = build_url(lincs_api_url, endpoint, method, params=params, key=user_key)
	data = post_request(request_url, quiet=quiet)
	if (nsigs > 0 and (len(data) == 0 or (same_plate and all([len(d["distil_id"]) < nsigs for d in data])) or (not same_plate and sum([len(d["distil_id"]) for d in data]) < nsigs))):
		if (not quiet):
			print("(1) No (enough) signatures (%d instead of min. %d) retrieved via LINCS L1000 API.\n%s\n" % (sum([len(d["distil_id"]) for d in data]), nsigs, request_url))
		raise ValueError#return []
	if (not quiet):
		print("- Number of data: "+str(len(data)))
	## Filter signatures
	if (len(filters) > 0):
		data = list(filter(lambda x : all([x.get(k, filters[k]+1) > filters[k] for k in list(filters.keys())]), data))
		if (not quiet):
			print("- Number of filtered data: "+str(len(data)))
	if ("pert_dose_unit" in params["fields"] and iunit):
		data = list(filter(lambda x : x["pert_dose_unit"] == iunit.decode("utf-8"), data))
	## Keep only valid id
	if ("distil_id" in params["fields"]):
		data_ = []
		for i in range(len(data)):
			distil_ids = list(filter(lambda x : len(x.split(":")) == 2, data[i]["distil_id"]))
			if (len(distil_ids) > 0):
				d = data[i]
				d["distil_id"] = distil_ids
				data_.append(d)
		data = data_
	if (same_plate and "distil_id" in params["fields"] and len(data) > 0):
		## Select the "best" signature in terms of "distil_ss" value with enough DISTINCT replicates
		data = list(filter(lambda x : len(list(set(x["distil_id"]))) >= nsigs, data))
		argmax_rank = int(np.argmax([float(x.get(selection, 0)) for x in data]))
		data_rank = data[argmax_rank]
		if (not quiet):
			print(selection+" = "+str(data_rank[selection]))
		data = [{"distil_id": distil_id} for distil_id in list(set(data_rank["distil_id"]))]
		for f in list(filter(lambda x : x != "distil_id", params["fields"])):
			for i in range(len(data)):
				data[i].setdefault(f, data_rank[f])
		if (not quiet):
			print("-- Number of same-plate data: "+str(len(data)))
	elif ("distil_id" in params["fields"]):
		for i in range(len(data)):
			data[i]["distil_id"] = data[i]["distil_id"][0]
	return data

## Using PMC5192966 shRNA/cDNA interference scale
#' @param sigs Pandas DataFrame of signatures
#' @param samples list of 1: control, 2: treated for each columns of sigs
#' @param entrez_id integer 
#' @param pert_type string character
#' @param quiet boolean
#' @return interference scale for the input experiment
def compute_interference_scale(sigs, samples, entrez_id, pert_type, quiet=True, eps=0.):
    assert len(sigs.columns) == len(samples)
    entrez_id = int(entrez_id)
    if (entrez_id not in sigs.index):
        return 0.
    ## ""Housekeeping genes"" selected in PMC5192966
    treated = sigs.columns[np.array(samples)==2]
    control = sigs.columns[np.array(samples)==1] 
    exp_ratio = np.mean(sigs[treated].loc[entrez_id])/float(np.mean(sigs[control].loc[entrez_id])+eps)
    ctl_genes = [int(g) for g in lincs_specific_ctl_genes if (g in sigs.index)]
    assert len(ctl_genes) > 0
    hk_inv_ratios = [np.mean(sigs[control].loc[g])/float(np.mean(sigs[treated].loc[g])+eps) for g in ctl_genes]
    ## The most stable (ratio ~ 1) between control and treated groups
    best_hk_gene_id = int(np.argmin(np.abs(np.array(hk_inv_ratios)-1)))
    iscale = exp_ratio*hk_inv_ratios[best_hk_gene_id]
    if (pd.isna(iscale)):
        return 0.
    if (pert_type in ["trt_sh", "trt_sh.cgs", "trt_xpr"]):
        iscale = 1-iscale
    else:
        iscale = iscale-1
    if (not quiet):
        print((entrez_id,exp_ratio))
        print((ctl_genes[best_hk_gene_id],1./hk_inv_ratios[best_hk_gene_id]))
        print(("interference scale", iscale))
        print("\n\n* Interference scale for the experiment is: "+str(iscale)+"\n\n")
    return iscale

## Get treated/control dataset + annotation (nsigs treated, nsigs untreated; from the same plate if same_plate=True)
#' @param treatment Python character string
#' @param pert_type Python character string
#' @param cell Python character string
#' @param filters Python dictionary
#' @param genes Python list of character strings
#' @param add_untreated Python bool
#' @param selection Python character string
#' @param dose OPTIONAL Python character string
#' @param iunit OPTIONAL Python character string
#' @param itime OPTIONAL Python character string
#' @param nsigs Python integer
#' @param same_plate Python bool
#' @param quiet Python bool
#' @param path Python character string
#' @param trim_w_interference_scale Python bool
#' @param return_metrics Python character string list
#' @return sigs Pandas DataFrame of signatures
def get_treated_control_dataset(treatment, pert_type, cell, filters, genes, entrez_id=None, add_untreated=False, selection="distil_ss", dose=None, iunit=None, itime=None, which_lvl=[[3,5][1]], nsigs=2, same_plate=True, quiet=False, trim_w_interference_scale=True, return_metrics=[], path=path_to_lincs, iscale_thres=0.):
    assert all([x in [3,5] for x in which_lvl])
    genes = list(set(genes))
    ## Filters
    where = {"cell_id": cell, "pert_type": pert_type}
    if (len(treatment) > 0):
        where.setdefault("pert_iname", treatment)
    if (str(dose) != "None"):
        where.setdefault("pert_dose", dose)
    if (str(itime) != "None"):
        where.setdefault("pert_itime", itime)
    # https://www.biostars.org/p/211896/
    cid = ("sig_id" if (which_lvl == 5) else "distil_id")
    fields = list(set(list(filters.keys())+[cid]+([selection] if (len(selection) > 0) else [])+return_metrics))
    params = {	
        "where": where,
        "fields": fields+["pert_dose_unit"]
    }
    ## Looking for the ids for profiles treated with treatment
    data_treated = select_best_sig(params, filters=filters, selection=selection, nsigs=nsigs, same_plate=same_plate, quiet=quiet, iunit=iunit)
    if (len(data_treated) == 0):
        if (not quiet):
            print("(1) No treated signature available.\n")
        return None
    where = {"cell_id": cell, "pert_type": ("ctl_vehicle" if (pert_type == "trt_cp") else "ctl_vector")}
    if (same_plate):
        where.setdefault("brew_prefix", str(data_treated[0]["brew_prefix"][0]))
    ## Looking for the ids for profiles control (in the same plate if specified)
    params = {	
        "where": where,
        "fields": fields
    }
    data_control = select_best_sig(params, filters=filters, selection=selection, nsigs=nsigs, same_plate=same_plate, quiet=quiet)
    if (len(data_control) == 0):
        if (not quiet):
            print("(2) No control signature available.\n")
        return None
    ## Add (if needed) untreated profiles
    if (add_untreated):
        where = {"cell_id": cell, "pert_type": "ctl_untrt"}
        if (same_plate):
            where.setdefault("brew_prefix", str(data_treated[0]["brew_prefix"][0]))
        params = {	
            "where": where,
            "fields": fields
        }
        ## no replicate for untreated experiments
        data_untreated = select_best_sig(params, filters=filters, selection=selection, nsigs=1, same_plate=False, quiet=quiet)
        if (len(data_untreated) == 0):
            if (not quiet):
                print("(3) No untreated signature available.\n")
            return None
    else:
        data_untreated = []
    sig_ids = list(map(lambda x : str(x[cid]), data_treated+data_control+data_untreated))
    if (not quiet):
        print("Got signature(s) "+reduce(lambda x,y : x+","+y, sig_ids))
    if (trim_w_interference_scale and pert_type != "trt_cp"):
        genes_ = list(set(lincs_specific_ctl_genes+genes))
    else:
        genes_ = genes
    sigs = create_restricted_drug_signatures(sig_ids, cell, genes_, which_lvl=which_lvl, path=path, strict=True)
    if (str(sigs) == "None"):
        if (not quiet):
            print("(4.a) No signature available from LINCS data file.\n")
        return None
    if ((pert_type != "trt_cp") and (not int(entrez_id) in sigs.index)):
        if (not quiet):
            print("(4.b) Treatment "+treatment+" not getable from data (due to missing gene in LINCS L1000).\n")
        return None
    ## Create noisy replicate (if needed) for untreated/drug control experiment
    if (len(data_untreated) > 0 and nsigs > 1):
        if (not quiet):
            print("* Replication step for experiment "+data_untreated[0][cid])
        repl_sig = sigs[[data_untreated[0][cid]]]
        vals = repl_sig.values.flatten().tolist()
        lambda_ = 1/float(len(vals))*float(np.var(vals))
        for n in range(1, nsigs):
            noisy_replicate = repl_sig.values+np.random.normal(np.zeros(np.shape(repl_sig.values)), lambda_)
            noisy_column = "rep_X_"+str(n)+"_"+repl_sig.columns[0]
            sigs.insert(nsigs+1, noisy_column, noisy_replicate, True)
            from copy import deepcopy
            data_replicate = deepcopy(data_untreated[0])
            data_replicate[cid] = noisy_column
            if (len(return_metrics) > 0):
                for metric in return_metrics:
                    data_replicate[metric] = np.nan
            data_untreated.append(data_replicate)
    ## Annotations
    signame_trt, signame_ctl, signame_untrt = [[str(s[cid]) for s in x] for x in [data_treated, data_control, data_untreated]]
    samples = [(2 if (c in signame_trt) else (1 if (c in signame_ctl) else (3 if (c in signame_untrt) else 0))) for c in sigs.columns]
    assert 0 not in samples
    sigs.loc["annotation"] = samples
    if (len(return_metrics) > 0):
        for metric in return_metrics:
            sigs = sigs.append(pd.Series([data.get(metric, np.nan) for data in data_treated+data_control+data_untreated], index=sigs.columns, name=metric), ignore_index=False)
    ## Using PMC5192966 shRNA/cDNA interference scale
    if (trim_w_interference_scale and pert_type != "trt_cp"):
        select = list(filter(lambda col : sigs.loc["annotation"][col] in [1, 2], sigs.columns))
        iscale = compute_interference_scale(sigs[select].drop(["annotation"]), sigs[select].loc["annotation"], entrez_id, pert_type, quiet=quiet)
        sigs.loc["interference_scale"] = [iscale]*sigs.shape[1]
    sigs = sigs.loc[list(filter(lambda x : x in list(sigs.index), genes))+["annotation"]+return_metrics+([] if (not trim_w_interference_scale) else ["interference_scale"])]
    return sigs

def get_untreated_dataset(cell, filters, genes, which_lvl=[[3,5][1]], quiet=False, return_metrics=[], path=path_to_lincs):
    assert all([x in [3,5] for x in which_lvl])
    genes = list(set(genes))
    where = {"pert_type": "ctl_untrt"}
    cid = ("sig_id" if (which_lvl == 5) else "distil_id")
    fields = list(set(list(filters.keys())+[cid]))
    params = {	
            "where": where,
            "fields": fields
    }
    ## no replicate for untreated experiments
    data_untreated = select_best_sig(params, filters=filters, selection="", nsigs=0, same_plate=False, quiet=quiet)
    if (len(data_untreated) == 0):
        if (not quiet):
            print("(3) No signature available.\n")
        return None
    sig_ids = list(map(lambda x : str(x[cid]), data_untreated))
    if (not quiet):
        print("Got signature(s) "+reduce(lambda x,y : x+","+y, sig_ids))
    sigs = create_restricted_drug_signatures(sig_ids, cell, genes, which_lvl=which_lvl, path=path)
    samples = [1]*len(sigs.columns)
    sigs_cols = [si for si, s in enumerate(sig_ids) if (s in sigs.columns)]
    ids_ = np.argwhere(np.array([str(data["cell_id"]) for data in [data_untreated[i] for i in sigs_cols]]) == cell)
    for i in ids_:
        samples[i] = 2
    ## Annotations
    sigs = sigs.append(pd.Series(samples, index=sigs.columns, name="annotation"), ignore_index=False)
    sigs = sigs.loc[list(filter(lambda x : x in list(sigs.index), genes))+["annotation"]]
    return sigs
