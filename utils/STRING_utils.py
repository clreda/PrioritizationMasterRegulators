#coding: utf-8

import sys
import pandas as pd

sys.path.insert(1, './credentials/')

from credentials_STRING import *

from time import sleep
import requests

string_api_url = "https://string-db.org/api"

#' @param genes Python list of character strings
#' @param score Python integer
#' @param taxon_id Python integer
#' @param fname Python character string
#' @return network Pandas DataFrame of undirected edges + scores
def get_network_from_STRING(genes, score=0, taxon_id=10090, fname="", quiet=True):
	assert score >= 0 and score <= 1000
	if (not quiet):
		print("* Getting the STRING name mapping for genes")
	output_format = "tsv"
	method = "get_string_ids"
	params = {
	    "identifiers" : "\r".join(genes), # your protein list
	    "species" : taxon_id, # species NCBI identifier
	    "limit" : 1, # only one (best) identifier per input protein
	    "echo_query" : 1, # see your input identifiers in the output
	    "caller_identity" : app_name # your app name
	}
	request_url = "/".join([string_api_url, output_format, method])
	results = requests.post(request_url, data=params).text
	sleep(1)
	with open(fname, "w+") as f:
		f.write(results)
	results = pd.read_csv(fname, sep="\t")
	id_di = {}
	my_genes = []
	for line in range(len(results.index)):
		input_identifier, string_identifier = results["queryItem"][line], results["stringId"][line]
		my_genes.append(string_identifier)
		id_di.setdefault(string_identifier, input_identifier)
	if (not quiet):
		print("* Getting the STRING network interactions")
	output_format = "tsv"
	method = "network"
	request_url = "/".join([string_api_url, output_format, method])
	params = {
		"identifiers" : "%0d".join(my_genes), # your protein
		"species" : taxon_id, # species NCBI identifier 
		"required_score" : score, # in 0 - 1000, 0 : get all edges
		"caller_identity" : app_name # your app name
	}
	if (not quiet):
		print("* Getting the STRING network interactions")
	response = requests.post(request_url, data=params).text
	with open(fname, "w+") as f:
		f.write(response)
	for line in response.strip().split("\n")[1:]:
		l = line.strip().split("\t")
		p1, p2 = l[2], l[3]
		## filter the interaction according to experimental score
		experimental_score = float(l[10])
		if experimental_score > 0.4:
			## print 
			if (not quiet):
				print("\t".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))
	sleep(1)
	network = pd.read_csv(fname, sep="\t")
	network["preferredName_A"] = [id_di.get(x, x) for x in list(network["preferredName_A"])]
	network["preferredName_B"] = [id_di.get(x, x) for x in list(network["preferredName_B"])]
	network.to_csv(fname, sep="\t")
	return network
