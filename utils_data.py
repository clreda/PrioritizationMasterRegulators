#coding: utf-8

import numpy as np
import pandas as pd
import json
import subprocess as sb
from time import sleep
import os
from glob import glob
from qnorm import quantile_normalize as qnm

import io

from params import *

## From qnorm module
quantile_normalize = lambda df : qnm(df, axis=1, ncpus=njobs)

#https://biodbnet-abcc.ncifcrf.gov/webServices/RestWebService.php
#' @param probes string list
#' @param from_ string: recognized by BioDBnet
#' @param to_ string: recognized by BioDBnet
#' @param taxonId integer: recognized by NCBI taxon numbers
#' @param chunksize integer
#' @return result Pandas DataFrame: first column query, second column matches
def request_biodbnet(probes, from_="Illumina ID", to_="Gene Symbol", taxonId=10090, chunksize=500):
    chunk_probes=[probes[i:i+chunksize] for i in range(0,len(probes),chunksize)]
    result = []
    url="https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json"
    for i, chunk in enumerate(chunk_probes):
        print("%d/%d" % (i, len(chunk_probes)))
        args_di = {
            "method":"db2db",
            "format": "row",
            "inputValues":",".join([x.upper() for x in chunk]),
            "input":from_,
            "outputs":to_,
            "taxonId":str(taxonId),
        }
        #https://biodbnet-abcc.ncifcrf.gov/webServices/RestSampleCode.php
        query=url+"?"+"&".join([k+"="+args_di[k] for k in args_di])
        result += json.loads("; ".join(sb.check_output("wget -qO- \""+query+"\"", shell=True).decode("utf-8").split("//")))
        sleep(1)
    result = pd.DataFrame(result)
    result.index = result["InputValue"]
    result = result.drop(columns=["InputValue"])
    return result
