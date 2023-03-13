#coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sb
import os
import json
from time import sleep

from utils_data import request_biodbnet
from utils_grn import solution2influences
from params import *

##################
### PARAMETERS ###
##################

rfolder = root_folder+name+"/rankings/"
ddfolder = rfolder+"data/"
vfolder = rfolder+"values/"
pfolder = rfolder+"plots/"
## Threshold for building ORA gene lists
thres = {'pLI': 0.9999, "Spread": 0.01, "CC": 1, "Influence": 0}
sb.call("mkdir -p "+vfolder+" "+pfolder+" "+ddfolder, shell=True)
french=False

##################
## GENE VALUES  ##
##################

## Consider ALL nodes in the network
network_name=solution_folder+"solution.bnet"
with open(network_name, "r") as f:
    all_genes = [line.split(" <- ")[0] for line in f.read().split("\n")]

## Influence TF [Nicolle, Radvanyi, and Elati, 2015] in CoRegNet
with open(network_name, "r") as f:
    lines = f.read().split("\n")
    genes = [line.split(" <- ")[0] for line in lines]
    solution = pd.DataFrame([line.split(" <- ")[-1] for line in lines], index=genes, columns=["solution"])
influences = solution2influences(solution["solution"])
initial_states = pd.read_csv(path_to_initial_states, index_col=0)
def influence_t_per_state(t, state):
    regulated_t = list(influences.loc[t])
    state = state.fillna(0)
    if (len(set(regulated_t))==1 and regulated_t[0]==0):
        return 0
    A_t, I_t = [[g for g in list(influences.index[np.argwhere(np.array(regulated_t)==v).flatten().tolist()]) if (g in state.index)] for v in [1,-1]]
    if (len(A_t)==0 and len(I_t)==0):
        return 0#-float("inf")
    if (len(A_t)>0):
        P_At = state.loc[A_t].values
        A = np.mean(P_At)
        a = np.var(P_At)/len(A_t)
    else:
        A = 0
        a = 0
    if (len(I_t)>0):
        P_Bt = state.loc[I_t].values
        B = np.mean(P_Bt)
        b = np.var(P_Bt)/len(I_t)
    else:
        B = 0
        b = 0
    assert a+b>0
    score = (A-B)/np.sqrt(a+b)
    return score
def influence_t(t):
    return influence_t_per_state(t, initial_states)

influences_t = pd.DataFrame([[influence_t(g)] for g in all_genes], index=all_genes, columns=["Influence"]).sort_values(by="Influence", ascending=False)
influences_t.to_csv(vfolder+"influences_t_values.rnk", sep="\t")
## ORA
I = influences_t.loc[influences_t["Influence"].abs()>thres["Influence"]].index
influence_list = pd.DataFrame(I, index=range(len(I)))
influence_bklist = pd.DataFrame(influences_t.index, index=range(influences_t.shape[0]))
influence_list.to_csv(vfolder+"influences_t_list.txt", sep="\t", header=None, index=None)
influence_bklist.to_csv(vfolder+"infleunces_t_bkgnd.txt", sep="\t", header=None, index=None)

cytoscape_name = root_folder+name+"/cytoscape/solution.sif default node.csv"
infl_HIPP_epileptic_name = root_folder+name+"/spread_values.rnk"

## Process Cytoscape results (Cytoscape)
cytoscape = pd.read_csv(cytoscape_name)
cytoscape.index = cytoscape[["shared name"]].values.flatten().tolist()
cytoscape["Degree"] = cytoscape["Indegree"]+cytoscape["Outdegree"]
cytoscape = cytoscape[["ControlCentrality","MDS","Degree","Outdegree"]].astype(float)
for g in [gg for gg in all_genes if (gg not in cytoscape.index)]:
    cytoscape.loc[g] = 0

def list2ranking(measure, thres=thres["CC"]):
    centrality = cytoscape[[measure]].astype(np.float64).sort_values(by=measure, ascending=False)
    centrality.index = list(map(str,centrality.index))
    centrality.to_csv(vfolder+measure.lower()+"_values.rnk", sep="\t", header=None)
    ## ORA
    C = centrality.loc[centrality[measure]>thres].index
    centrality_list = pd.DataFrame(C, index=range(len(C)))
    centrality_bklist = pd.DataFrame(centrality.index, index=range(centrality.shape[0]))
    centrality_list.to_csv(vfolder+measure.lower()+"_list.txt", sep="\t", header=None, index=None)
    centrality_bklist.to_csv(vfolder+measure.lower()+"_bkgrnd.txt", sep="\t", header=None, index=None)
    return centrality

def im2ranking(nm, thres=thres["Spread"]):
    if (nm == "random"):
        fn = infl_random_name
    elif (nm == "hipp"):
        fn = infl_HIPP_name
    else:
        fn = infl_HIPP_epileptic_name
    ranking = pd.read_csv(fn, sep="\t", index_col=0, header=None)
    ranking.index = list(map(str,ranking.index))
    ranking.columns = [nm[0].upper()+nm[1:]+"_Spread"]
    ranking.to_csv(vfolder+nm+"_values.rnk", sep="\t", header=None)
    ## ORA lists
    R = ranking.loc[ranking[ranking.columns[0]]>thres].index
    ranking_list = pd.DataFrame(R, index=range(len(R)))
    ranking_bklist = pd.DataFrame(ranking.index, index=range(ranking.shape[0]))
    ranking_list.to_csv(vfolder+nm+"_list.txt", sep="\t", header=None, index=None)
    ranking_bklist.to_csv(vfolder+nm+"_bkgrnd.txt", sep="\t", header=None, index=None)
    return ranking

## Create ranking with ControlCentrality (Cytoscape)
controlcentrality = list2ranking("ControlCentrality")

## Create ranking with MDS (minimum driver set)
mds = list2ranking("MDS")

## Create ranking with degree (Cytoscape)
degree = list2ranking("Degree")

## Create ranking with outdegree (Cytoscape)
outdegree = list2ranking("Outdegree")

## Create ranking with HIPP_epileptic_Spread (Spread)
hipp_epileptic_im = im2ranking("hipp_epileptic")
hipp_epileptic_im.columns = ["Spread"]

## Create ranking with pLI scores
## https://gnomad.broadinstitute.org/downloads#v2-lof-curation-results
#intolerance_url = "gs://gcp-public-data--gnomad/legacy/exacv1_downloads/release0.3.1/cnv"
## https://cloud.google.com/storage/docs/gsutil_install
## sudo apt-get install -y apt-transport-https ca-certificates gnupg
## echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
## curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
## sudo apt-get update && sudo apt-get install google-cloud-sdk (unfinished)
##
## sudo apt-get install gcc python3-dev python3-setuptools libffi-dev
## sudo apt-get install python3-pip
## python3 -m pip install gsutil
ffolder = ddfolder+"gnomad/"
ffile = "gnomad.v2.1.1.lof_metrics.by_gene.txt"
sb.call("mkdir -p "+ffolder, shell=True)
sb.call("wget -nc -c -q -O "+ffolder+ffile+".bgz https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/"+ffile+".bgz", shell=True)
sb.call("gunzip -c "+ffolder+ffile+".bgz > "+ffolder+ffile, shell=True)
## https://www.biostars.org/p/9475852/

pli = pd.read_csv(ffolder+ffile, index_col=0, sep="\t")[['pLI']].loc[all_genes].astype(np.float64).sort_values(by="pLI", ascending=False)
pli.index = list(map(str,pli.index))
pli.to_csv(vfolder+"pli_values.rnk", sep="\t", header=None)
P = pli.loc[pli["pLI"]>thres["pLI"]]
assert (P.values>thres["pLI"]).all()
pli_list = pd.DataFrame(P.index, index=range(len(P.index)))
pli_bklist = pd.DataFrame(pli.index, index=range(pli.shape[0]))
pli_list.to_csv(vfolder+"pli_list.txt", sep="\t", header=None, index=None)
pli_bklist.to_csv(vfolder+"pli_bkgrnd.txt", sep="\t", header=None, index=None)

## EDS
ffolder=ddfolder+"EDS/"
ffile = "mmc2.xlsx"
eds_url="https://www.cell.com/cms/10.1016/j.ajhg.2020.01.012/attachment/16ff2f31-e4aa-46b0-a384-a242843ac763/"
sb.call("mkdir -p "+ffolder, shell=True)
sb.call("wget -nc -c -q -O "+ffolder+ffile+" '"+eds_url+ffile+"'", shell=True)
## pip install xlsx2csv
if (not os.path.exists(ffolder+ffile.split("xlsx")[0]+"csv")):
    sb.call("xlsx2csv "+ffolder+ffile+" > "+ffolder+ffile.split("xlsx")[0]+"csv", shell=True)
mmc = pd.read_csv(ffolder+ffile.split("xlsx")[0]+"csv", sep=",", index_col=0)
probes = list(mmc.index)
if (not os.path.exists(ffolder+"matches.csv")):
    matches = request_biodbnet(probes, from_="Ensembl Gene ID", to_="Gene Symbol", taxonId=9606, chunksize=500)
    matches.to_csv(ffolder+"matches.csv")
matches = pd.read_csv(ffolder+"matches.csv", index_col=0)
matches["Ensembl"] = matches.index
## Using GeneCards
manual_di = {
        'ATP5G3': "ENSG00000154518",
        'CACNA1C': "ENSG00000151067",
        'CHML': "ENSG00000203668",
        'FAM49A': "ENSG00000197872",
        'GUCY1B3': "ENSG00000061918",
        'PAK7': "ENSG00000101349",
        'PDE10A': "ENSG00000112541",
        'SEPT6': "ENSG00000125354",
        'TMEM35': "ENSG00000126950",
}
missing_genes = matches.loc[[manual_di[k] for k in manual_di]]
missing_genes.index = [k for k in manual_di]
matches.index = matches["Gene Symbol"]
matches = matches.loc[[g for g in all_genes if (g in matches.index)]]
matches = pd.concat((matches, missing_genes), axis=0)
matches = matches.loc[all_genes]
mmc2 = mmc.loc[matches["Ensembl"]]
mmc2.index = all_genes

def mmc2ranking(nm, thres=thres["pLI"]):
    erk = mmc2[[nm]].astype(np.float64).sort_values(by=nm, ascending=False)
    erk.to_csv(vfolder+nm.lower()+"_values.rnk", sep="\t", header=None)
    ## convert to ORA gene lists
    E = erk[erk>thres].index
    erk_glist = pd.DataFrame(E, index=range(len(E)))
    erk_bklist = pd.DataFrame(erk.index, index=range(erk.shape[0]))
    erk_glist.to_csv(vfolder+nm.lower()+"_list.txt", sep="\t", header=None, index=None)
    erk_bklist.to_csv(vfolder+nm.lower()+"_bckgrnd.txt", sep="\t", header=None, index=None)
    return erk 

eds = mmc2ranking("EDS")
csv_nucl = mmc2ranking("ActivityLinking_Conserved_nt_count")
csv_nucl.columns = ["Conserved_nt_counts"]
rvis = mmc2ranking("RVIS")
## RVIS > 0: tolerance// RVIS < 0: intolerance
##http://epilepsygenetics.net/2013/10/06/mutation-intolerance-why-some-genes-withstand-mutations-and-others-dont/
rvis = -rvis
rvis.columns=["-RVIS"]

## Compare rankings
all_rankings = [controlcentrality, outdegree, hipp_epileptic_im, pli, eds, rvis, influences_t]
##all_rankings = [controlcentrality,degree,outdegree,hipp_epileptic_im,pli]
###all_rankings = [controlcentrality,hipp_epileptic_im,pli,rvis]

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib_venn import venn3

method="spearman"#"pearson"

fig, ax = plt.subplots(figsize=(12,12))
ffsize=35
##ffsize=35
###ffsize=47
res = pd.concat(tuple([x.loc[all_genes] for x in all_rankings]), axis=1)
M = res.corr(method=method)
M.columns = ["Control\nCentrality", "Outgoing", "Spread", "pLI", "EDS", "-RVIS", "TIF"]
#M.columns = ["Control\nCentrality","Outdegree","Spread","pLI","EDS","-RVIS","Influence"]
##M.columns = ["Control\nCentrality","Degree","Outdegree","Spread","pLI"]
###M.columns = ["Control\nCentrality","Influence","pLI","-RVIS"]
M.index = M.columns
im = ax.imshow(M, vmin=0, vmax=1, interpolation='None')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(40)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        text = ax.text(j, i, np.round(M.iloc[i, j],2), ha="center", va="center", color="w" if (M.iloc[i,j]<0.5) else "k",fontsize=ffsize)
ax.set_xticks([0]+list(range(M.shape[1])))
ax.set_xticklabels([""]+list(M.columns), rotation=45,fontsize=ffsize)
ax.set_yticks([0]+list(range(M.shape[0])))
ax.set_yticklabels([""]+list(M.index), rotation=0,fontsize=ffsize)
plt.savefig(pfolder+"correlation_heatmap.pdf", bbox_inches="tight")
plt.close()

## Values for HIPP Spread/pLI for genes
res_hipp_pli = pd.concat(tuple([x.loc[all_genes] for x in all_rankings]), axis=1)
imname = "Spread" 
res_hipp_pli = res[["ControlCentrality", imname,"pLI"]]
ttt = 0.01
ffsize=40
ffsize=43
res_hipp_pli = res_hipp_pli.loc[res_hipp_pli[imname]>ttt].sort_values(by=imname, ascending=False)
res_hipp_pli.columns = ["CC", "Spread", "pLI"]
##res_hipp_pli.columns = ["CC","Infl.","pLI"]
#res_hipp_pli.index = ["*"+x if (mds.loc[x]["MDS"]==1) else x for x in list(res_hipp_pli.index)]
cm = ["Blues", "Oranges", "Greens", "Purples", "winter", "bone", 'Reds']
f, axs = plt.subplots(res_hipp_pli.columns.size, 1, gridspec_kw={'wspace': 0},figsize=(50,res_hipp_pli.shape[1]))
for i, (s, a, c) in enumerate(zip(res_hipp_pli.columns, axs, cm)):
    im = a.imshow(np.array([res_hipp_pli[s].values]), interpolation="None", cmap=c, aspect='auto')
    for u in range(res_hipp_pli.shape[0]):
        text = a.text(u, 0, np.round(res_hipp_pli[s].values.flatten().tolist()[u],3), ha="center",va="center",color="w" if ((s=="pLI" and res_hipp_pli[s].values.flatten().tolist()[u]>0.75) or (s in ["Spread","Infl."] and res_hipp_pli[s].values.flatten().tolist()[u]>0.03) or (s=="CC" and res_hipp_pli[s].values.flatten().tolist()[u]>8)) else "k",fontsize=ffsize)
    a.set_xticks([0,0])
    a.set_xticklabels(["",res_hipp_pli.index[i]])
    divider = make_axes_locatable(a)
    ca = divider.append_axes("right",size="2%",pad=0.)
    cbar = plt.colorbar(im,cax=ca)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(ffsize*2/3)
    a.set_yticks([0])
    a.set_yticklabels([s],fontsize=ffsize)
    if (i>res_hipp_pli.shape[1]-2):
        a.set_xticks([0]+list(range(res_hipp_pli.shape[0])))
        a.set_xticklabels([""]+list(res_hipp_pli.index), rotation=0,fontsize=ffsize)
    else:
        a.set_xticks([])
        a.set_xticklabels([])
plt.savefig(pfolder+"comparaison_values.pdf",bbox_inches="tight")
plt.close()

## Venn diagram for pLI>t, HIPP_Spread>0, ControlCentrality>0
im = hipp_epileptic_im
t_pLI,t_HIPPSpread,t_CC=1,0,0
pli_set = set(list(pli.loc[pli["pLI"]>=t_pLI].index))
hipp_im_set = set(list(im.loc[im[imname]>t_HIPPSpread].index))
controlcentrality_set = set(list(controlcentrality.loc[controlcentrality["ControlCentrality"]>t_CC].index))
print("pLI set = %d" % len(pli_set))
print(imname.lower()+" set = %d" % len(hipp_im_set))
print("cc set = %d" % len(controlcentrality_set))
Abc = len(pli_set.difference(hipp_im_set).difference(controlcentrality_set))
aBc = len(hipp_im_set.difference(pli_set).difference(controlcentrality_set))
ABc = len(pli_set.intersection(hipp_im_set).difference(controlcentrality_set))
abC = len(controlcentrality_set.difference(pli_set).difference(hipp_im_set))
AbC = len(pli_set.intersection(controlcentrality_set).difference(hipp_im_set))
aBC = len(controlcentrality_set.intersection(hipp_im_set).difference(pli_set))
ABC = controlcentrality_set.intersection(hipp_im_set).intersection(pli_set)
fig, ax = plt.subplots(figsize=(6,6))
venn3(subsets = (Abc, aBc, ABc, abC, AbC, aBC, len(ABC)), set_labels = (r'pLI$\geq %.2f$'%t_pLI, r'%s $> %d$' % (imname,t_HIPPSpread), r'CC $> %d$' % t_CC))
plt.title("Venn diagram of gene lists for measures pLI, "+imname+", CC")
plt.savefig(pfolder+"venn_diagram.png",bbox_inches="tight")
plt.close()

print(list(ABC))

## Compare the values of influence maximization and the size of the (weakly) connected component from this gene

## DFS on undirected network
def get_weakly_connected(edges, genes, gene):
    N = len(genes)
    adjacency = np.zeros((N,N))
    ## build DIRECTED adjacency matrix
    for g1,g2 in edges:
        ig1,ig2 = genes.index(g1),genes.index(g2)
        adjacency[ig1,ig2] = 1
    ## DFS
    to_visit = [gene]
    component = []
    while (True):
        node = to_visit.pop()
        if (node not in component):
            component.append(node)
        children = np.argwhere(adjacency[genes.index(node),:]==1).flatten().tolist()
        children = [genes[c] for c in children]
        to_visit = [child for child in children if (child not in component)]+to_visit
        if (len(to_visit)==0):
            break
    return component

def get_edges_from_solution(sol):
    edges = []
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
        edges += [[r,k] for r in regulators]
    return edges

with open(solution_folder+"solution.bnet", "r") as f:
    g = f.read().split("\n")
    genes = [x.split(" <- ")[0] for x in g]
    grfs = [x.split(" <- ")[-1] for x in g]
    edges = get_edges_from_solution({g: grfs[ig] for ig,g in enumerate(genes)})

im = hipp_im if (imname == "Hipp_Spread") else hipp_epileptic_im
hipp_genes_CC_size = [len(get_weakly_connected(edges, genes, g)) for g in im.index]
hipp_genes_CC = pd.DataFrame(hipp_genes_CC_size, index=im.index, columns=["CC size"])

from scipy.stats import pearsonr, spearmanr
test="spearmanr"
stat, p = eval(test)(list(im[im.columns[0]]), list(hipp_genes_CC["CC size"]))
print((r"Spearman $\rho$" if (test == "spearmanr") else r"Pearson's $r$")+" correlation %s (p=%s)" % (stat, p))

## ORA enrichment barplot
from glob import glob
enrichment_files=[x.split("/")[-1].split(".zip")[0] for x in glob(rfolder+"enrichments/ORA_*_HippE.zip")]
enrichment_files+=[x.split("/")[-1].split(".zip")[0] for x in glob(rfolder+"enrichments/ORA_*_pLI.zip")]
fsize=40
thres=0.05
ncategories=10

#' @param normalized Python bool
#' @return Axes of corplot
def plot_barplot(name_list, gene_list, nes_list, fdrs, overlap, order, enrichment, pvalues=[], thres=thres):
    '''Plots barplot to sum up enrichment results'''
    #colors  UP-FDR < thres   UP-FDR >= thres  DN-FDR < thres  DN-FDR >= thres
    colors = {'DN-FDR-': (0, 76/255., 117/255.), 'DN-FDR+': (0, 166/255., 1.), 'UP-FDR-': (166/255., 33/255., 0), 'UP-FDR+': (1, 51/255., 0)}
    plt.figure(figsize=(23, 11))
    if (enrichment_type == "ORA"):
        fdr_colors = [colors['UP-FDR-'] if (float(x) < thres) else colors['UP-FDR+'] for x in fdrs]
        color_mapper = list(map(str, fdr_colors))
        unique_colors = list(set(color_mapper))
        handle_makers, handles, handle_labels = [color_mapper.index(x) for x in unique_colors], [], []
        fdr_labels = [((r'BH-adjusted $p<%s$' if (not french) else r'$p$ corrigée $<%s$') % thres) if (float(x) < thres) else ((r'BH-adjusted $p\geq%s$' if (not french) else r"$p$ corrigée $\geq%s$") % thres) for x in fdrs]
        gene_list_ordered = gene_list
    if (enrichment_type == "GSEA"):
        fdr_colors, fdr_labels = [], []
        for i, x in enumerate(fdrs):
            head = "UP" if (nes_list[i] > 0) else "DN"
            tail = "-" if (float(x) < thres) else "+"
            fdr_colors.append(colors[head+"-FDR"+tail])
            fdr_labels.append("FDR "+("<" if (float(x) < thres) else ">=")+" "+str(thres)+" ("+head+")")
        color_mapper = list(map(str, fdr_colors))
        unique_colors = list(set(color_mapper))
        handle_makers, handles, handle_labels = [color_mapper.index(x) for x in unique_colors], [], []
        gene_list_ordered = []
        for i, ls in enumerate(gene_list):
            gene_order = list(sorted(range(len(ls)), key=lambda j : pvalues[i][j]))
            gene_list_ordered.append([ls[k] for k in gene_order])
    shifts = [(0.05, 0.41)]*len(name_list) if ("ORA" == enrichment_type) else [(0.05, 0.05) if (nes_list[i] > 0) else (np.min(nes_list)+0.25, 0.05) for i in order]
    for id_, i in enumerate(order):
        h = plt.barh(id_, nes_list[i], height=0.7, color=fdr_colors[i])
        plt.text(shifts[i][0], id_-0.15, name_list[i], fontsize=fsize, color="white" if (float(fdrs[i])<thres) else "black")
        genes = reduce(lambda x,y: x+", "+y, gene_list_ordered[i][:5])+(", ..." if (len(gene_list_ordered[i]) > 5) else "")
        if (i in handle_makers):
            handles.append(h)
            handle_labels.append(fdr_labels[i])
    if (False):
        if (any([x > (4 if ("ORA" == enrichment_type) else 2) for x in nes_list])):
            h = plt.plot([(4 if ("ORA" == enrichment_type) else 2)]*2, [-1, len(name_list)], "r--")
            handles.append(h[0])
        if (any([x < -2 for x in nes_list])):
            if ("GSEA" == enrichment_type):
                h = plt.plot([-2]*2, [-1, len(name_list)], "r--")
                handles.append(h[0])
        handle_labels += ["High "+("NES" if ("GSEA" == enrichment_type) else "ER")+" (|.| > "+("4" if ("ORA" == enrichment_type) else "2")+")"]
    plt.yticks([])
    plt.xticks(fontsize=fsize)
    plt.xlabel(("NES" if ("GSEA" == enrichment_type) else ("Odds Ratio" if (not french) else "Rapport des chances")),fontsize=fsize)
    plt.ylabel(r"Annotations by BH-adjusted $p$-value" if (not french) else ("Annotations classées par "+(r'$p$')+"-valeur")+"\ncorrigée croissante",fontsize=fsize)
    plt.legend(handles, handle_labels, loc='upper right', bbox_to_anchor=(1, 1), fontsize=fsize)
    return h

from functools import reduce
plot_fname="barplot"
for ef in enrichment_files:
    enrichment_type = ef.split("_")[0]
    sb.call("rm -rf "+enrichment_type+"/run/", shell=True)
    sb.call("mkdir -p "+enrichment_type+"/run/", shell=True)
    ## Only plot TOP-ncategories FDR
    print(ef)
    sb.call("cp "+rfolder+"enrichments/"+ef+".zip "+enrichment_type+"/run/", shell=True)
    sb.call("unzip -q -u "+enrichment_type+"/run/"+ef.split("/")[-1]+" -d "+enrichment_type+"/run/", shell=True)
    result = glob(enrichment_type+"/run/enrichment_results_*.txt")
    assert len(result) == 1
    res_df = pd.read_csv(result[0], sep="\t")
    try:
        name_list = list(res_df["description"])
    except:
        name_list = list(res_df["geneSet"])
    gene_list = [x.split(";") for x in list(res_df["userId"])]
    fdrs = list(res_df["FDR"])
    order = list(sorted(range(len(list(res_df.index))), key=lambda i : fdrs[i]))[:ncategories]
    fname_header = plot_fname+"-"+ef
    def round_fdr(f, n=2):
        if ("e" in str(f)):
            subzero = format(float(f), '.'+str(int(str(f).split("e-")[-1])+n)+'f')
        else:
            if (f == 1):
                return "1.00"
            else:
                subzero = float(f)
        subzero = str(subzero).split("0.")[-1]
        stop, nb_subzeros, i = False, 0, 0
        while (not stop):
            if (int(subzero[i]) == 0):
                nb_subzeros += 1
            i += 1
            stop = (int(subzero[i]) != 0)
        nonzero_subzeros = subzero[i:]
        if (len(nonzero_subzeros) > n):
            nonzero_subzeros = nonzero_subzeros if (int(nonzero_subzeros[n]) < 5) else str(int(nonzero_subzeros[:n])+1)
        nonzero_subzeros = nonzero_subzeros[:n]
        return "0."+("0"*nb_subzeros)+nonzero_subzeros
    if (enrichment_type == "ORA"):
        nes_list = list(map(lambda x : round(float(x), 2), list(res_df["enrichmentRatio"])))
        overlap = [r"$p\approx$"+str(round_fdr(res_df.loc[x]["FDR"])) for x in list(res_df.index)]
        pvalues = [1.]*len(nes_list)
    if (enrichment_type == "GSEA"):
        nes_list = list(map(lambda x : round(float(x), 2), list(res_df["normalizedEnrichmentScore"])))
        overlap = [str(res_df.loc[x]["leadingEdgeNum"])+"/"+str(res_df.loc[x]["size"])+", fdr="+str(round_fdr(res_df.loc[x]["FDR"])) for x in list(res_df.index)]
        data = df[df.columns[0]]
        ## Select a single value for multiple same indices
        pvalues = [[data.iloc[list(data.index).index(g)] for g in ls] for ls in gene_list]
        ## Gets enrichment plots for significantly enriched gene categories in GSEA (FDR < 0.05 & |NES| > 2)
        if (not os.path.exists(fname_header)):
            sb.call("mkdir "+fname_header, shell=True)
        from numpy import argwhere, array
        significantly_enriched = argwhere(array([int(abs(float(nes_list[i])) > 2 and float(fdrs[i]) < 0.05) for i in range(len(nes_list))]) == 1).flatten().tolist()
        name_SE = [list(res_df["geneSet"])[i].split(":") for i in significantly_enriched]
        name_SE = [x[0]+"_"+x[1] for x in name_SE]
        for i, name in enumerate(name_SE):
            ls = glob(enrichment_type+"/run/Project_*_GSEA/"+name+".png")
            assert len(ls) == 1
            sb.call("cp "+ls[0]+" "+fname_header+"/"+reduce(lambda x,y : x+"_"+y, name_list[significantly_enriched[i]].split(" "))+".png", shell=True)
    sb.call("mkdir -p "+enrichment_type+"/run/", shell=True)
    plot_barplot(name_list, gene_list, nes_list, fdrs, overlap, order, "/".join(fname_header.split("/")[-1].split("-")[:2]), pvalues=pvalues)
    plt.savefig(pfolder+fname_header+".pdf", bbox_inches="tight")
    sb.call("rm -rf "+enrichment_type+"/", shell=True)
