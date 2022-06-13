#coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt

from params import *

vals = pd.read_csv(file_folder+"score_STRING_results.csv", index_col=0)
thres = args["score_STRING"]

plt.figure(figsize=(5,5))
plt.plot(list(vals.index), list(vals["#edges"]), "b+")
#plt.plot(list(vals.index), list(vals["#edges"]), "k-")
plt.plot([100,thres],[vals.loc[thres]["#edges"]]*2, "r--")
plt.plot([thres]*2, [100,vals.loc[thres]["#edges"]], "r--")
plt.text(0, vals.loc[thres]["#edges"], str(int(vals.loc[thres]["#edges"])), c="red")
plt.title(r"Number of edges in PPI depending on $\eta$")
plt.xlabel(r"Value of $\eta$")
plt.ylabel("Number of filtered edges in PPI")
plt.xlim(100,1000)
plt.ylim(0,8000)
plt.savefig(plot_folder+"selection_eta.pdf", bbox_inches="tight")
plt.close()
