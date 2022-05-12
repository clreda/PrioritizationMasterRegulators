#coding;: utf-8

from multiprocessing import cpu_count
njobs=max(1, cpu_count()-2)

## TO CHANGE
##############################################
root="/media/kali/1b80f30d-2803-4260-a792-9ae206084252/Code/M30/"
root_folder = root+"PrioritizationMasterRegulators/"
##############################################

## Global variables (to modify if needed)
seed_number = 0 # replicability
data_folder = root+"data/"

## Model-related files
## TO CHANGE
###############################################################
dataset_folder = root_folder+"Refractory_Epilepsy_Data/"
disease_cid = "C0014544"
disease_name = "Epilepsy"
disease_cell_line = ["NPC","SHSY5Y"]
taxon_id=9606 # human
## List of HUGO gene symbols (one / line) 
path_to_genes = dataset_folder+"S1_Delahayeetal2016_M30genes.txt"
path_to_initial_states = dataset_folder+"EMTAB3123.csv"
name = "refractory_epilepsy"
###############################################################

file_folder = root_folder+name+"/FILES/"
plot_folder = root_folder+name+"/PLOTS/"
solution_folder = root_folder+name+"/"

## PARAMETERS ## TO CHANGE IF NEEDED, especially those with the *TO CHANGE* annotation
args = {
        ## Experiments
        "pert_types" : ["trt_sh", "trt_oe", "trt_xpr"], # Types of perturbations
        "cell_lines" : disease_cell_line,               # Selection of cell lines (empty list if automatic selection in LINCS)
        "selection" : "distil_ss",                  # Criterion to maximize when selecting experiments
        "thres_iscale" : 0,                             # "Interference scale" to select experiments
        ## Visualization
        "plot_it" : True,                               # Generation of plots
        ## GRN
        "score_STRING": 400, #*TO CHANGE*               # Threshold on STRING score for extracted PPI
        "full_graph": False,                            # Gives as input to the solver the complete graph (of N^2 edges for N nodes)
        "direct_graph": False,                           # Gives as input to the solver a directed graph using CBDN
        "cor_method": "pearson",                        # Use as correlation method to sign edges
        "exact": False,                                 # If set to True, no selection of edges in the input graph is done ( = all edges preserved in the input network )
        ## Generation of solutions
        "limit":1,                                # Limit on the number of solutions to generate
        "niterations": 25, #*TO CHANGE*                   # Number of iterations of the procedure
        ## Edge filtering
        "beta": 1,                                      # Power used to filter out edges based on their correlation
        "tau": 0.4,  #*TO CHANGE*                      # Threshold used to filter out edges based on their correlation
        "use_diverse": True,                            # For reproductibility, set it to False (BoneSiS parameter)
        "add_other_genes": False,                       # If a set of genes is not present in the correlation matrix, should they be added to the input network (potential edges go from these genes to all the genes in the network)? 
        ## Processing of experimental data
        "bin_thres": 0.265,  #*TO CHANGE*                # Threshold on the binarization of expression data in [0,0.5] 
        "bin_method": "binary",                         # Binarization method
        "use_CD": False,                  # Use CD-binarization for final signatures
        "use_subset_exp": None,             # Use subset of experiments
        # For the desirability criterion: DS: network density, CL: clustering coefficient, Centr: centralization, GT: heterogeneity, RD: redundancy
        "weights": {"DS":3, "CL": 3, "Centr":3, "GT": 1},
        }

disgenet_args = {
        "min_score":0, 
        "min_ei":0, 
        "min_dsi":0.25, 
        "min_dpi":0,
}

maboss_params = {
    'sample_count': 1000,
    'use_physrandgen': 0,
    'thread_count': njobs,
    'max_time': 50,
    'time_tick': 1,
}

im_params = {
    "k":1, #*TO CHANGE*
    "njobs": min(5,njobs),
    "state_window": 100,
    "method": "greedy",
}

for a in args:
    globals()[a] = args[a]
