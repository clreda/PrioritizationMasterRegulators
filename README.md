# Prioritization of master regulators through a Boolean network
(c) Clémence Réda, 2022.

Due to the presence of copyrighted databases, the license for this code is [Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## Environment

Create two Conda environments using the two scripts *install_env_application.sh* (environment "application\_envir") and *install_env_synthesis.sh* (environment "synthesis\_envir")

## Data (refractory epilepsy)

Import the initial states from Mirza et al., 2017 and the M30 genes from Delahaye-Duriez et al., 2016

```bash
conda activate synthesis_envir
python3 download_Refractory_Epilepsy_Data.py
conda deactivate
```

## Building a Boolean network

You need to register to the [LINCS L1000 database](https://clue.io/developer-resources#apisection) and the [DisGeNet database](https://www.disgenet.org/) and write up the corresponding credentials and API keys to files *utils/credentials/credentials_LINCS.py* and *utils/credentials_DISGENET*. Write up your name in *utils/credentials/credentials_STRING.py*.

Then check out the *params.py* file, where some parameters might need to be changed (documentation therein). Then run the following commands

```bash
conda activate synthesis_envir
python3 iter_main.py ## takes some time
python3 model_selection.py ## final model selection and network robustness plots
conda deactivate
```

The final solution is *solution.bnet* (with the initial values in *params.py*, this file is located at *refractory_epilepsy/*).

## Detection of master regulators

Check out *params.py* to change parameter values, if needed.

```bash
conda activate application_envir
python3 application_regulators.py
conda deactivate
```

The result files are named *application_regulators.csv*, *spread_values.csv* and *spread_values.rnk* (with the initial values in *params.py*, these files are located at *refractory_epilepsy/*).

## Network analysis with Cytoscape

Network analyses are performed with Cytoscape 3.8.0. You need to download the module CytoCtrlAnalyser (version 1.0.0). Then run

```bash
conda activate application_envir
python3 bonesis2cytoscape.py # converts graph files to Cytoscape-readable formats
conda deactivate
```

Compute at least MDS, ControlCentrality (in CytoCtrlAnalyzer) and Outdegree, Indegree (NetworkAnalyzer, built-in module in Cytoscape). You can run the following script to get source evidence annotations on edges from the STRING database

```bash
conda activate application_envir
python3 annotate_edges.py
conda deactivate
```
## Visualization of gene measures & nrichment analyses 

First run

```bash
conda activate application_envir
python3 compare_rankings.py
conda deactivate 
```

Then, use the online tool [WebGestalt](http://webgestalt.org/). Select the following parameters

- Organism of interest: Homo sapiens
- Method of interest: ORA
- Functional Database: disease > DisGeNet
- Select Gene ID Type: Gene symbol
- Upload Gene List: upload the list at *refractory_epilepsy/rankings/values/hipp_epileptic_list.txt*
- Upload User Reference Set: Gene symbol
- File and Select ID type: *refractory_epilepsy/rankings/values/hipp_epileptic_bkgrnd.txt*

In "Advanced parameters", leave default parameters, except for 
- Significance Level: FDR > 0.2

Download the corresponding .zip file and put it in *refractory_epilepsy/rankings/enrichments/* by using the same type of filename syntax as shown in the examples (note: if the enrichment files cannot be found, have a look at *compare_rankings.py*).

## Running on the set of 50 solutions

Copy the solutions in folder *refractory_epilepsy/25additionalsolutions/* in *refractory_epilepsy/*, and run (NOTE THAT THIS ERASES FOLDER *PLOTS/* AND FILES *SOLUTION.BNET*, *COLLAPSED_MODEL.CSV*)

```bash
conda activate synthesis_envir
python3 model_selection.py
conda deactivate
```

## Pull requests, issues, suggestions?

clemence.reda@inserm.fr
