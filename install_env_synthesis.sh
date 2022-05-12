#!/bin/bash

conda create -n synthesis_envir -c potassco clingo python=3.8.5
conda activate synthesis_envir
apt-get install -y graphviz
python3 -m pip install git+https://github.com/bioasp/bonesis.git@64e88178816f86ff112cd9c9e423cf40e5029c5c
python3 -m pip install cmapPy==4.0.1 matplotlib==3.3.4 scikit-learn==0.24.2
python3 -m pip install qnorm==0.5.1 maboss==0.8.1 #tqdm==4.62.3 mygene=3.2.2
python3 -m pip install openpyxl==3.0.9
## for profile_binr (not used here)
#R -e "for (pkg in c('mclust', 'diptest', 'moments', 'magrittr')) install.packages(pkg);"
#R -e "for (pkg in c('tidyr', 'dplyr', 'tibble', 'bigmemory')) install.packages(pkg);"
#R -e "for (pkg in c('doSNOW', 'foreach', 'glue')) install.packages(pkg);"
conda deactivate
