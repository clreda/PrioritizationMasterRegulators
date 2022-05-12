#!/bin/bash

conda create -n application_envir -c colomoto pymaboss python=3.9.6
conda activate application_envir
python3 -m pip install git+https://github.com/bioasp/bonesis.git@64e88178816f86ff112cd9c9e423cf40e5029c5c
python3 -m pip install cmapPy==4.0.1 qnorm==0.8.0 matplotlib-venn==0.11.6 seaborn==0.11.2 patsy==0.5.2 statsmodels==0.13.2
#python3 -m pip install tqdm==4.62.3
conda deactivate
