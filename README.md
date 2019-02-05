# nmd_nar
Autoregulation via AS-NMD

# RBP autoregulation via NMD
This repository contains scripts and pipelines for "Integrative transcriptomic analysis suggests new autoregulatory splicing events coupled with nonsense-mediated mRNA decay" paper by Pervouchine et al.

## Installation
After git clone https://github.com/pervouchine/nmd_nar, enter the directory and run 'make all'. During the execution of make, it will download the processed data and produce all figures and tables
Dependencies:
# R packages: ggplot2, ggrepel, plyr, dplyr, ggsignif
# bedtools v2.25.0
# bedToBigBed
# devtools
# jdstorey/qvalue (install_github("jdstorey/qvalue"))

## Data structure
There are three data repositories in this project
1. eCLIP
2. shRNA-KD followed by RNA-seq
3. shRNA-codepletion of UPF1/XRN1 and SMG6/XRN1
After unpacking, they will be stored in data/ directory
