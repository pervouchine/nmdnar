# nmd_nar
Autoregulation via AS-NMD

# RBP autoregulation via NMD
This repository contains scripts and pipelines for "Integrative transcriptomic analysis suggests new autoregulatory splicing events coupled with nonsense-mediated mRNA decay" paper by Pervouchine et al.

## Installation
After git clone https://github.com/pervouchine/nmd_nar, enter the directory and run 'make all'. During the execution of make, it will download the processed data and produce all figures and tables

Dependencies:
* R packages: ggplot2, ggrepel, plyr, dplyr, ggsignif
* bedtools v2.25.0
* bedToBigBed
* devtools
* jdstorey/qvalue (install_github("jdstorey/qvalue"))

## Data structure
There are three data repositories in this project
1. eCLIP
2. shRNA-KD followed by RNA-seq
3. shRNA-codepletion of UPF1/XRN1 and SMG6/XRN1

They have been pre-processed by IPSA pipeline (https://github.com/pervouchine/ipsa) and deposited in compressed form at our local data storage (http://arkuda.skoltech.ru/~dp/papers/nmd/data.tar.gz). The pipeline will automatically download and unpack it. After unpacking, the data will be stored in data/ directory.

## Execution
When in the main directory, type 'make all'. It will first download the data and then run the analysis steps on them. For the explanation of pipeline steps please refer to the makefile, which also has comments explaining each step among necessary parts such as rules and recipies.

## Using the track hub
To use the UCSC Genome Browser track hub system, go the page https://genome.ucsc.edu/cgi-bin/hgHubConnect and insert the following link https://raw.githubusercontent.com/pervouchine/nmd_nar/master/hub/hub.txt
In case if it doesn't work for some reason, open hub/hub.txt in this repository, press "Raw" button and copy and paste the URL into the Genome Browser field. Note that the track hub is currently available only for GRCh37 assembly of the human genome and the latest GENCODE transcript annotation (v.19) available for that release.
