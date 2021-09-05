# TMI - ttest module identification
Codebase for TMI and GLADIATOR analyses.
TMI is a new method for inferring cell type-specific enhancer-promoter (EP) mediating proteins

# Getting TMI R objects (*.rds) and data
You can access the data in the following paths:
1) /home/gaga/tomhait/projects/ePPI/rds
2) /home/gaga/tomhait/projects/ePPI/data
3) Raw data - /home/gaga/html/ct-focs/

# Dependencies
Please make sure you have the following pre-installed R packages:
* parallel
* igraph
* httr
* stringi
* dplyr
* plyr
* BSgenome.Hsapiens.UCSC.hg19
* topGO
* seqinr
* clusterProfiler

Also, make sure you have the MEME suite installed on Linux platform:
https://meme-suite.org/meme/doc/download.html

# Running TMI/GLADIATOR to infer EP mediating proteins
Please note that you can run TMI/GLADIATOR only under linux platform

# To run GLADIATOR:
python -W ignore GLADIATOR.py -o data/GLADIDATOR_modules.txt -n data/Interactome.tsv -s data/SeedPS.tsv -p data/PhenSimMat.tsv

The file data/GLADIDATOR_modules.txt, data/Interactome.tsv, data/SeedPS.tsv and data/PhenSimMat.tsv are generated using the MAIN_GLADIATOR_alg.R script.

# Within the scripts folder you will find the relevant scripts for the analysis:

Filename | Description  
--- | ---  
MAIN_TF_pair_identification.R | Identification of TF-TF pairs in EP links.
FUNCTIONS_motif_finding.R | Auxiliary functions for motif finding.
MAIN_String_analysis.R | Calculates normalize betweenness values per cell type.
MAIN_GLADIATOR_alg.R | Creates the input files for GLADIATOR program.
MAIN_TMI_alg.R | Applies the TMI method on the obtained betweenness values from MAIN_String_analysis.R.
MAIN_enrichment_analysis.R | Performs GO enrichment analyses on predicted cell type-specific PPI modules.

Please note that the generated code is a line by line running code.

Additional instructions: sthait at gmail dot com; amitlevon at gmail dot com
