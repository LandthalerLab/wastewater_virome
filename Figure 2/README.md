This part describes the workflow to generate the heatmaps in Figure 2. Based on the kaiju metagenomics pipeline output, the process is as following:

1. The kaiju pipelines identifies taxonomy IDs found in the sequencing data. The first step is to download all the fasta files belonging to the taxonomy IDs in the kaiju output file from a specific datasets. These output files for the Berlin wastewater RNA and DNA data are on the GEO entry (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228220). The kaiju output files of the three other datasets, California (…), India (Stockdale et al. 2023), and Baseline (Nieuwenhuijse et al. 2020) will be available there soon. This is done using the bash script mapvirus_part1.sh, which requires next to the location of the kaiju output the NCBI entrez direct tools.