# Pipelines and Scripts for Wyler, E. et al. 2024

Wyler, E., Lauber, C., Manukyan, A., Deter, A., Quedenau, C., Alves, L. G. T., ... & Landthaler, M. (2024). Pathogen dynamics and discovery of novel viruses and enzymes by deep nucleic acid sequencing of wastewater. Environment International, 108875.

**DOI**: [https://doi.org/10.1016/j.envint.2024.108875](https://doi.org/10.1016/j.envint.2024.108875)

This repository incorporates two pipelines for RNA and DNA samples generated from wastewater: 

* **Kaiju Pipeline** for taxonomy classification:
  - Duplicate removal of reads with [CD-HIT](https://sites.google.com/view/cd-hit).
  - Taxonomy classification and annotation with [Kaiju](https://kaiju.binf.ku.dk/).
  - Custom R scripts for summarizing annotated reads per sample.
  
* **CCTyper Pipeline** for predicting cas proteins:
  - Divide assemblies into chunks with [pyfasta](https://anaconda.org/bioconda/pyfasta).
  - CRISPR-Cas gene detection with [CRISPRCasTyper](https://anaconda.org/russel88/cctyper).

In addition to preprocessing pipelines, we also provide access to custom scripts for the downstream analysis for both taxonomy classification and Cas protein detection:

* **Metagenomic Analysis** on taxonomically classified reads:
  - taxonomy ranks and lineage with [taxizedb](https://cran.r-project.org/web/packages/taxizedb/index.html). 
  - Processing and PCA of taxonomy counts.
  - Visualize heatmaps with [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html).
  
* **Downstream Analysis** on reads associated with CRISPR-Cas genes:
  - Clustering open reading frames (ORFs) with [CD-HIT](https://sites.google.com/view/cd-hit).
  - Aligning ORFs to NR database in NCBI with [rBLAST](https://github.com/mhahsler/rBLAST).
  - Protein embeddings of ORFs using [ProtTrans](https://github.com/agemagician/ProtTrans).
