# Tree construction

The viral genome sequences were mapped to the reference sequence of Wuhan-Hu-1 (GenBank accession number: NC_045512.2) using Minimap2 v2.17 (Li, 2018) and subsequently converted to a multiple sequence alignment according to the GISAID phylogenetic analysis pipeline (https://github.com/roblanf/sarscov2phylo). The alignment sites corresponding to the 1–265 and 29674–29903 positions in the reference genome were masked (i.e., converted to NNN). Alignment sites at which >50% of sequences contained a gap or undetermined/ambiguous nucleotide were trimmed using trimAl v1.2 (Capella-Gutierrez et al., 2009). Phylogenetic tree construction was performed via a three-step protocol: i) the first tree was constructed; ii) tips with longer external branches (Z score > 4) were removed from the dataset; iii) and the final tree was constructed. Tree reconstruction was performed by RAxML v8.2.12 (Stamatakis, 2014) under the GTRCAT substitution model. The node support value was calculated by 100 times bootstrap analysis.

## Contents:
*  **tree_construction/batch_pipeline.sh:** the main bash script
*  **tree_construction/replace_N2hyphen.py:** a python script replacing "N" (undetermined nucleotide) to "-" (gap)
*  **tree_construction/detect_outlier_OTU.R:** an R script removing OTUs with longer external branches
*  **tree_construction/NC_045512.fas:** the reference genome sequence of SRAS-CoV-2 (Wuhan-Hu-1, a lineage B isolate)
