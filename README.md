# Background

This repository’s code and data accompany the manuscript, _Codon usage and gene expression patterns associated with five-year survival rates in cancer patients_. In this study, novel quantitative techniques were performed on four cancer types from The Cancer Genome Atlas: Kidney Renal Cell Carcinoma, Head and Neck Squamous Cell Carcinoma, Lung Squamous Cell Carcinoma, and Brain Lower Grade Glioma. This code and data were used to identify new prognostic markers for these cancer types and explore the survival mechanisms.

This study utilizes both previously published and novel analytical techniques. Some examples of these techniques are Kaplan-Meier curves, the recently published betaMix method, and the differential variation and expression analysis method. More information about the techniques is provided in the manuscript.

# Dependencies

Download R’s dependencies using RStudio’s built in functionality or your preferred method from The Comprehensive R Archive Network (<https://cran.r-project.org/>), Bioconductor (<https://bioconductor.org/>), and/or GitHub (<https://github.com/>).

Our list of canonical transcripts is included in the publication’s supplemental files. Although, a different list of transcripts could be used. This data should be structured as a N by 3 matrix. Where N is the number of desired canonical genes.

Lists of oncogenes and tumor suppressor genes were obtained from the Ongene (<https://bioinfo-minzhao.org>), and TSgene (<https://bioinfo.uth.edu/TSGene>) databases respectfully. The codon usage of each oncogene and tumor suppressor gene was matched by its ENSG identifier to those of the canonical genes selected for this study.

Several different types of TCGA expression and mutation files are needed from the National Cancer Institute’s Genetic Data Commons for each cancer type’s analysis.

For the expression data download the tar balls: Spliced Transcripts Alignment to a Reference (.star), clinical data, biospecimen data, sample sheet, and metadata. After untarring the files, change sample sheet.tsv to expression_sample_sheet.csv and metadata.json to expression_metadata.json.

For the mutational data download the tar balls: Mutation Annotation Format (.maf), sample sheet, and metadata. After untarring the files, change sample sheet .tsv to mutation_sample_sheet.csv and metadata.json to mutation_metadata.json.

# Archived Dataset

1. Supplemental File 6.xlsx

Tables of codon changes associated with each cancer type’s observed SNVs. The MAF files were filtered to only include SNVs from canonical protein-coding genes and had the associated codon change provided. The tables also include supplemental information about each mutation. For example, the gene a mutation occurred in, the mutation’s impact, and whether the gene was DE, DV, or survival associated. KIRC, LGG, LUSC, and HNSC had 15,551, 21,378, 114,780, and 58,574 observed mutations respectively.

# Archived Scripts

The scripts are listed by the order they should be ran. Some of the scripts are computationally intensive and will have the best performance on a high-performance computing (HPC) platform. Computationally intensive scripts have “HPC” in their names, and the rest of the scripts have “Local” in their names.

1. HPC Gene Reference Codon Usage.R

This script calculates reference codon usage values for each canonical protein-coding gene.

1. HPC Transcriptomic-Weighted Codon Usage and Survival.R

This script calculated the transcriptomic weighted codon usage of the samples for each cancer type.

1. Local Preparation for Differential Gene Expression.R

This script prepares the data for differential expression analysis using the DVX package (available for download from <https://haim-bar.uconn.edu/software/dvx/>).

1. HPC TCGA Mutation Tables.R

This script creates TCGA SNV mutation summary tables for each cancer type.

1. HPC Figure Pre-processing.R

This script prepares data for post-processing and plotting.

1. Local Five-Year Survival.R

This script identifies codons with survival associated transcriptomic-weighted codon usage and/or age of diagnosis.

1. Local Analyses and Plots.R

This script computes RNA sequence frequency-weighted relative synonymous codon usage and differential gene expression/variability -weighted relative synonymous codon usage for each gene in each cancer. It also creates several plots and heatmaps to summarize these results and others from earlier in the pipeline.

# Authors

Nathan J. Clement <sup>1, +</sup>, Haim Bar <sup>2, +</sup>, Ryan C. Hunt <sup>1+</sup>, Douglas Meyer <sup>1</sup>, Anton A. Komar <sup>3</sup>, Michael DiCuccio <sup>4</sup>, and Chava Kimchi-Sarfaty <sup>1, \*</sup>

# Affiliations

<sup>1</sup> Hemostasis Branch 1, Division of Hemostasis, Office of Plasma Protein Therapeutics CMC, Office Therapeutic Products, Center for Biologics Evaluation and Research, Food and Drug Administration, Silver Spring, USA.

<sup>2</sup> Department of Statistics, University of Connecticut, Storrs, CT, USA.

<sup>3</sup> Center for Gene Regulation in Health and Disease, Department of Biological, Geological and Environmental Sciences, Cleveland State University, Cleveland, OH, USA.

<sup>4</sup> Rockville, MD, USA.

<sup>+</sup> Equal contributor

<sup>\*</sup> Corresponding author

Please direct questions to Nathan Clement ([nathan.clement@fda.hhs.gov](mailto:nathan.clement@fda.hhs.gov)) and Chava Kimchi-Sarfaty ([chava.kimchi-sarfaty@fda.hhs.gov](mailto:chava.kimchi-sarfaty@fda.hhs.gov)).

# Open-Source Disclaimer

This open-source software has been released as-is without expectation of support or future releases.
