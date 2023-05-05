# Arabidopsis_ETI_transcriptome_response

R scripts and data sets used in the study. Run script 1 before script 2 as script 2 uses some intermediate .RData files generated by script 1.

1. “Script_Ed-AvrRpt2_ms_230503f.r” is a script to generate all results including all figures, except flg22-PTI-related results and Fig 6.

Input 1: “final_selected_ETI_libraries_gene_counts.R_input_matrix.txt” is the RNA-seq read count data set from NCBI-GEO GSE196892. This file is essential.

Inputs 2 and 3: “Yang.2020.circadian.159584_1_supp_509367_q8lbct.csv” and “Yang.2020.diurnal.159584_1_supp_509370_q8lbct.csv” are the gene lists including the model fit parameters for the circadian and the diurnal response genes from the supplemental tables in Yang et al. 2020: https://doi.org/10.1074/jbc.RA120.013513. They were converted from .xlsx to .csv. These files are required only for the part generating figures for Text S1.

Outputs include: Figs 1, 2, 3, 4, 6, S1, S2, S3, and S6. Table S1. Multiple intermediate .RData files


2. “Script_flg22.230503f.r” is a script to generate flg22-PTI-related results and Fig 6. 

Input 1: “flg22_raw_gene_counts.RData” includes the RNA-seq read count data set from NCBI-GEO GSE78735. The R matrix object name for the dataset in this file is “raw_data_corrected_sample_labels_nonzero_rows”. This file is essential.

Intermediate 1 and 2: Two “.RData” files generated by script 1, “NRAM.algorithm.w.mean.se.RData” and “ETI.gfa.sorted.genes.v230327.RData” are required to run this script.

Outputs include: Figs 5, S4, and S5. Table S2. Multiple intermediate .RData files


Input data sets are in ./data/

Outputs subfolder is required to run the scripts: ./outputs/
