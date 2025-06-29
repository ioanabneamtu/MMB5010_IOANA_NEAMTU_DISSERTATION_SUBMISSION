Differential Gene Expression Analysis - README
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Overview

	This project analyzed differential gene expression in Media vs. LPS treated monocytes, and in LPS-treated monocytes where TNF-alpha release was measured, using DESeq2. 
	This analysis was performed in RStudio Editor using R version 4.0.3 and through using specific versions of R packages from CRAN and Bioconductor to ensure reproducibility. 
	Below is the setup required to install the necessary libraries and run the analysis. 
	Reproducing the analysis requires not only the same libraries, but also the use of specific package versions and adherence to the documented workflow in this README. 
	By following the exact steps, input files, and software environment described—including R version 4.0.3 and listed package versions—users can reliably replicate the 
	original results.	
	The libraries used were crucial for the initial work, and loading these into the R environment will allow the same workflow to be reproduced seamlessly.
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
File with R code: Differential_Gene_Expression.Rmd
HTML Format Notebook Output: Differential_Gene_Expression.nb.html
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Working Directory Setup

    	Manual adjustment of file paths is needed to the locaion of the "DGE_Analysis.Rmd" file containing the R script.
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Libraries Employed (CRAN and Bioconductor)

	To install the required libraries, use the following commands:

	# Installing Bioconductor and required packages
	if (!requireNamespace("BiocManager", quietly = TRUE))
    		install.packages("BiocManager")

	# Installing necessary packages from Bioconductor (BiocManager version 3.17)
	BiocManager::install(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "EnhancedVolcano", "pheatmap"))

	# Installing necessary packages from CRAN
	install.packages(c("tidyverse", "ggrepel", "stringr"))
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Loading Libraries into the R Environment

	After installation, the following libraries can be loaded using library(). 

		Libraries used:
        		DESeq2 (1.40.2)
      			tidyverse (2.0.0)
        		ggrepel (0.9.5)
       			pheatmap (1.0.12)
        		EnhancedVolcano (1.18.0)
        		AnnotationDbi (1.66.0)
        		org.Hs.eg.db (3.17.0)
        		clusterProfiler (4.8.3)
        		stringr (1.5.1)
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
The Differential Gene Expression Pipeline described in the notebook follows these steps:

1. Environment Setup

    -Libraries from CRAN and Bioconductor were installed and loaded: DESeq2, tidyverse, ggrepel, pheatmap, EnhancedVolcano, clusterProfiler, stringr and annotation packages (org.Hs.eg.db, AnnotationDbi)

2. Data Loading and Preprocessing

    -Count data (LPSvsMEDIA_counts.csv, TNFalpha_counts.csv) and sample metadata (LPSvsMEDIA_sample_info.csv, TNFalpha_sample_info.csv) were loaded
    -Data cleaning was performed by removing rows with low counts

3. Building the DESeqDataSet and Normalization

    -A DESeqDataSet was created using the cleaned count data and metadata
    -Size factors were estimated using estimateSizeFactors() to normalize the data for sequencing depth variations

4. Quality Control

    -Exploratory data analysis included:
        -Principal Component Analysis (PCA) using variance-stabilized data to visualize sample groupings
        -Scree plots to assess variance explained by each principal component
        -Heatmaps using pheatmap to cluster samples based on their conditions
    -Data sparsity was visualized using plotSparsity()

5. Differential Expression Analysis

    -DESeq2 was used to fit a Negative Binomial General Linear Model to the count data
    -The DESeq() function performed the differential expression analysis, followed by Wald tests for significance
    -Raw and normalized counts were extracted, and total counts per sample were plotted

6. Dispersion and Shrinkage Estimation

    -Dispersion estimates were visualized using plotDispEsts()
    -Log2 Fold Change (LFC) Shrinkage was applied using lfcShrink() with various methods:
        -Normal shrinkage (default)
        -Ashr (adaptive shrinkage) - retained
        -Apeglm (Bayesian shrinkage)

7. Visualization of Results

    -MA plots were generated to visualize LFC against mean counts for each shrinkage method
    -A Volcano Plot was created using EnhancedVolcano to highlight significant differentially expressed genes
    -Histograms were plotted to show the distribution of raw and adjusted p-values for significant genes

8. Gene Ontology (GO) Enrichment Analysis

    -GO Enrichment was performed using clusterProfiler for both upregulated and downregulated genes, focusing on:
        -Biological Process (BP)
        -Cellular Component (CC)
        -Molecular Function (MF)
    -Results were visualized using bar plots for the top 20 enriched GO terms

9. Pathway Analysis

    -Gene symbols were mapped to ENSEMBL IDs using AnnotationDbi and org.Hs.eg.db
    -Upregulated and downregulated gene lists were exported for further analysis, such as STRING interaction network analysis  
    -STRING Results were provided using Permalinks accessible from the HTML Notebook

10. Saving Results

    -Differential expression results, including shrunken log2 fold changes and significant genes, were saved as text and CSV files
    -Plots, Figures and other Results Vizualisation are also provided

This structured approach enabled a comprehensive analysis of differential gene expression, from raw count data through quality control, statistical testing, and pathway analysis, providing insights into biological functions and pathways affected by the experimental conditions
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Input Files

	LPSvsMEDIA_counts.csv
	LPSvsMEDIA_sample_info.csv
	TNFalpha_counts.csv
	TNFalpha_sample_info.csv
	
Output Files

	DGE_DOWN_LPSvsMEDIA.csv
	DGE_DOWN_LPSvsMEDIA_gene_ensemblID.txt
	DGE_DOWN__LPSvsMEDIA_gene_symbols.txt
	DGE_DOWN_TNF-alpha.csv
	DGE_DOWN_TNF-alpha_gene_ensemblID.txt
	DGE_DOWN_TNF-alpha_gene_symbols.txt
	DGE_results_high_vs_low_LPS_TFN-alpha.txt
	DGE_results_LPS_vs_media_Type_ashr.txt
	DGE_UP_LPSvsMEDIA.csv
	DGE_UP_LPSvsMEDIA_gene_ensemblID.txt
	DGE_UP_LPSvsMEDIA_gene_symbols.txt
	DGE_UP_TNF-alpha_gene_symbols.txt
	DGE_UP_TNF-alpha.csv
	DGE_UP_TNF-alpha_gene_ensemblID.txt
	LPS_TNFalpha_normalised_counts.txt
	LPS_TNFalpha_unnormalised_counts.txt
	LPSvsMedia_normalised_counts.txt
	LPSvsMedia_PCA_loadings_PC1.csv
	LPSvsMedia_unnormalised_counts.txt
	sigDE_padj<=0.01_LPS_TNF-alpha_table_condition_Type_ashr.txt
	sigDE_padj<=0.01_LPS_vs_media_table_condition_Type_ashr.txt
	sig_DGE_results_LPS_TNF-alpha_annotated.txt
	sig_DGE_results_LPS_vs_media_annotated.txt
	TNFalpha_PCA_loadings_PC1.csv
	TNFalpha_PCA_loadings_PC2.csv
	 
Output Figues and Plots

	LPS_TNFalpha_Heatmap_responses-samples.png
	LPS_TNFalpha_MA_Plot_Ash_Shrinkage_withContrast_SexCovariate.png
	LPS_TNFalpha_PCA_Plot.png
	LPS_TNFalpha_Positively_Enriched_GO_CC.png
	LPS_TNFalpha_Positively_Enriched_GO_MF.png
	LPS_TNFalpha_Raw_Counts.png
	LPS_TNFalpha_Size_Factors.png
	LPS_TNFalpha_Variance_Scree_Plot.png
	LPS_TNFalpha_Volcano_Plot_DEGs.png
	LPSvsMedia_Dispersion_Estimates.png
	LPSvsMedia_Heatmap_conditions-samples.png
	LPSvsMedia_MA_Plot_Ash_Shrinkage_withContrast.png
	LPSvsMedia_MA_Plots_Shrinkage.png
	LPSvsMedia_Negatively_Enriched_GO_BP.png
	LPSvsMedia_Negatively_Enriched_GO_CC.png
	LPSvsMedia_Negatively_Enriched_GO_MF.png
	LPSvsMedia_Normalized_Counts.png
	LPSvsMedia_PCA_Plot.png
	LPSvsMedia_Positively_Enriched_GO_BP.png
	LPSvsMedia_Positively_Enriched_GO_CC.png
	LPSvsMedia_Positively_Enriched_GO_MF.png
	LPSvsMedia_Raw_Counts.png
	LPSvsMedia_Size_Factors.png
	LPSvsMedia_Sparsity_Plot.png
	LPSvsMedia_Variance_Scree_Plot.png
	LPSvsMedia_Volcano_Plot_DEGs.png
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

The input files require further manual file manipulation and handling, after generation through code in R, thus, final GSEA input files (expression_GSEA_Media_vs_LPS.txt.gct; expression_GSEA_LPS_TNFalpha.txt.gct; phenotype_labels_GSEA_file_Media_vs_LPS.cls; phenotype_labels_GSEA_file_LPS_TNFalpha.cls) are provided for review and perusal. 

Non Pre-Ranked Gene Set Enrichment Analysis (Linux v4.3.3) 

To run GSEA, use the following bash command in the terminal:

bash

(base) ioana@ioana:~/Downloads/GSEA_Linux_4.3.3$ ./gsea.sh
