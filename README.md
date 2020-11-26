# ITH-code
Code written to perform all analysis involved in my Master's thesis about intratumor heterogeneity (ITH). This analysis is the first bioinformatics project that I have authored and serves as my introduction to R and statistical modeling. 

I explore the association between genetic ITH of and the gene expression of primary tumors by using a linear model and the available data from The Cancer Genome Atlas (TCGA). In other words, I look for genes whos expression might be associated with the level of ITH and could thereby serve as a potential biomarker of ITH. Two different computational methods, EXPANDS and PhyloWGS, are used for ITH measures and the whole analysis is performed with both separately and the results compared at each step. 

The analysis involves three essential steps: 
1) measure ITH of tumor samples, separately for each cancer type, by using a computational method called EXAPNDS. (PhyloWGS ITH measures are brought from literature)
2) reduce the search space of the transcriptome to a few differentially expressed genes by performing Differential Expression Analysis between samples of low- and high ITH
3) apply a linear model with lasso regularization on each sample, using the respective ITH measure as resopnse and the selected genes' expression data as predictor variables. 

The code is split into several .R scripts and are used as follows:

1)
ExpandsRun.R - Using EXPANDS to measure ITH in all available samples with sufficient data. The measure of ITH used in this analysis is the number of genetically distinct cellular subpopulations (SPs) withing the tumor. Expands uses Single Nucleotide Variants derived from WES data and Copy Number Variants derived from microarrays through the TCGA computational pipelines. 

2)
DEA_input.R - Split samples into groups of low- and high ITH based on EXPANDS and PhyloWGS separately and for each cancer type separately

DEA_run.R - Run Diffefential Expression Analysis between groups of low- and high ITH

DE_filter.R - Select differentially expressed genes for modeling

3) 
LMFIT_allGenes_bycanc.R - linear model selection with lasso regularization, using differentially expressed genes and predictor variables and the sample ITH measure as response

