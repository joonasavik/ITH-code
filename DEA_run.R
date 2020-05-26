library(edgeR)
library(dplyr)
library(data.table)
library(limma)

## ITH table:
#source("DEA_input.R")
set.seed(666)

## DEA loop
Gene_names <- fread("mart_uniq.txt");   
Gene_names <- data.frame(Gene_name = Gene_names$`Gene name`, row.names = Gene_names$`Gene stable ID`)
exp_data.expan <- list()
exp_data.phylo <- list()

exp_data.phylo.merge <- list()
exp_data.expan.merge <- list()
colname <- "expands_ITHclass"
for (i in names(MergeSPs.bycanc)){   # samples from expands-phylo merger
  #ITH table input
  ITH_set <- MergeSPs.bycanc[[i]] 
  #RNA-seq data, TCGA_counts: samples as cols, genes as rows
  expr <- read.table(paste("Pancancer_GDC/Transcriptome/TCGA-",i,"_htseq_counts.tab",
                           sep=""),
                     sep="\t", header=T, quote = "", comment.char = '', stringsAsFactors = F)     
  #format TCGA_counts input file
  expr <- expr[-c(1:5),]          # remove first 5 rows of quality data
  rownames(expr) <- expr[,1]      # set genes, 1st col, as rownames
  expr <- expr[-c(1)]             # remove 1st col
  ## Gene name conversion
  row.names(expr) <- substr(rownames(expr), 1, 15)
  expr <- merge(Gene_names, expr, by = 0)
  expr <- expr[!duplicated(expr[,2]), ]
  row.names(expr) <- expr[,2]; expr <- expr[-c(1,2)]

  # RUN DEA if >5 samples per group
  if (length(ITH_set[[colname]][ITH_set[[colname]] == "low"])<=5 |  
      length(ITH_set[[colname]][ITH_set[[colname]] == "high"])<=5){
    next 
  }
  ITH_set[["patient"]] <- gsub("-", ".", ITH_set$patient)
  ITH_low <- ITH_set[ITH_set[colname] == "low",]
  ITH_high<- ITH_set[ITH_set[colname] == "high",]
  ## RNAseq samples
  RNA_set <- data.frame("samples" =colnames(expr)[grep(".01A", colnames(expr))],
                        "patient" = substr(colnames(expr)[grep(".01A", colnames(expr))], 1, 12),
                        stringsAsFactors = F)  # only tumor samples "01A"
  # RUN DEA if count data exists for >5 patients
  if (nrow(RNA_set[RNA_set$patient %in% ITH_low$patient,])<=5 |  
      nrow(RNA_set[RNA_set$patient %in% ITH_high$patient,])<=5){
    next
  }
  RNA_low <- RNA_set[RNA_set$patient %in% ITH_low$patient ,]    # RNA samples, low ITH
  RNA_low$ITHclass <- "low"  
  RNA_high<- RNA_set[RNA_set$patient %in% ITH_high$patient ,]   # RNA samples, high ITH
  RNA_high$ITHclass <- "high"
  expr_matched <- as.data.frame(expr[,c(RNA_low$samples,  # here, the above, stringsAsFactors = F, is crucial
                                        RNA_high$samples)])             # read counts: low ITH vs high ITH
  
  colData <- data.frame(rbind(RNA_low, RNA_high), stringsAsFactors = F)
  rownames(colData) <- colData$samples
  colData$samples <- NULL
  colData$ITHclass <- factor(colData$ITHclass, levels = c("low","high"))

  ## edgeR: normalizing with Counts Per Million (CPM)
  y <- DGEList(counts=expr_matched)
  #Maybe choose another critera for filtering low TCGA_countsessed genes... since we have so many samples
  keep <- filterByExpr(y)                         # Genes that have sufficiently large counts to be retained in a statistical analysis
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)                         # Calculate normalization factors to scale the raw library sizes.
  cpm <- cpm(y, normalized.lib.sizes = TRUE)
  # With so many samples, we can use other methods to make the statistical test...
  # Voom, or a simple wilcoxon!
  
  ## limma    
  design <- model.matrix(~ ITHclass, data=colData)    
  y <- voom(y, design, plot = T)
  fit <- lmFit(y, design)                         # fit where variance of low exp genes is shrunk towards fit with increasing penalty
  tmp <- eBayes(fit)                              # DEA
  top.table <- topTable(tmp, sort.by = "P", n = Inf, coef="ITHclasshigh")   # see limma guide 
  
  exp_data.expan.merge[[i]][["cpm"]] <- cpm
  exp_data.expan.merge[[i]][["results"]] <- top.table
  exp_data.expan.merge[[i]][["metadata"]] <- colData

  print(i)
}

##checking DEA direction:
lowsamples <- exp_data.eheat$LUSC$metadata[exp_data.eheat$LUSC$metadata[["ITHclass"]]=="low",]
lowsamples <- row.names(lowsamples)
highsamples <- exp_data.eheat$LUSC$metadata[exp_data.eheat$LUSC$metadata[["ITHclass"]]=="high",]
highsamples <- row.names(highsamples)
WBP1L_data <- exp_data.eheat$LUSC$cpm["WBP1L",]
WBP1L_low <- WBP1L_data[names(WBP1L_data)%in%lowsamples]
WBP1L_high <- WBP1L_data[names(WBP1L_data)%in%highsamples]

  
  


  
  
  
  
  
  
