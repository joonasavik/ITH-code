library(ISLR)
library(glmnet)
library(pheatmap)
library(edgeR)
library(ggplot2)


Gene_names <- fread("mart_uniq.txt");   
Gene_names <- data.frame(Gene_name = Gene_names$`Gene name`, row.names = Gene_names$`Gene stable ID`)

# normalize expr data for modelling (CPM)
CPMs.expan <- list()
SPs.expan <- list()
CPMs.phylo <- list()
SPs.phylo <- list()

CPMs.merge <- list()
SPs.merge <- list()
cancers <- union(names(exp_data.expan.merge), names(exp_data.phylo.merge))
for(i in "STAD"){  
  SPs <- MergeSPs.bycanc[[i]][,c("patient","phylowgs_SPs","expands_SPs")] 
  exprRaw <- read.table(paste("Pancancer_GDC/Transcriptome/TCGA-",i,"_htseq_counts.tab",
                              sep=""),
                        sep="\t", header=T, quote = "", comment.char = '', stringsAsFactors = F)     # RNA-seq data, pats as colnames
  expr <- exprRaw[-c(1:5),]                                                                       # remove first 5 rows of quality data
  rownames(expr) <- expr[,1]      # set genes, 1st col, as rownames
  expr <- expr[-c(1)]             # remove 1st col
  # Gene name conversion
  row.names(expr) <- substr(rownames(expr), 1, 15)
  expr <- merge(Gene_names, expr, by = 0)
  expr <- expr[!duplicated(expr[,2]), ]
  row.names(expr) <- expr[,2]; expr <- expr[-c(1,2)]
  # normal tissue samples removed
  expr <- expr[,colnames(expr)[grep(".01A", colnames(expr))]] 
  # select expr samples matching SP samples (ITH)
  SPs[[1]] <- gsub("-", ".", SPs[[1]])          #matching sample names in ITH data
  colnames(expr) <- substr(colnames(expr), 1, 12)     #matching sample names in RNA data
  expr <- expr[,colnames(expr)%in%SPs[[1]]]        #selecting 
  # normalize TMM, CPM with edgeR
  y <- DGEList(counts=expr)    # not filtering low expression, was done in DEA      
  keep <- filterByExpr(y)                         # Genes that have sufficiently large counts to be retained in a statistical analysis
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)      # scale the raw library sizes
  expr <- cpm(y, normalized.lib.sizes = TRUE)  # Counts per million reads mapped
  print(i)
  
  CPMs.merge[[i]] <- expr
  SPs.merge[[i]] <- SPs
  rm(expr)
  rm(y)
}

modMats.expan.top6 <- list()
testMats.expan.top6 <- list()
modMats.phylo.top6 <- list()
testMats.phylo.top6 <- list()

modMats.phylo.top100 <- list()
testMats.phylo.top100 <- list()
modMats.expan.top100 <- list()
testMats.expan.top100 <- list()

cancs <- names(top100.expan.merge)  # which cnacers?        
for (i in cancs){  
  genes <- top100.expan.merge[[i]]
  SPs <- SPs.merge[[i]][,c("patient","expands_SPs")]      #!!!
  expr <- CPMs.merge[[i]][row.names(CPMs.merge[[i]])%in%genes,]    #!!!
  expr <- as.data.frame(t(expr))
  expr$patient <- row.names(expr)
  #merge SPs and CPMs
  modMat <- merge(SPs, expr, by = "patient")
  row.names(modMat) <- modMat[[1]]
  modMat <- modMat[ ,-1]
  modMat[ ,-1] <- scale(modMat[ ,-1])
  colnames(modMat)[1] <- "SPs"
  #test set
  if (nrow(modMat) > 150){             # test/training split for larger sample sizes
    indx <- sample(rownames(modMat), 0.25*nrow(modMat))
    testSet <- modMat[rownames(modMat) %in% indx, ]
    testMats.expan.top100[[i]] <- testSet         #!!!
    modMat <- modMat[!(rownames(modMat) %in% indx), ]
  }
  modMats.expan.top100[[i]] <- modMat             #!!!
  print(i)
}


## FUNCTION: MODEL FIT
#output: 
# model (list) the model
# Coefs (list) vector of coefficients of the model


model.expan.top6 <- list()
Coefs.expan.top6 <- list()
model.phylo.top6 <- list()
Coefs.phylo.top6 <- list()

model.expan.top100 <- list()
Coefs.expan.top100 <- list()
model.phylo.top100 <- list()
Coefs.phylo.top100 <- list()
#input: modMats, the model matrix
for (i in names(modMats.expan.top100)){
  data <- modMats.expan.top100[[i]]
  x <- model.matrix(SPs~.-1, data=data)  
  y <- data$SPs
  # glmnet does a shrinkage on the variables and performs subset selection by shrinking the coefficients towards zero
  cv.lasso <- cv.glmnet(x,y,nfolds = 10)               # cv.glmnet does CROSS VALIDATION to pick an appropriate lambda (shrinkage factor)
  model.expan.top100[[i]] <- glmnet(x,y, lambda = cv.lasso$lambda.1se)    # lasso summary
  tmp_coefs <- coef(cv.lasso)
  coefs <- data.frame(name = tmp_coefs@Dimnames[[1]][tmp_coefs@i + 1], coefficient = tmp_coefs@x)   #coefs extracted
#  if (nrow(coefs) > 1){
    Coefs.expan.top100[[i]] <- coefs      # coeficients for genes from shirkage, from this build COEF MATRIX
#  }
}

#intersecting coefs
intcoefs <- list()
for (i in intersect(names(Coefs.expan.top100), names(Coefs.phylo.top100))){
  intcoefs[[i]] <- intersect(as.character(Coefs.expan.top100[[i]][["name"]][-1]), as.character(Coefs.phylo.top100[[i]][["name"]][-1]))
}

## TEST RUN models
rmse.expan.top6 <- c()
pred.expan.top6 <- list()
rmse.phylo.top6 <- c()
pred.phylo.top6 <- list()

rmse.expan.top100 <- c()
pred.expan.top100 <- list()
rmse.phylo.top100 <- c()
pred.phylo.top100 <- list()


#calc RMSE (rmse...) and a list of predictions per cancer (pred...)
for (i in names(testMats.expan.top100)){
  data <- testMats.expan.top100[[i]]  #!!!
  x.test <- model.matrix(SPs~.-1, data=data)
  y.test <- data$SPs
  pred <- predict(model.expan.top100[[i]], x.test)  #!!! model[[i]] is the model for cancer i, used to predict SPs from test data
  rmse <- sqrt(apply((y.test-pred)^2,2,mean))
  
  rmse.expan.top100 <- append(rmse.expan.top100, rmse) #!!!
  compare <- data.frame(row.names = row.names(data), expands_SPs = data$SPs, preds = pred)
  colnames(compare)[2] <- "preds"
  pred.expan.top100[[i]] <- compare
}
names(rmse.expan.top100) <- names(testMats.expan.top100)

## cor test
cor.phylo.top100 <- list()
cor.expan.top100 <- list()
for (i in names(pred.expan.top100)){
  cor.expan.top100[[i]] <- cor.test(x = pred.expan.top100[[i]][["expands_SPs"]], y = pred.expan.top100[[i]][["preds"]])
}

# compile corr estimates and p values

final.cancs <- intersect(names(cor.expan.top100), names(cor.phylo.top100))
cors.p <- numeric(length = 4); names(cors.p) <- final.cancs
cors.e <- numeric(length = 4); names(cors.p) <- final.cancs
pval.p <- numeric(length = 4); names(cors.p) <- final.cancs
pval.e <- numeric(length = 4); names(cors.p) <- final.cancs
for (i in final.cancs){
  cors.p[[i]] <- cor.phylo.top100[[i]][["estimate"]]
  cors.e[[i]] <- cor.expan.top100[[i]][["estimate"]]
  pval.p[[i]] <- cor.phylo.top100[[i]][["p.value"]]
  pval.e[[i]] <- cor.expan.top100[[i]][["p.value"]]
}
## TABLE compare TEST ERRORS and CORRELATION
TEST_table <- data.frame(row.names = final.cancs, 
                         p_RMSE=signif(rmse.phylo.top100[final.cancs],3), 
                         e_RMSE=signif(rmse.expan.top100[final.cancs],3),
                         p_Pearsoncors=signif(cors.p,3),
                         e_Pearsoncors=signif(cors.e,3),
                         p_pval=signif(pval.p,3),
                         e_pval=signif(pval.e, 3))
write.table(TEST_table, "COMPARISON_table.txt", sep="\t", quote = FALSE)



## RUN DOTPLOTS for ease of interpretations
pdf(paste(getwd(),"/Plots/DOTPLOT_PHYLO_TOP100.pdf", sep=""))    #change name of file accordingly
for(i in names(testMats.phylo.top100)){   #names(pred.phylo.top6)
  preds <- pred.phylo.top100[[i]] 
  stats <- cor.phylo.top100[[i]]
  #dotplot
  par(fig=c(0.05,0.95,0,0.9))
  plot(preds$phylowgs_SPs, preds$preds,
       xlab = "PhyloWGS SPs",    
       ylab = "model predictions",
       type = "p")
  #boxplot
  par(fig=c(0.05,0.95,0.62,0.98),new=TRUE)
  boxplot(preds$phylowgs_SPs, horizontal=TRUE, axes=FALSE)
  title(paste(i,": PhyloWGS ITH predicted, top DE genes"))                      #change title accordingly
  mtext(paste("Pearson's estimate: ",round(stats[["estimate"]],3),", p-value: ", signif(stats[["p.value"]],3)))
}
dev.off() 



## DOTPLOT, predictions vs SPs, COMPARISON
png("LUAD_COMP.png")
for(i in final.cancs[4]){   #intersect(names(pred.expan.top100),names(pred.phylo.top100))
  dotphylo <- pred.phylo.top100[[i]]
  dotexpan <- pred.expan.top100[[i]] 
  # this is how I add the title to the "doubleplot"
  par(fig=c(0,1,0,1))
  plot(1,1,type="n",axes=FALSE, ann=FALSE)
  title(paste(i),
        adj  = 0.5,
        line = -11)
  # boxplots of SPs input
  par(fig=c(0,0.5,0.31,0.68),new=TRUE)
  boxplot(dotphylo$phylowgs_SPs, horizontal=TRUE, axes=FALSE)
  par(fig=c(0.5,1,0.31,0.68), new=TRUE)
  boxplot(dotexpan$expands_SPs, horizontal=TRUE, axes=FALSE)
  # dotplots of SPs vs predictions
  par(fig=c(0,0.5,0,0.6),new=TRUE)
  plot(dotphylo$phylowgs_SPs, dotphylo$preds,
       xlab = "PhyloWGS SPs", 
       ylab = "model predictions", 
       type = "p")
#  title(sub = paste("corr P:", signif(cors.p[[i]],3), 
#                    ", pval:", signif(pval.p[[i]],3), 
#                    ", RMSE:", signif(rmse.phylo.top100[[i]],3)))
#
  par(fig=c(0.5,1,0,0.6),new=TRUE)
  plot(dotexpan$expands_SPs, dotexpan$preds,
       xlab = "expands SPs",    
       ylab = "model predictions",
       type = "p")
#  title(sub = paste("corr P:",signif(cors.e[[i]],3), 
#                    ", pval:", signif(pval.e[[i]],3), 
#                    ", RMSE:", signif(rmse.expan.top100[[i]],3)))
}
dev.off()






