library(dbplyr)
library(dplyr)
library(data.table)

# filter significant MAKE FUNC
expan.SDEGs <- list()
phylo.SDEGs <- list()

phylo.merge.SDEGs <- list()
expan.merge.SDEGs <- list()

expan.SDEGs.merge <- list()
for (i in names(exp_data.expan.merge)){
  approach <- exp_data.expan.merge[[i]][["results"]]
  if (nrow(approach[approach["adj.P.Val"]<0.05,])>=1){     # filter significant DE genes
    expan.SDEGs.merge[[i]] <- approach[approach["adj.P.Val"]<0.05,]
  }
}
#get intersections 
SDEG.cancers <- intersect(names(expan.SDEGs.merge), names(phylo.SDEGs.merge))

eSDEGs <- c()
pSDEGs <- c()
intersecting <- c()
samples <- c()
elow <- c()
ehigh <- c()
plow <- c()
phigh <- c()
int.top100 <- c()
for (i in SDEG.cancers){
  eSDEGs <- append(eSDEGs, nrow(expan.SDEGs.merge[[i]]))
  pSDEGs <- append(pSDEGs, nrow(phylo.SDEGs.merge[[i]]))
  intersecting <- append(intersecting, length(intersect(rownames(expan.SDEGs.merge[[i]]), 
                                                        rownames(phylo.SDEGs.merge[[i]]))))
  samples <- append(samples, ncol(CPMs.merge[[i]]))
  egroups <- exp_data.expan.merge[[i]][["metadata"]]["ITHclass"]
  pgroups <- exp_data.phylo.merge[[i]][["metadata"]]["ITHclass"]
  elow <-  append(elow, length(egroups[egroups$ITHclass == "low",]))
  ehigh <- append(ehigh, length(egroups[egroups$ITHclass == "high",]))
  plow <-  append(plow, length(pgroups[egroups$ITHclass == "low",]))
  phigh <- append(phigh, length(pgroups[egroups$ITHclass == "high",]))
  int.top100 <- append(int.top100, length(intersect(top100.expan.merge[[i]], 
                                             top100.phylo.merge[[i]])))
}
SDEG.table <- data.frame(cancers = SDEG.cancers,
                         samples = samples,
                         exp.low = elow,
                         exp.high = ehigh,
                         phy.low = plow,
                         phy.high = phigh,
                         exp.SDEGs = eSDEGs,
                         phy.SDEGs = pSDEGs,
                         inters = intersecting,
                         int.top100 = int.top100)
write.table(SDEG.table, file = "SDEG_table.txt",sep = "\t", quote = FALSE,row.names = FALSE)

top100.expan.merge

#for run.table
for (i in names(expan.SDEGs)){
  nrow(expan.SDEGs[[i]])
}




sorter.expan <- numeric(length = length(names(expan.merge.SDEGs))); 
names(sorter.expan) <- names(expan.merge.SDEGs)
for(i in names(expan.merge.SDEGs)){
  sorter.expan[[i]] <- nrow(expan.merge.SDEGs[[i]])
}
sorter.expan <- rev(sort(sorter.expan))


#length(intersect(row.names(top100.expan[["BRCA"]]), row.names(top100.phylo[["BRCA"]])))

##filter FUNCTION: filter DE genes by absolute logFC and avegare expression 
#input: DEA output table (limma)
filter.DEgenes <- function(DEgenes){
  #order by abs(logFC)
  DEgenes <- DEgenes[order(-abs(DEgenes[["logFC"]])),]
  #select top 500, abs(logFC)
  DEgenes <- DEgenes[1:500,]
  DEgenes <- DEgenes[order(-abs(DEgenes[["AveExpr"]])),]
  #select top 100, aveExpr
  DEgenes <- DEgenes[1:100,]
}
# DE output (list):
#exp_data.phylo.Qex
#exp_data.expan.Qex
# SDEGs (list):
#phylo.Qex.SDEGs
#expan.Qex.SDEGs
top100.phylo.merge <- list()
top100.expan.merge <- list()
for (i in names(sorter.expan)){
  SDEGs <- expan.merge.SDEGs[[i]]
  if (nrow(SDEGs)<500){
    next
  }else{
    top100 <- filter.DEgenes(SDEGs)
    top100.expan.merge[[i]] <- rownames(top100)
  }
}

## top100.phylo and top100.expan have the genes for the model :)
# 25 are in common. Could be interesting to look at them more closely. Also look at GSEA pathways that are in common
length(intersect(row.names(top100.expan[["BRCA"]]), row.names(top100.phylo[["BRCA"]])))

top5.expan <- Reduce(intersect, list(rownames(expan.SDEGs$LUAD),
                                     rownames(expan.SDEGs$STAD),
                                     rownames(expan.SDEGs$HNSC),
                                     rownames(expan.SDEGs$BLCA),
                                     rownames(expan.SDEGs$BRCA),
                                     rownames(expan.SDEGs$UCEC)))
                                   

top6cancers <- c("LUAD","STAD","BRCA","HNSC","UCEC","BLCA")


# filter common genes
top6.phylo <- Reduce(intersect, list(rownames(phylo.SDEGs$THYM), 
                               rownames(phylo.SDEGs$BRCA), 
                               rownames(phylo.SDEGs$PRAD), 
                               rownames(phylo.SDEGs$LUAD), 
                               rownames(phylo.SDEGs$STAD), 
                               rownames(phylo.SDEGs$KIRC)))

top6.phylo.cancers <- c("THYM","BRCA","PRAD","LUAD","STAD","KIRC")
top5.phylo.cancers <- c("THYM","BRCA","PRAD","LUAD","STAD")
