library(mclust)
library(ggplot2)

#---------------------------------------plots and extra stuff----------------------------------------------

## -----BOXPLOT of SPs [COULD MAKE FUNCTION]---------------------------------------
# sorter: vector of cancers sorted by cohort size
sorter <- numeric(length = length(unique(MergeSPs$cancer))); names(sorter) <- unique(MergeSPs$cancer)
for(i in unique(MergeSPs$cancer)){
  sorter[[i]] <- nrow(MergeSPs[MergeSPs$cancer%in%i,])
}
sorter <- sort(sorter)
merg.samples <- data.frame()
# merg.samples: cancers sorted by cohort size
for(i in names(sorter)){
  merg.samples <- rbind(merg.samples, MergeSPs[MergeSPs$cancer%in%i,])
} 

png("BOXplot_new_ITHmethods.png")
pbox <- ggplot(merg.samples, aes(y=cancer, x=phylowgs_SPs)) +   #how to generalize name of third column for y (i.e. phylowgs_SPs or expands_SPs)?
  geom_boxplot() +
  scale_x_reverse(limits = c(15,1),name = "PhyloWGS SPs") +
  scale_y_discrete(limits = names(sorter), name=NULL) +
  labs(title = "ITH measured with PhyloWGS") 
#  print(pbox)
ebox <- ggplot(merg.samples, aes(y=cancer, x=expands_SPs)) +   #how to generalize name of third column for y (i.e. phylowgs_SPs or expands_SPs)?
  geom_boxplot() +
  scale_x_continuous(limits = c(1,15), name = "expands SPs") +
  scale_y_discrete(limits = names(sorter), position = "right",name=NULL) +
  labs(title = "ITH measured with expands")
#  print(ebox)
pebox <- grid.arrange(pbox, ebox, nrow=1) 
#  labs(title="ITH per cancer type as SP ranges")
print(pebox)
dev.off()

#save alternative...
ggsave(
  "Boxplot_ITH.pdf",
  plot = pebox,
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 16,
  height = 13,
  units = c("cm"),
  dpi = 300,
  limitsize = TRUE
)

## -----BARPLOT of GROUPS [THIS NEEDS TO BA A FUNCTION......]---------------------------------
H <- rep("high", 3)
M <- rep("moderate", 3)
L <- rep("low", 3)
#rev(names(merged.SPs))
png("BARplot_STAD_groups_new.png")
for (i in "STAD"){        #rev(names(sorter))
  df <- MergeSPs.bycanc[[i]]
  #Phylo Bars
  if (length(summary(as.factor(df$expands_ITHclass))) == 3){ #if "low", "moderate", "high" groups exist
    PHigh <- df[df$phylowgs_ITHclass == "high",]
    PHigh_ex <- numeric(length=3);names(PHigh_ex) <- c("low","moderate","high")
    PHigh_ex["low"] <- nrow(PHigh[PHigh$expands_ITHclass == "low",])
    PHigh_ex["moderate"] <- nrow(PHigh[PHigh$expands_ITHclass == "moderate",])
    PHigh_ex["high"] <- nrow(PHigh[PHigh$expands_ITHclass == "high",])
    
    PMod <- df[df$phylowgs_ITHclass == "moderate",]
    PMod_ex <- numeric(length=3);names(PMod_ex) <- c("low","moderate","high")
    PMod_ex["low"] <- nrow(PMod[PMod$expands_ITHclass == "low",])
    PMod_ex["moderate"] <- nrow(PMod[PMod$expands_ITHclass == "moderate",])
    PMod_ex["high"] <- nrow(PMod[PMod$expands_ITHclass == "high",])
    
    PLow <- df[df$phylowgs_ITHclass == "low",]
    PLow_ex <- numeric(length=3);names(PLow_ex) <- c("low","moderate","high")
    PLow_ex["low"] <- nrow(PLow[PLow$expands_ITHclass == "low",])
    PLow_ex["moderate"] <- nrow(PLow[PLow$expands_ITHclass == "moderate",])
    PLow_ex["high"] <- nrow(PLow[PLow$expands_ITHclass == "high",])
    
    PplotH <- data.frame(P_class=H,
                         E_class=names(PHigh_ex),
                         cases=PHigh_ex,
                         row.names = NULL)
    PplotM <- data.frame(P_class=M,
                         E_class=names(PMod_ex),
                         cases=PMod_ex,
                         row.names = NULL)
    PplotL <- data.frame(P_class=L,
                         E_class=names(PLow_ex),
                         cases=PLow_ex,
                         row.names = NULL)
    Pplot <- rbind(PplotH,PplotM, PplotL)

    EHigh <- df[df$expands_ITHclass == "high",]
    EHigh_ph <- numeric(length=3);names(EHigh_ph) <- c("low","moderate","high")
    EHigh_ph["low"] <- nrow(EHigh[EHigh$phylowgs_ITHclass == "low",])
    EHigh_ph["moderate"] <- nrow(EHigh[EHigh$phylowgs_ITHclass == "moderate",])
    EHigh_ph["high"] <- nrow(EHigh[EHigh$phylowgs_ITHclass == "high",])
    
    EMod <- df[df$expands_ITHclass == "moderate",]
    EMod_ph <- numeric(length=3);names(EMod_ph) <- c("low","moderate","high")
    EMod_ph["low"] <- nrow(EMod[EMod$phylowgs_ITHclass == "low",])
    EMod_ph["moderate"] <- nrow(EMod[EMod$phylowgs_ITHclass == "moderate",])
    EMod_ph["high"] <- nrow(EMod[EMod$phylowgs_ITHclass == "high",])
    
    ELow <- df[df$expands_ITHclass == "low",]
    ELow_ph <- numeric(length=3);names(ELow_ph) <- c("low","moderate","high")
    ELow_ph["low"] <- nrow(ELow[ELow$phylowgs_ITHclass == "low",])
    ELow_ph["moderate"] <- nrow(ELow[ELow$phylowgs_ITHclass == "moderate",])
    ELow_ph["high"] <- nrow(ELow[ELow$phylowgs_ITHclass == "high",])
    
    EplotH <- data.frame(E_class=H,
                         P_class=names(EHigh_ph),
                         cases=EHigh_ph,
                         row.names = NULL)
    EplotM <- data.frame(E_class=M,
                         P_class=names(EMod_ph),
                         cases=EMod_ph,
                         row.names = NULL)
    EplotL <- data.frame(E_class=L,
                         P_class=names(ELow_ph),
                         cases=ELow_ph,
                         row.names = NULL)
    if (nrow(EMod) > nrow(PMod)){
      yscale <- nrow(EMod)
    } else {
      yscale <- nrow(PMod)
    }
    Eplot <- rbind(EplotH,EplotM, EplotL)
    ebar <- ggplot(data=Eplot, aes(x=E_class, y=cases, fill=P_class)) +
      geom_bar(stat="identity") +
      ggtitle(paste(i,": expands SPs", sep="")) +
      xlab("expands groups") + ylab("Number of Cases") +
      ylim(c(0,yscale)) +
      labs(fill = "PhyloWGS groups")
    pbar <- ggplot(data=Pplot, aes(x=P_class, y=cases, fill=E_class)) +
      geom_bar(stat="identity") +
      ggtitle(paste(i,": PhyloWGS SPs", sep="")) +
      xlab("PhyloWGS groups") + ylab("Number of Cases") +
      ylim(c(0,yscale)) +
      labs(fill = "expands groups")
    
    pebar <- grid.arrange(pbar, ebar, nrow=1) 
    print(pebar)
  }
}
dev.off()


## -----ADJUSTED RAND INDEX--------------------------------------------------------------------
cancs <- c()
inds <- c()
for (i in rev(names(sorter))){
  p.labels <- MergeSPs.bycanc[[i]][["phylowgs_ITHclass"]]
  e.labels <- MergeSPs.bycanc[[i]][["expands_ITHclass"]]
  ind <- adjustedRandIndex(p.labels, e.labels)
  cancs <- append(cancs, i)
  inds <- append(inds, ind)
}
mean(inds)
RAND_table <- data.frame(cancer=cancs,Adj.Rand.Index=signif(inds,2))
write.table(RAND_table, file = "RAND_table.txt",sep = "\t", quote = FALSE,row.names = FALSE)
lowRand <- RAND_table[RAND_table[[2]]<0.1,]

## -----NON MATCHING SAMPLES (expands & Phylowgs)-----------------------------------------------
## split by cancer 
eSPs.bycanc <- table.to.list(table = expandSPs, colname = "cancer")
pSPs.bycanc <- table.to.list(table = phyloSPs, colname = "cancer")
## adding the classifier lables to the tables, write FUNCTION
for (i in names(eSPs.bycanc)){
  tab <- eSPs.bycanc[[i]]
  tab$expands_ITHclass <- classifier(table=tab,by="expands_SPs")
  eSPs.bycanc[[i]] <- tab
}
for (i in names(pSPs.bycanc)){
  tab <- pSPs.bycanc[[i]]
  tab$phylowgs_ITHclass <- classifier(table=tab,by="phylowgs_SPs")
  pSPs.bycanc[[i]] <- tab
}
common.cancs <- intersect(names(eSPs.bycanc),names(pSPs.bycanc))
eSPs.bycanc <- eSPs.bycanc[names(eSPs.bycanc)%in%common.cancs]
pSPs.bycanc <- pSPs.bycanc[names(pSPs.bycanc)%in%common.cancs]
## NON-matches lists
eSPs <- data.frame()
pSPs <- data.frame()
for (i in names(eSPs.bycanc)){
  eSPs <- rbind(eSPs, eSPs.bycanc[[i]])
  pSPs <- rbind(pSPs, pSPs.bycanc[[i]])
}
## compare NON-match samples
no.match <- as.data.frame(merge(eSPs, pSPs, by="patient", all=TRUE))
no.match <- no.match[is.na(no.match$expands_SPs)|is.na(no.match$phylowgs_SPs),]
e.unmatched <- summary(as.factor(no.match["expands_ITHclass"][!is.na(no.match["expands_ITHclass"])]))
p.unmatched <- summary(as.factor(no.match["phylowgs_ITHclass"][!is.na(no.match["phylowgs_ITHclass"])]))



## DOTPLOT of predictions vs SPs, both methods (merged)
png("DOTplot_pred_LUAD.png")
for(i in c("LUAD")){
  dotphylo <- pred.phylo[[i]]
  dotexpan <- pred.expan[[i]] 
  # this is how I add the title
  par(fig=c(0,1,0,1))
  plot(1,1,type="n",axes=FALSE, ann=FALSE)
  title(paste(i),
        adj  = 0.5,
        line = -11)
  # boxplots of SPs input
  par(fig=c(0,0.5,0.31,0.68),new=TRUE)
  boxplot(dotphylo$phylo_SPs, horizontal=TRUE, axes=FALSE)
  par(fig=c(0.5,1,0.31,0.68), new=TRUE)
  boxplot(dotexpan$expan_SPs, horizontal=TRUE, axes=FALSE) 
  # dotplots of SPs vs predictions
  par(fig=c(0,0.5,0,0.6),new=TRUE)
  plot(dotphylo$phylo_SPs, dotphylo$preds,
       xlab = "PhyloWGS SPs", 
       ylab = "model predictions", 
       type = "p")
  par(fig=c(0.5,1,0,0.6),new=TRUE)
  plot(dotexpan$expan_SPs, dotexpan$preds,
       xlab = "expands SPs",    
       ylab = "model predictions",
       type = "p")
}
dev.off()

## ggplot attempt Dotplot
pdf(file = "DOTplot_predictions")
for(i in names(pred.phylo)){
  dotphylo <- pred.phylo[[i]]
  dotexpan <- pred.expan[[i]]
  dp.phy <- ggdotplot(dotphylo, x = "phylo_SPs", y = "preds",xlab = "PhyloWGS SPs", ylab = "model predictions") #ADD RMSE somehow
  dp.exp <- ggdotplot(dotexpan, x = "expan_SPs", y = "preds",xlab = "expands SPs", ylab = "model predictions")       
  double <- ggarrange(dp.exp, dp.phy + title(main=i),
                      labels = c("expands", "PhyloWGS"),
                      ncol = 2)
  print(double)
}
dev.off()
