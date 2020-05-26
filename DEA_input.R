library(readxl)
library(ggplot2)
library(data.table)
library(gridExtra)
library(mclust)

# -----ITH results, output of either ITH method as SPs----------------------------------------

## PhyloWGS SPs from (Raynaud et al., 2018), download link: https://doi.org/10.1371/journal.pgen.1007669.s008
Raynaud2018.PhyloWGSOutput <- as.data.frame(read_excel("journal.pgen.1007669.s008.xlsx"))
phyloSPs <- Raynaud2018.PhyloWGSOutput[, c(2,1,7)]   # subset of data frame
names(phyloSPs) <- c("cancer", "patient", "phylowgs_SPs")
phyloSPs$patient <- substr(phyloSPs$patient, 1,12)
## expands SPs from: 
#source("ExpandsRun.R") # double check & polish
#expandSPs <- as.data.frame(fread("EXPANDS_results_nozero.tab"))
#colnames(expandSPs) <- c("cancer", "patient", "expands_SPs")

#new expands run, due to earlier omission 
expands.newSPs.all <- as.data.frame(fread("TCGA_expands_14_05_2020.tab"))
colnames(expands.newSPs.all) <- c("cancer", "patient", "expands_SPs")
expands.newSPs <- expands.newSPs.all[!expands.newSPs.all$expands_SPs == 0, ]
#involves zeroes for samples not captured




# split by="colname" function 
table.to.list <- function(table, by){
  ## first FUNCTION:  
  # input: table: table to subset into list; by: specification of the column, by values of which to subset, e.g. "cancer"
  # function: make list, rows subset by unique values in one column of the table (e.g. table$cancer)
  # output: list of subsets
  seplist <- list()
  for (i in unique(table[[by]])){
    sepd <- table[table[[by]] == i,]
    seplist[[i]] <- sepd
  }
  return(seplist)
}
# classifier function (quartile borders excluded)
classifier <- function(table, by){
  ##classifier FUNCTION: create column with row lables: "low", "moderate", "high". 
  #lable by column quartiles as given by summary()
  #input: table
  #output: column of lables
  lables <- c()   
  for (i in table[[by]]){                                                                        # write as function, and apply on phyilist
    if (i < summary(table[[by]])["1st Qu."] | 
        i == summary(table[[by]])["Min."]){   # "low" if <= 1st Quartile
      lables <- append(lables, "low")
    } else if (i > summary(table[[by]])["3rd Qu."]){                                             # "high" if > 3rd Quartile
      lables <- append(lables, "high")
    } else {                                                                                                    # "moderate" else
      lables <- append(lables, "moderate")
    }
  }
  return(lables)
}

#expands alone .all includes zeroes
expands.newSPs.bycanc.all <- table.to.list(table = expands.newSPs.all,by = "cancer")
expands.newSPs.bycanc <- table.to.list(table = expands.newSPs,by = "cancer")

#record how many samples per cancer were captured
cancers <- intersect(unique(expands.newSPs$cancer),unique(phyloSPs$cancer))
samples.run <- numeric(length = length(cancers))
names(samples.run) <- cancers
samples.success <- numeric(length = length(cancers))
names(samples.success) <- cancers
for (i in cancers){
  samples.run[[i]] <- nrow(expands.newSPs.bycanc.all[[i]])
  samples.success[[i]] <- nrow(expands.newSPs.bycanc[[i]])
}

#order by samples captured by expands
cancers <- cancers[order(-samples.success)]
#table of samples run successfully
run.table <- data.frame(row.names = cancers, 
                        "nr samples" = as.numeric(samples.run[cancers]),
                        "nr success" = as.numeric(samples.success[cancers]))
#add labels
for (i in cancers){
  tab <- expands.newSPs.bycanc[[i]]
  tab$expands_ITHclass <- classifier(table = tab, by="expands_SPs")
  expands.newSPs.bycanc[[i]] <- tab
}
#phylo alone
phyloSPs.bycanc <- table.to.list(table = phyloSPs,by = "cancer")
#add labels
for (i in names(phyloSPs.bycanc)){
  tab <- phyloSPs.bycanc[[i]]
  tab$phylowgs_ITHclass <- classifier(table = tab, by="phylowgs_SPs")
  phyloSPs.bycanc[[i]] <- tab
}

## MERGE: for method comparison, use only samples for which both methods have given results
MergeSPs <- merge(phyloSPs, expands.newSPs, by = "patient")
colnames(MergeSPs)[4] <- c("cancer"); MergeSPs <- MergeSPs[-2]

overall.cor <- cor.test(MergeSPs$phylowgs_SPs, MergeSPs$expands_SPs)

## split by cancer
MergeSPs.bycanc <- table.to.list(table = MergeSPs, by = "cancer")
## add labels with classifier
for (i in names(MergeSPs.bycanc)){
  tab <- MergeSPs.bycanc[[i]]
  tab$expands_ITHclass <- classifier(table=tab,by="expands_SPs")
  MergeSPs.bycanc[[i]] <- tab
}
for (i in names(MergeSPs.bycanc)){
  tab <- MergeSPs.bycanc[[i]]
  tab$phylowgs_ITHclass <- classifier(table=tab,by="phylowgs_SPs")
  MergeSPs.bycanc[[i]] <- tab
}
