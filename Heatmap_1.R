library(gplots) # bluered
library(pheatmap)

# build heatmap matrix:
# Rows: top8 genes; Cols: Cancers

## heatmap matrix [genes, cancers]
genes <- top6.expan             #genes used for model
##coefs from model fit
coef_df <- data.frame(matrix(nrow = length(genes), 
                             ncol = length(names(Coefs.expan.top6))),
                             row.names = genes); 
colnames(coef_df) <- names(Coefs.expan.top6)
## fill in coefs
for (i in colnames(coef_df)){ 
coefs <- Coefs.expan.top6[[i]][-1,]       #take coefs+gene names from results, -intercept
coef_df[as.character(coefs[[1]]), i] <- as.numeric(coefs[[2]])
}
#remove NAs
coef_df <- coef_df[rowSums(is.na(coef_df)) != ncol(coef_df),]
coef_df <- coef_df[,colSums(is.na(coef_df)) != nrow(coef_df)]

extreme_point <- max(abs(coef_df), na.rm = T)
png(paste(getwd(),"/Plots/Heatmaps/expands_top6_new_heatmap.png", sep=""))
eheat.top6.new <- pheatmap(
  coef_df,
  cluster_rows = F,
  cluster_cols = F,
  color = bluered(9),
  breaks = c(seq(-extreme_point, 0, length = 5), seq(0.000001, extreme_point, length = 5)),
  na_col = 0,
  angle_col = 315,
  fontsize = 8,
#  annotation_row = as.data.frame(Adjusted_R_Squared),
#  display_numbers = switch(order + 1, F, vorder),
  number_color = 1)
print(eheat.top6.new)
dev.off()

#with separate genes

all.genes <- Reduce(union, list(row.names(top100.expan$LUAD), 
                   row.names(top100.expan$LUSC),
                   row.names(top100.expan$STAD),
                   row.names(top100.expan$BRCA),
                   row.names(top100.expan$HNSC),
                   row.names(top100.expan$COAD),
                   row.names(top100.expan$BLCA),
                   row.names(top100.expan$PRAD),
                   row.names(top100.expan$UCEC),
                   row.names(top100.expan$ESCA),
                   row.names(top100.expan$LIHC),
                   row.names(top100.expan$CESC)))





