########################
## This program is a RUVg function that provides 
## 1) (optional) non-HK gene normalized correlation plot of ESP1 and GATA3
## 2) PCA plot of designated K value
## Input file/data:
## raw gene data, phenotype data, feature data, flag
########################


pca_corr_kfold <- function(study_name, corr_plot = F, 
                           data, phenodata, featuredata, 
                           category, k_fold){
  
  # Set work direction (where you want to put plots and tables)
  # setwd("C:/Users/Amy/Desktop/Lab/Ghana/Tables")
  
  #### Start of program #####
  
  set<-newSeqExpressionSet(as.matrix(data),
                           phenoData=phenodata,
                           featureData=featuredata)
  
  # Between-lane normalization for sequencing depth and 
  # possibly other distributional differences between lanes
  # Upper quartile normalization (Bullard 2010, BMC Bioinformatics)
  set <- betweenLaneNormalization(set, which="upper")
  
  # Isolate housekeeping genes.
  cIdx <- rownames(set)[fData(set)$Gene_Class == "Housekeeping"]
  
  # Run RUVg with the housekeeping genes and k dimensions of UV; K cannot exceed the number of housekeeping genes
  set <- RUVg(set, cIdx, k=k_fold) #you can experiment with other levels of k 
  
  # Build a data structure with DESeq2 for intermediate calculations
  dds <- DESeqDataSetFromMatrix(counts(set),colData=pData(set),design=~1)
  rowData(dds) <- featuredata
  
  # Estimate size factors by the median ratio method (Anders 2010, Genome Biology)
  dds <- estimateSizeFactors(dds)
  
  # Check for overdispersion
  dds <- estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  
  # Estimate dispersion per gene
  disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  
  # Assign dispersion to metadata of dds
  mcols(dds)$dispGeneEst <- disp
  
  # Correct for overdispersion
  dds <- estimateDispersionsFit(dds, fitType="mean")
  
  # Induce homoskedasticity based on fitted dispersion-mean relations
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  
  # Remove UV
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat <- removeBatchEffect(mat, covariates=covars)
  
  # Write out file of normalized data (if you need to write table, unquote)
  table_name <- paste(study_name,"_RUVg_k",as.character(k_fold),".txt",sep = "")
  write.table(mat, table_name, col.names = NA, sep="\t")
  
  
  # Bind with mat
  mat_cat <- rbind(mat,category)
  mat_cat <- mat_cat[c(nrow(mat_cat),1:(nrow(mat_cat)-1)),]
  data_name <- paste(study_name,"_RUVg_PCA_ER_k",as.character(k_fold),".txt",sep = "")
  write.table(mat_cat, data_name, col.names = NA, sep="\t")
  
  ### Above writes table of all genes, if want to do endo and hk separately, use the following
  ## Specify endogenous gene and HK gene if you want to do PCA plots separately
  # endo <- rownames(featuredata[which(featuredata$Gene_Class=="Endogenous"),])
  # hk <- rownames(featuredata[which(featuredata$Gene_Class=="Housekeeping"),])
  # 
  # mat_endo <- matcat[c(1,endo),]
  # data_name <- paste(study_name,"_Endo_RUVg_PCA_ER_k",as.character(k_fold),".txt",sep = "")
  # write.table(mat_endo, data_name, col.names = NA, sep="\t")
  # 
  # mat_hk <- matcat[c(1,hk),]
  # data_name <- paste(study_name,"_HK_RUVg_PCA_ER_k",as.character(k_fold),".txt",sep = "")
  # write.table(mat_hk, data_name, col.names = NA, sep="\t")
  
  
  ######Correlation Straight from Normalization
  if (!corr_plot){
  NormCounts_RUVg <- as.data.frame(t(mat))
  
  # levels(NormCounts_RUVg$ER.Status)
  
  BioCo_RUV <- cor.test(NormCounts_RUVg$ESR1,NormCounts_RUVg$GATA3, method = "spearman")[[4]]
  
  # Correlation plot
  corr_name <- paste("RUVg_corr_k",as.character(k_fold),".png",sep = "")
  png(filename=corr_name)
  
  p <- plot(NormCounts_RUVg$ESR1,NormCounts_RUVg$GATA3, pch=19,
            col=rgb(t(col2rgb(c( 'red', 'steelblue1','seagreen2'))), alpha=225, maxColorValue = 255)[NormCounts_RUVg$ER.Status], 
            main='Correlation Between Replication Genes',  cex.main = 1.2, mgp = c(3,.8,0),las=1,
            cex.axis = 1, cex.lab = 1.2, cex = 1.5, ylab="",
            xlab = expression(paste(italic(ESR1) , " (", log[2] , " normalized counts" , ")")))
  mtext(expression(paste(italic(GATA3) , " (", log[2] , " normalized counts" , ")")), 2, line = 2.5, cex = 1.2)
  legend("top",bty="n",paste('r = ', round(BioCo_RUV, digits = 2)))
  legend("right",bty="n", col=rgb(t(col2rgb(c('red', 'steelblue1','seagreen2'))), alpha=225, maxColorValue = 255), pch=19, legend=c("Negative","Positive","Unknown"))
  
  dev.off()
  }
  
  # PCA plot
  source("C:/Users/Amy/Desktop/Lab/Ghana/arrayTools_6_22.R") #PC
  
  ## Specify the plot you want to do (all/endo/hk)
  data_name <- paste(study_name,"_RUVg_PCA_ER_k",as.character(k_fold),".txt",sep = "")
                              # add _HK or _Endo here
  
  x<-readarray(data_name,hr=2) #first row is sample names, second row are the classes you are coloring by (for example: QC FLag or Tumor Subtype)
  
  pca_plot <- paste("RUVg_PCA_k",as.character(k_fold),".png",sep = "")
               # add _HK or _Endo here
  
  #Change numbers in "startPC" and "stopPC" to designate how many principle components you want to view. This example will display principle components 1-3
  loadings<- pca(x$xd,x$classes[,1],startPC=1,stopPC=2,
                 jitt=5,size=1,legendloc="topleft",returnLoadings=T)
  
  dev.copy(png,pca_plot)
  dev.off()
  
  pca_name <- paste("PC.Loadings_",study_name,"_RUVg_PCA_k",as.character(k_fold),".txt",sep = "")
                                         # add _HK or _Endo here
  write.table(loadings, pca_name, sep='\t', col.names=NA) # Read out of PCA scores
  
}

#### Things modified:
#   1. No need to change names for tables and plots, but need study name in the function
#   2. Need to specify if correlation plot is needed (default: no correlation plot)
#   3. No need to change sequence of mat_cat by hand
#   4. Add options for endo/hk
#################################################


###### Example ##################################
# Import ER status
ER_stat <- as.matrix(t(read.xlsx("C:/Users/Amy/Desktop/Lab/Ghana/Ghana_Pilot_ERstatus.xlsx", 2)))
colnames(ER_stat) <- ER_stat[1,]
ER_stat <- as.data.frame(ER_stat[2,,drop=F])

for(k in 1:6){
  pca_corr_kfold(Ghana_raw, Ghana_Col, Ghana_Genes,ER_stat,k)
}

