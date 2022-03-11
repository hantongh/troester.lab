########################## 
##   This program is designed to clean and organize raw data from nSolver
##   Required data sets: raw count, phenodata(QC,content flags etc.)
##   Output: Missingness set to 0; log2 of Raw counts; organized phenoData (with flags and outliers)
##
##   Functions:
##   1. missingness(raw_data)
##        Purpose: Missingness set to 0
##                Calculates the average of Negative genes and compare it to endogenous and HK genes
##                If Endo & HK < Avg NEG, set to zero. Count # and percentage of zeros
##   2. Set raw count to log2
##   3. Find outlier for log2
##        Can be endo or HK, can exclude unwanted samples
##   4. Organize PhenoData Flag
##########################

### All outputs can append to existing sheet or write to different file
## Append
# write.xlsx(output, original_file, sheetName = "Missing_Set0", append=TRUE)
## Write new file
# write.xlsx(output, file)
####################################

### Missingness ####################################
## Input: (data frame) raw data
## raw data format:
# 
#              Class_name         Sample1  Sample2 ...
# Probe1      Endogenous           xxx
# Probe2      Endogenous           xxx
# Probe3      Housekeeping         xxx
# Probe4      Negative             xxx
# ...         ...                  ...
# 
## Output: (data frame) raw data with missingness set to 0
## Syntax: output <- missingness(input)

missingness <- function(raw_data){
  ######### Missingness set to 0 #########
  ### 1. Find average NEG
  # Use the raw count data here
  raw <- raw_data
  
  # Calculate average for NEG
  neg <- raw[which(raw$Class_Name=="Negative"),c(-1)]
  average_neg <- apply(neg, 2, mean)
  
  # Need to keep element type as double; also need column number equal to raw data
  # Add 2 NAs before average
  average_neg <- c(NA,average_neg)
  
  # Combine average to raw
  addAVG <- rbind(raw,"Average"=average_neg)
  ###
  
  ### 2. Set Endo & HK < Avg NEG to 0
  for(i in c(2:ncol(addAVG))){
    addAVG[c(1:nrow(raw)),i][addAVG[c(1:nrow(raw)),i]<average_neg[i]] <- 0
  }
  ###
  
  ### 3. Count missingness
  # Split endogenous and HK
  
  # Endogenous count
  endo <- addAVG[which(addAVG$Class_Name=="Endogenous"),c(-1)]
  endo_count <- apply(endo, 2, function(x) sum(x==0))
  endo_count <- c(NA,endo_count)
  # percent
  endo_pct <- endo_count/nrow(endo)*100
  
  # HK count
  hk <- addAVG[which(addAVG$Class_Name=="Housekeeping"),c(-1)]
  hk_count <- apply(hk, 2, function(x) sum(x==0))
  hk_count <- c(NA,hk_count)
  # percent
  hk_pct <- hk_count/nrow(hk)*100
  ###
  
  ### 4. Bind all together
  Missingness_setZero <- rbind(addAVG,"Num_EndM"=endo_count,"PCT_EndM"=endo_pct,
                               "Num_HKM"=hk_count,"PCT_HKM"=hk_pct)
  
  return(Missingness_setZero)
}
#####################################################################

### Set log2 #####################

## Input: (data frame) raw data
## raw data format:
# 
#              Class_name         Sample1  Sample2 ...
# Probe1      Endogenous           xxx
# Probe2      Endogenous           xxx
# Probe3      Housekeeping         xxx
# Probe4      Negative             xxx
# ...         ...                  ...
# 
## Output: (data frame) log2 for every numeric value
## Syntax: output <- set_log2(input)

### Set all counts to log2
set_log2 <- function(raw_data){
  
  # Duplicate raw data
  raw_log2 <- raw_data
  
  # Use loop to log2() each value
  for(i in c(2:ncol(raw))){
    for(j in c(1:nrow(raw))){
      raw_log2[j,i] <- log2(raw[j,i])
    }
  }
  
  return(raw_log2)
}
######################################


### Calculate outliers ##################################
## Input: (data frame) Log2 data, (string) specify outlier for hk or endo
## log2 data format:
# 
#              Class_name         Sample1  Sample2 ...
# Probe1      Endogenous           xx.xx
# Probe2      Endogenous           xx.xx
# Probe3      Housekeeping         xx.xx
# Probe4      Negative             xx.xx
# ...         ...                  ...
# 
## Output: (data frame) Subset of log2 data with median and outlier flag
## Syntax: output <- outlier(data = input, class_name = "")  If need HK outlier, class = "Housekeeping". If need Endo outlier, class = "Endogenous"

outlier <- function(log2data, class_name){
  raw_log2 <- log2data
  
  log2_sub <- as.data.frame(t(raw_log2[which(raw_log2$Class.Name==class_name),c(-1)]))
    
  ## If need to exclude any samples, specify here
  # log2_sub[c(unwanted_samples)] <- NULL
    
  # Find median of log2
  median <- apply(log2_sub, 1, median)
    
  # quantile function output a matrix like this
  # quartile  0%        25%          50%      ... 
  # data      min   1st quartile    median     ...
    
  q1 <- quantile(median)[2]               # Find 1st and 3rd quartile
  q3 <- quantile(median)[4]
  
  # Check if it is outlier, returns list
  # function: outlier if median<q1-1.5*(q3-q1) or median>q3+1.5*(q3-q1)
  outlier_check <- lapply(median, function(x) (x<q1-1.5*(q3-q1)) | (x>q3+1.5*(q3-q1)))
  outlier <- as.integer(unlist(outlier_check))
  # Bind gene, median, outlier (True or false) into data set
  sample_outlier <- cbind(log2_sub,Median=median,Outlier=outlier)
  
  return(sample_outlier)
}

#######################################################################################


### PhenoData Flag ############################################
## Input: (data frame) Phenotype data, Missingness data, Outlier data, HKM threshold
## phenotype data format:
# 
#              QC.Flag    Blank.Lane.Flag  mRNA.Positive.Normalization.Flag ...
# Sample1       YES           NO
# Sample2       NO            NO
# Sample3       NO            NO
# Sample4       YES           NO
# ...           ...          ...
# 
## missingness data format:
# 
#              Class_name     Sample1    Sample2 ...
# ...           ...             ...
# Average       NA             xx.xx
# Num_EndM      NA              0
# PCT_EndM      NA              0
# Num_HKM       NA              1
# PCT_HKM       NA             9.09
# 
## Output: (data frame) Original Phenotype data with
## Syntax: output <- pheno_flag(pheno = pheno_data, missingness = missingness_data, outlier = outlier_data, hkm_threshold = xxx)  HKM>xxx, then HKM flag, by default xxx=0

pheno_flag <- function(pheno, missingness, outlier, hkm_threshold=0){
  
  ## 1. Original flag (YES/NO value) into 0/1
  pheno$QC.Flag <- as.integer(pheno$QC.Flag)-1
  pheno$Blank.Lane.Flag <- as.integer(pheno$Blank.Lane.Flag)-1
  pheno$mRNA.Positive.Normalization.Flag <- as.integer(pheno$mRNA.Positive.Normalization.Flag)-1
  pheno$mRNA.Content.Normalization.Flag <- as.integer(pheno$mRNA.Content.Normalization.Flag)-1
  
  ## 2. NEG, HKM
  # Use Missingness data produced from missingness function
  neghkm <- t(missingness[c((nrow(missingness)-3):nrow(missingness)),c(-1)]) # Use the last 4 rows and add to pheno data
  pheno <- cbind(pheno,neghkm)
  
  # Check if it is neg flag, return list
  # function: neg flag if pct_end > 30
  neg_check <- lapply(pheno$PCT_EndM, function(x) x>30) 
  neg_flag <- as.integer(unlist(neg_check))
  
  # Check if pct_hk>0 and pct_end out of 2 sd
  hkm_check <- lapply(pheno$PCT_HKM, function(x) x>hkm_threshold)
  hkm_flag <- as.integer(unlist(hkm_check))
  
  pheno <- cbind(pheno, NEG_Flag=neg_flag, HKM_Flag=hkm_flag)
  
  ## 3. Outlier (can omit this step if didn't check outlier)
  pheno <- cbind(pheno, Outlier=outlier)
  
  ## 4. Conclude flag
  flag <- rep("Pass",times = nrow(pheno)) # First assume every sample is "PASS"
  
  # If sample has EQC flag, set flag to EQC (empty flag); if there is 1 or more than 1 flag for the other flags, set flag to "" (empty flag)
  for (i in c(1:nrow(pheno))){
    if (pheno$QC.Flag[i]==1){
      flag[i] <- "EQC"
    } else if (pheno$Blank.Lane.Flag[i]+pheno$mRNA.Positive.Normalization.Flag[i]+
               pheno$mRNA.Content.Normalization.Flag[i]+pheno$NEG_Flag[i]+pheno$HKM_Flag[i]+pheno$HK_Outlier[i]>0){
      flag[i] <- ""  # If outlier is not checked, delete +pheno$HK_Outlier[i] from previous line
    }
  }
  
  # Check other flags, separate different flags by .
  for (i in which(flag=="")){
    if (pheno$Blank.Lane.Flag[i]==1){
      flag[i] <- paste(flag[4],"Blane",sep = ".")
    }
    
    if (pheno$mRNA.Positive.Normalization.Flag[i]==1){
      flag[i] <- paste(flag[i],"POS",sep = ".")
    }
    
    if (pheno$mRNA.Content.Normalization.Flag[i]==1){
      flag[i] <- paste(flag[i],"Content",sep = ".")
    }
    
    if (pheno$NEG_Flag[i]==1){
      flag[i] <- paste(flag[i],"NEG",sep = ".")
    }
    
    if (pheno$HKM_Flag[i]==1){
      flag[i] <- paste(flag[i],"HKM",sep = ".")
    }
    
    # Omit following step if outlier is not checked
    if (pheno$HK_Outlier[i]==1){
      flag[i] <- paste(flag[i],"Outlier",sep = ".")
    }
    
    if (substr(flag[i],1,1)=="."){
      flag[i] <- substring(flag[i],2)
    }
  }
  
  # Combine with pheno data set
  pheno <- cbind(pheno,Flag=flag)
  
  return(pheno)
}

#########################################################################################























