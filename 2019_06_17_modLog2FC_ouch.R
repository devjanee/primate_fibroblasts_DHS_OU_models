# 2019_06_18
# getting log2foldchange from DESeq2 in respect to Macaque
# log2foldchange for Macaque will always be 0

# load library

library(DESeq2)

# get traits
traits<-read.table("~/Documents/fibroblast/all_DHS_sites.glm_analysis.all_sites.all_information.txt",header=T)

# get raw counts
mydata<-traits[,44:58]

# make coldata
coldata<-data.frame(matrix(ncol = 1, nrow = 15))
rownames(coldata)<-colnames(mydata)
colnames(coldata)<-"species"
coldata[1:3,]<-"h"
coldata[4:6,]<-"c"
coldata[7:9,]<-"g"
coldata[10:12,]<-"o"
coldata[13:15,]<-"m"

# make DESeq data set (dds)
dds <- DESeqDataSetFromMatrix(countData = mydata, colData = coldata, design = ~ species)

# run DESeq to get log2FC
dds<-DESeq(dds)

# get moderated log2fold change

res_H <- lfcShrink(dds, contrast=c("species","h","m"))
res_C <- lfcShrink(dds, contrast=c("species","c","m"))
res_G <- lfcShrink(dds, contrast=c("species","g","m"))
res_O <- lfcShrink(dds, contrast=c("species","o","m"))

# combine the results and change column names

res_H.df<-as.data.frame(res_H$log2FoldChange)
res_C.df<-as.data.frame(res_C$log2FoldChange)
res_G.df<-as.data.frame(res_G$log2FoldChange)
res_O.df<-as.data.frame(res_O$log2FoldChange)
modLog2FC <- cbind(res_H.df,res_C.df)
modLog2FC <- cbind(modLog2FC,res_G.df)
modLog2FC <- cbind(modLog2FC,res_O.df)
modLog2FC$M <-0
colnames(modLog2FC)<-c("H","C","G","O","M")

write.table(modLog2FC,"~/Documents/fibroblast/2019_06_17_moderated_log2FC_DESeq2.txt",col.names = T, sep="\t",quote = F)

# filter out traits significant by glm gateway test
modLog2FC_sig<-modLog2FC[which(traits$gateway_adjusted_pvalue<0.01),]

#################################
# Use moderated log2FC in ouch  #
#################################

# load libraries

library(ape)
library(ouch)

# read in tree

primate.tree <-read.nexus("~/Downloads/consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")

# change ape tree to ouch tree
primate.tree.ouch<-ape2ouch(primate.tree,scale = F,branch.lengths = primate.tree$edge.length)

df_O <- data.frame(matrix(ncol = 8, nrow = 0))
x <- c("bm_sigsq","bmAIC","ou_sigsq","ou_alpha","ouAIC","lbm","lou","ou_bm")
colnames(df_O) <- x

for(i in 1:nrow(modLog2FC_sig)){
  print(i)
  test<-as.numeric(modLog2FC_sig[i,c(3,1,2,4,5)])
  names(test)<-primate.tree$tip.label
  ERreconstruction <- ace(test, primate.tree, type="continuous", model="ER")
  test.ouch<-c(test,ERreconstruction$ace)
  names(test.ouch)<-c("9","8","7","6","5","1","4","3","2")
  
  # run BM and OU models
  
  test.ouch.bm<-brown(test.ouch,primate.tree.ouch)
  #regimes_G<-as.factor(c("1","1","1","1","1","1","1","1","G"))
  # for each branch, define a different selection regime
  # regimes_C<-as.factor(c("1","1","1","1","1","1","C","1","1"))
  # regimes_H<-as.factor(c("1","1","1","1","1","1","1","H","1"))
   regimes_O<-as.factor(c("1","1","1","1","1","O","1","1","1"))
  # This should be more efficiently done in one for loop, but I didn't so here we are
  
  names(regimes_O)<-myoudataframe$nodes
  test.ouch.ou<-hansen(data = test.ouch,tree = primate.tree.ouch,sqrt.alpha = 1, sigma = 1,method="L-BFGS-B",regimes = regimes_O)
  
  # get sigsq, z0, AIC, AICc, and likelihoods from bm
  bm_sigsq <- summary(test.ouch.bm)$sigma.squared
  bm_alpha <- summary(test.ouch.bm)$alpha
  bmAIC <- summary(test.ouch.bm)$aic
  lbm<-summary(test.ouch.bm)$loglik
  
  # get sigsq, z0, AIC, AICc, and likelihoods from ou
  ou_sigsq <- summary(test.ouch.ou)$sigma.squared
  ou_alpha <- summary(test.ouch.ou)$alpha
  ouAIC <- summary(test.ouch.ou)$aic
  lou<-summary(test.ouch.ou)$loglik
  
  # Are the likelihoods different according to chi-sq?
  ou_bm<- pchisq(2*(lou-lbm), df=1, lower.tail=F)
  
  # populate the data frame
  df_O[i,]<-c(bm_sigsq,bmAIC,ou_sigsq,ou_alpha,ouAIC,lbm,lou,ou_bm)
}

write.table(df_G,"~/Documents/fibroblast/2019_06_17_ouch_G.txt",sep="\t",quote = F,col.names = T)
# write all the tables

##################################
# Filter data and compare to glm #
##################################

# pick the lowest AIC

for(i in 1:nrow(df)){if((which.min(c(df[i,2],df[i,6],df[i,11],df[i,16],df[i,21])))==1){df[i,24]="BM"}
                  else if((which.min(c(df[i,2],df[i,6],df[i,11],df[i,16],df[i,21])))==2){df[i,24]="H"}
                  else if((which.min(c(df[i,2],df[i,6],df[i,11],df[i,16],df[i,21])))==3){df[i,24]="C"}
                  else if((which.min(c(df[i,2],df[i,6],df[i,11],df[i,16],df[i,21])))==4){df[i,24]="G"}
                  else if((which.min(c(df[i,2],df[i,6],df[i,11],df[i,16],df[i,21])))==5){df[i,24]="O"}
                  else({df[i,24]="undetermined"})
}

# see if selective regime is significantly different than BM. If not - then "undetermined"

for(i in 1:nrow(df)){if(df[i,24]=="H"&&df[i,]$ou_bm_H<0.01){df[i,25]="H"}
  else if(df[i,24]=="C"&&df[i,]$ou_bm_C<0.01){df[i,25]="C"}
  else if(df[i,24]=="G"&&df[i,]$ou_bm_G<0.01){df[i,25]="G"}
  else if(df[i,24]=="O"&&df[i,]$ou_bm_O<0.01){df[i,25]="O"}
  else if(df[i,24]=="BM"){df[i,25]="BM"}
  else({df[i,25]="undetermined"})
}

for(i in 1:nrow(df)){if(df[i,]$significant_AIC=="BM"){df[i,]$classification="BM"}
  else if(df[i,]$significant_AIC=="H"&&df[i,]$H>0){df[i,]$classification="H_gain"}
  else if(df[i,]$significant_AIC=="H"&&df[i,]$H<0){df[i,]$classification="H_loss"}
  else if(df[i,]$significant_AIC=="C"&&df[i,]$C>0){df[i,]$classification="C_gain"}
  else if(df[i,]$significant_AIC=="C"&&df[i,]$C<0){df[i,]$classification="C_loss"}
  else if(df[i,]$significant_AIC=="G"&&df[i,]$G>0){df[i,]$classification="G_gain"}
  else if(df[i,]$significant_AIC=="G"&&df[i,]$G<0){df[i,]$classification="G_loss"}
  else if(df[i,]$significant_AIC=="O"&&df[i,]$O>0){df[i,]$classification="O_gain"}
  else if(df[i,]$significant_AIC=="O"&&df[i,]$O<0){df[i,]$classification="O_loss"}
  else({df[i,]$classification="undetermined"})
}



