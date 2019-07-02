# primate_fibroblasts_DHS_OU_models

2019_06_17_modLog2FC_ouch.R - R commands to get moderated log2 fold change (R package "DESeq2") from DHS-seq data and run Ornstein-Uhlenbeck models (R package "ouch"). DESeq2 is used to transform count data into a continuous distribution for OU models. OU models are used to fit selective regimes on the tips of the tree. Best model is chosen by AIC comparisons and for significant deviation from Brownian motion. 
