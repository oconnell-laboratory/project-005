library(data.table)
library(DESeq2)

meta_data<-read.csv(file="isolated_wbc_meta_data.csv", header=T, colClasses=c("character", "factor"))
rownames(meta_data)<-meta_data$Sample_id
meta_data$Sample_id<-NULL

data<-fread(file="isolated_wbc_tpm_values.csv")
ensembl_id<-data$Ensembl_id
data$Ensembl_id<-NULL
data<-as.matrix(data)
rownames(data)<-ensembl_id

data<-round(data*10000, 0)
mode(data)<-"integer"

DEseqDS<-DESeqDataSetFromMatrix(countData=data, colData=meta_data, design = ~ Cell_type)
sizeFactors(DEseqDS)<-rep(10000, dim(data)[2])
DEseqDS<-DESeq(DEseqDS)
results<-results(DEseqDS, independentFiltering=F)

write.csv(results, file="isolated_wbc_DESeq2_result.csv")

