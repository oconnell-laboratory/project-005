library(data.table)
library(DESeq2)

meta_data<-read.csv(file="whole_blood_meta_data.csv", header=T, colClasses=c("character", "factor", "numeric"))
rownames(meta_data)<-meta_data$Sample_id
meta_data$Sample_id<-NULL

raw_data<-fread(file="whole_blood_raw_counts.csv")
ensembl_id<-raw_data$Ensembl_id
raw_data$Ensembl_id<-NULL
raw_data<-as.matrix(raw_data)
rownames(raw_data)<-ensembl_id

DEseqDS<-DESeqDataSetFromMatrix(countData=raw_data, colData=meta_data, design = ~ NLR_category)
sizeFactors(DEseqDS)<-meta_data$Library_size/min(meta_data$Library_size) 
DEseqDS<-DESeq(DEseqDS)
results<-results(DEseqDS, independentFiltering=F)

write.csv(results, file="whole_blood_DESeq2_result.csv")