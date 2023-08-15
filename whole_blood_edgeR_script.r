library(data.table)
library(edgeR)

meta_data<-read.csv(file="whole_blood_meta_data.csv", header=T, colClasses=c("character", "factor", "numeric"))
raw_data<-fread(file="whole_blood_raw_counts.csv")

ensembl_id<-raw_data$Ensembl_id
raw_data$Ensembl_id<-NULL
raw_data<-as.matrix(raw_data)
rownames(raw_data)<-ensembl_id

design<-model.matrix(~meta_data$NLR_category)

dge<-DGEList(counts=raw_data, lib.size=meta_data$Library_size, norm.factors=rep(1,ncol(raw_data)), samples=meta_data$Sample_id, group=meta_data$NLR_category, remove.zeros=F)
dge<-estimateDisp(y=dge, design=design)
result<-exactTest(dge, pair=c(2,1))
result_table<-topTags(result, n=Inf)

write.csv(result_table, file="whole_blood_edgeR_result.csv")

