library(data.table)
library(limma)

meta_data<-read.csv(file="whole_blood_meta_data.csv", header=T, colClasses=c("character", "factor", "numeric"))
raw_data<-fread(file="whole_blood_raw_counts.csv")

ensembl_id<-raw_data$Ensembl_id
raw_data$Ensembl_id<-NULL
raw_data<-as.matrix(raw_data)
rownames(raw_data)<-ensembl_id

group<-group<-meta_data$NLR_category
design<-model.matrix(~ 0 + group)

v<-voom(raw_data, design=design , lib.size=meta_data$Library_size, normalize.method = "none")
fit<-lmFit(v, design)

contrast<-makeContrasts(groupHi-groupLow, levels=colnames(fit$coefficients))

result<-contrasts.fit(fit, contrast)
result<-eBayes(result)
result_table<-topTable(result, n=Inf)

write.csv(result_table, file="whole_blood_limma_result.csv")


