library(limma)
library(isva) 

dir <- "C:/Users/Lily"

setwd(dir)

data.dir <- "C:/Users/Lily/"

BMIQ <- read.csv(paste(data.dir,"/BMIQ_PD_CT_seed5000-run1_3-20-16.csv", sep=""), check.names=FALSE, row.names=1)

info<-read.csv(paste(data.dir,"/RESULT/data_import_data/annotation.csv", sep=""))
info_PD <- info[which(info$Sample_Group=="PD"),]
BMIQ_PD<- BMIQ[, info_PD$Int]

duration <- info_PD$AOD-info_PD$AOO
batch <- factor(info_PD$Sentrix_ID)

## first remove batch and age effects

des_trt <- model.matrix(~duration, data=info_PD) 
methyl_adj <- removeBatchEffect(BMIQ_PD, batch = batch, batch2 = NULL, covariates = info_PD$AOD, design= des_trt)
#write.csv(methyl_adj, file="BMIQ_PD_adj_aod_batch.csv")

## iSVA
set.seed(5000)
isva.o <- DoISVA(data.m=methyl_adj, pheno.v=duration ,pvthCF=0.01,th=0.05)
isva.o$isv
library(Hmisc)
rcorr(duration, isva.o$isv[,1], type="spearman")
plot(duration, isva.o$isv[,1], ylab="first independent surrogate variable")

rcorr(duration, isva.o$isv[,2], type="spearman")
plot(duration, isva.o$isv[,2], ylab="second independent surrogate variable")

rcorr(duration, isva.o$isv[,3], type="spearman")
plot(duration, isva.o$isv[,3], ylab="third independent surrogate variable")

### output 
cf <- isva.o$isv
colnames(cf)<- c("IC1", "IC2", "IC3") 
cov <- cbind(info_PD, cf)

write.csv(cov, file="sample_with_ICs.csv")

