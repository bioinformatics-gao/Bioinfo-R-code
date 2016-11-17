library(RnBeads)
library(limma)
library(edgeR)
#library(FDb.InfiniumMethylation.hg19 )

##################### DMV PREPROCESSING ######################
data.dir <- "C:/Users/Lily"
setwd(data.dir)
idat.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "DMV_Old.csv")
report.dir <- file.path(data.dir, "RESULT")
data.source <- c(idat.dir, sample.annotation)

########### 1. import 
rnb.options(identifiers.column="Int")
result <- rnb.run.import(data.source=data.source, data.type="infinium.idat.dir", dir.reports=report.dir)
rnb.set <- result$rnb.set

########### 2. filter
any.bad.p.val <- apply(dpval(rnb.set)>0.01, 1, any)
rnb.set2 <- remove.sites(rnb.set, any.bad.p.val)
nrow(meth(rnb.set2))
head(meth(rnb.set2, type = "sites", row.names = TRUE))

########## 3. normalization
set.seed (5000)
rnb.set2.bmiq <- rnb.execute.normalization(rnb.set2, method="bmiq", bgcorr.method="methylumi.noob")
BMIQ <- mval(rnb.set2.bmiq, type = "sites", row.names = TRUE)
head (BMIQ)

write.csv(BMIQ, file="BMIQ_PD_CT_seed5000-run2_3-20-16.csv")
#BMIQ <- read.csv(paste(data.dir,"/BMIQ_PD_CT_seed5000-run1_3-20-16.csv", sep=""),check.names=FALSE, row.names=1)

########## 4. read in sample info
info<-read.csv(paste(report.dir,"/data_import_data/annotation.csv", sep=""))

########## 5. linear model: M~ duration + aod + batch
info_PD <- info[which(info$Sample_Group=="PD"),]
duration <- info_PD$AOD-info_PD$AOO
batch <- factor(info_PD$Sentrix_ID)

design = model.matrix(~duration+AOD+batch, data = info_PD)

BMIQ_PD<- BMIQ[, info_PD$Int]

# linear model without empirical bayes moderation

lmF=function(methylation){  

  tmp=coef(summary(lm(methylation~AOD+duration+batch, data=info_PD)))
  return(tmp[2:3,]) 
}


res<-t(apply(BMIQ_PD,1, lmF))

result_AOD<-res[,c(1,3,5,7)]
colnames(result_AOD)<-c("Estimate", "Std", "tvalue", "pvalue")
result_duration<-res[,c(2,4,6,8)]
colnames(result_duration)<-c("Estimate", "Std", "tvalue", "pvalue")

write.csv(result_AOD, "result_LM_AOD_byR_seed5000_3-20-2016.csv")
write.csv(result_duration, "result_LM_duration_byR_seed5000_3-20-2016.csv")


########## 5. linear model: M ~ group +  duration + aod + batch + group x aod
batch <- factor(info$Sentrix_ID)
desig_dms <- model.matrix(~Sample_Group + AOD + batch + AOD*Sample_Group, data=info)

#small <- BMIQ [ which(rownames(BMIQ)=="cg07158503") ,]
#coef(summary(lm(as.numeric(small) ~ Sample_Group + AOD + batch + AOD*Sample_Group, data=info)))

lmF_dms =function(methylation){  
  tmp=coef(summary(lm(methylation~Sample_Group + AOD + batch + AOD*Sample_Group, data=info)))
  return(tmp[2,]) 
}

res_dms <- t(apply(BMIQ,1, lmF_dms))

write.csv(res_dms, "result_dms_byR_seed5000_3-20-2016.csv")

###### add FDRs
dms<- read.csv(paste(data.dir,"/result_LM_dms_byR_seed5000_3-20-2016.csv", sep=""))
result_dms<-cbind (dms, p.adjust(dms[ ,5], method="fdr") )
colnames(result_dms)<-c("ProbeID", "Estimate", "Std", "tvalue", "pvalue", "fdr")

write.csv (result_dms, "result_LM_dms_byR_seed5000_FDR_3-20-2016.csv")
