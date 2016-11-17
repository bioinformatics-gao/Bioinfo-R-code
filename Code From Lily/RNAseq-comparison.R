# DMV HIHG RNA-Seq analysis Duration studies.

library(limma)
library(edgeR)
location <-"C:/Users/Lily"
setwd (location)


x <- list()
x$samples <- read.delim("DMVR_PD_info.txt") 
x$counts <- read.delim("DMVR_PD_RAW.txt")
PD <- x
PDdes <- model.matrix(~DIFF+AOD, data = PD$samples)
PDdge <- DGEList(counts= PD$counts)
PDdge <- calcNormFactors(PDdge)
PDv <- voom(PDdge, PDdes)
PDfit <- lmFit(PDv, PDdes)
PDfit <- eBayes(PDfit)

# duration of disease
res_duration <- topTable(PDfit, coef=2, adjust.method="BH", n=Inf, sort.by="P")
write.csv(res_duration, file="DMV_duration.csv")

# aod 
res_aod <- topTable (PDfit, coef=3, adjust.method="BH", n=Inf, sort.by="P")
write.csv(res_aod, file="DMV_aod.csv")

# file with both duration and aod
res_both <-cbind(PDfit$t[,2], PDfit$p.value[,2], PDfit$t[,3], PDfit$p.value[,3])
colnames(res_both)<-c("DIFF_t", "DIFF_p", "AOD_t", "AOD_p")
write.csv(res_both, file="DMV_duration_aod.csv")

