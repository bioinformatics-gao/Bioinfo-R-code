library(limma)
library(edgeR)
library(DESeq2)
library("pheatmap")


setwd("S:/pd/LW-PCA_calculation_7-19-2016")

targets = read.csv("Total_info_noPARK2.csv")
rownames(targets) <- targets[,11]
das <- readDGE(targets, skip=5, comment.char = "!")
keep <- rowSums(das$counts >= 10)>= 15 ##ten or more copies per individual in 15 or more individuals
TOT <- das[keep,]
dim(TOT) # 19320   31


pdf(file="pcaplots_heatmap.pdf", paper="USr")


dds <- DESeqDataSetFromMatrix(countData=data.frame(TOT$counts), colData=TOT$samples, design=~1)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup="Batch2")
data <- plotPCA(vsd, intgroup=c("Batch2","Tissue","Description","Day"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Batch2, shape=Tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

ggplot(data, aes(PC1, PC2, color=Batch2, shape=Description)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

ggplot(data, aes(PC1, PC2, color=Batch2, shape=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

ggplot(data, aes(PC1, PC2, color=Description, shape=Tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

ggplot(data, aes(PC1, PC2, color=Description, shape=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



## Heatmap of the count matrix
dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Batch2","Tissue")])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df)



graphics.off()



