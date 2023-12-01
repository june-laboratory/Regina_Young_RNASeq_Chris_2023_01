library(tidyverse)
library(readxl)
library(data.table)
library(ggrepel)
library(edgeR)
library(DESeq2)
library(RColorBrewer)
library(EDASeq)
library(RColorBrewer)
library(data.table)
library(plyr)
library(plotrix)

#genecount_test <- readr::read_delim('/Users/linhuichen/Desktop/projects/Neil_Sheppard_RNASeq_Berjis_2022_12/data/summary_merged_count.bowtie.txt')
setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report')
genecount <-  read.csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/raw_counts.csv')
genecount <- genecount %>% column_to_rownames(var='ID')

#sample list 1 tumor/product ## 1801 = Δ133p53-P2A-CTL119 1924 = mCherry2-P2A-CTL119
#sample.list1 <-data.frame('sample'=colnames(genecount)[1:13],'groups'=c(rep('tumor',6),rep('product',6)),row.names=NULL)
sample.list2 <- data.frame('sample'=colnames(genecount)[1:6],'groups'=rep(c('mCherry2 tumor','Δ133 tumor'),3))
sample.list3 <- data.frame('sample'=colnames(genecount)[7:12],'groups'=rep(c('Δ133 product','mCherry2 product'),3))
sample.list <- rbind(sample.list2,sample.list3)
sample.list$groups <- as.factor(sample.list$groups)


## Building matrix 
c <- as.matrix(genecount[,1:12])
#row.names(c) <- genecount$Gene.name[which(row.names(c) == row.names(genecount))]
dds <- DESeqDataSetFromMatrix(c, sample.list, ~ groups)
dds <- DESeq(dds)
c <- DGEList(counts=genecount[,1:12], group=sample.list$groups)
c$counts <- counts(dds, normalized=TRUE)
c <- cpm(c,log=TRUE)
c <- as.data.frame(c)


c <- DGEList(counts=genecount[,1:12], group=sample.list$groups)
c$counts <- counts(dds, normalized=TRUE)
c <- cpm(c,log=TRUE)
c <- as.data.frame(c)
## ref and cmp

refs = c('Δ133 tumor', 'Δ133 product','Δ133 product', 'mCherry2 product')

cmps = c('mCherry2 tumor','mCherry2 product','Δ133 tumor','mCherry2 tumor')

## visualization

fdr.cutoff=0.05			# 5% FDR cutoff
logfc.cutoff=0



for(i in 1:length(refs)){	
  i=4
  ref=refs[i]
  cmp=cmps[i]
  res <- results(dds,contrast=c("groups",cmp,ref))
  if (min(na.omit(res$padj)) >fdr.cutoff){
    o <- order(abs(res$log2FoldChange),decreasing=TRUE)[1:5]
    x <- res$baseMean
    y <- res$log2FoldChange
    G <- row.names(res)
  
  }else{
    o <- which(res$padj <= fdr.cutoff)
    x <- res$baseMean
    y <- res$log2FoldChange
    G <- row.names(res)}
  dev.new()
  pdf(paste("idx DEGs ", paste(cmp, "vs", ref, sep="_"), "pdf",sep="."))
  
  plotMA(res, alpha=0.05, ylim=range(res$log2FoldChange[!is.na(res$log2FoldChange)]))
  res$gene <- genecount$Gene.name[which(row.names(res) == row.names(genecount))]
  #G[o] <- genecount$Gene.name[which(row.names(genecount)==G[o])]
  #text(x[o],y[o], labels=G[o],cex=0.4,col=490)
  hist(res[,"pvalue"], breaks=40, probability=T, plot=T, main="Distribution of P-Values idx", xlab="Gene P-values", ylab="Density", col="lightgrey")
  hist(res[,"padj"], breaks=40, probability=T, plot=T, main="Distribution of adjusted P-Values idx", xlab="Gene P-values", ylab="Density", col="lightgrey")
## volcano
## volcano plot 


# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > 1 & res$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res$diffexpressed[res$log2FoldChange < -1 & res$padj < 0.05] <- "DOWN"
  res$delabel <- NA
  res$gene <- genecount$Gene.name[which(row.names(res) == row.names(genecount))]
  res$delabel[res$diffexpressed != "NO"] <-res$gene[which(res$diffexpressed != "NO")]
  res1 <- as.data.frame(res)
  ggplot(data=res1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text_repel(size = 3, box.padding = 0.01,max.overlaps = 10) +
    ggtitle('abs(log2FoldChange) > 1  & padj < 0.05')
  dev.off()  
  res$DE.Call <- "Same"
  res[(res[,"padj"] <= fdr.cutoff) & (res[,"log2FoldChange"] >= logfc.cutoff) & !is.na(res$padj) & !is.na(res$log2FoldChange),"DE.Call"] <- "Up"
  res[(res[,"padj"] <= fdr.cutoff) & (res[,"log2FoldChange"] <= (-1 * logfc.cutoff)) & !is.na(res$padj) & !is.na(res$log2FoldChange),"DE.Call"] <- "Down"
  res <- as.data.frame(res)
  res$Gene <- rownames(res)
  res <- res[,c(8,1:7)]
  res <- cbind(res,c)
  DESeq2file = paste("DEGs", paste(cmp, "vs", ref, sep="_"), "csv",sep=".")	
  write.csv(res,file=DESeq2file,row.names=T)
}	



## volcano plot 
res <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/DEGs.mCherry2\ tumor_vs_mCherry2\ product.csv')
alpha <- 0.05
res$sig <- -log10(res$padj)
sum(is.infinite(res$sig))
res[is.infinite(res$sig),"sig"] <- 350
genes.to.plot <- !is.na(res$pvalue)
range(res[genes.to.plot, "log2FoldChange"])

cols <- densCols(res$log2FoldChange, res$sig)
cols[res$pvalue ==0] <- "purple"
res$pch <- 19
res$pch[res$pvalue ==0] <- 6
plot(res$log2FoldChange, 
     res$sig, 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=res$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")


gn.selected <- abs(res$log2FoldChange) > 5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.6)


