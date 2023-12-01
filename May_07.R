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




##### May 06 #######
####3 product 
setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report')
genecount <-  read.csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/raw_counts.csv')
genecount <- genecount %>% column_to_rownames(var='ID')

#sample list 1 tumor/product ## 1801 = Δ133p53-P2A-CTL119 1924 = mCherry2-P2A-CTL119
#sample.list1 <-data.frame('sample'=colnames(genecount)[1:13],'groups'=c(rep('tumor',6),rep('product',6)),row.names=NULL)
sample.list2 <- data.frame('sample'=colnames(genecount)[1:6],'groups'=rep(c('mCherry2 tumor','Δ133 tumor'),3))
sample.list3 <- data.frame('sample'=colnames(genecount)[7:12],'groups'=rep(c('Δ133 product','mCherry2 product'),3))
sample.list <- rbind(sample.list2,sample.list3)
sample.list$groups <- as.factor(sample.list$groups)

product <- c('Δ133 product','mCherry2 product')

sample.list.g <- sample.list[sample.list$groups %in% product,]

c <- as.matrix(genecount[,sample.list$sample[sample.list$groups %in% product]])
#row.names(c) <- genecount$Gene.name[which(row.names(c) == row.names(genecount))]
dds <- DESeqDataSetFromMatrix(c, sample.list.g, ~ groups)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
c$counts <- counts(dds, normalized=TRUE)
write.csv(c$counts, file=paste(cmp, "vs", ref, ".norm_counts.csv",sep="_"))


cmps = c('Δ133 product')

refs = c('mCherry2 product')


fdr.cutoff=0.05			# 5% FDR cutoff
logfc.cutoff=0



for(i in 1:length(refs)){	
  i=1
  ref=refs[i]
  cmp=cmps[i]
  res <- results(dds,contrast=c("groups",cmp,ref))
  res <- res[is.na(res$padj)==FALSE,]
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
  res$gene <- genecount$Gene.name[match(row.names(res),row.names(genecount))]
  #G[o] <- genecount$Gene.name[which(row.names(genecount)==G[o])]
  #text(x[o],y[o], labels=G[o],cex=0.4,col=490)
  hist(res[,"pvalue"], breaks=40, probability=T, plot=T, main="Distribution of P-Values idx", xlab="Gene P-values", ylab="Density", col="lightgrey")
  hist(res[,"padj"], breaks=40, probability=T, plot=T, main="Distribution of adjusted P-Values idx", xlab="Gene P-values", ylab="Density", col="lightgrey")
  ## volcano
  ## volcano plot 
  
  
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res$diffexpressed[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
  res$delabel <- NA
  res$delabel[res$diffexpressed != "NO"] <-res$gene[which(res$diffexpressed != "NO")]
  res1 <- as.data.frame(res)
  ggplot(data=res1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text_repel(size = 3, box.padding = 0.01,max.overlaps = 10) +
    ggtitle('abs(log2FoldChange) > 0  & padj < 0.05')
  dev.off()  
  res$DE.Call <- "Same"
  res[(res[,"padj"] <= fdr.cutoff) & (res[,"log2FoldChange"] >= logfc.cutoff) & !is.na(res$padj) & !is.na(res$log2FoldChange),"DE.Call"] <- "Up"
  res[(res[,"padj"] <= fdr.cutoff) & (res[,"log2FoldChange"] <= (-1 * logfc.cutoff)) & !is.na(res$padj) & !is.na(res$log2FoldChange),"DE.Call"] <- "Down"
  res <- as.data.frame(res)
  res$Gene <- rownames(res)
  res <- res[,c(8,1:7)]
  res <- cbind(res,c$counts[match(row.names(res),row.names(c$counts)),])
  DESeq2file = paste("DEGs", paste(cmp, "vs", ref, sep="_"), "csv",sep=".")	
  write.csv(res,file=DESeq2file,row.names=T)
}	


####3 tumor
setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report')
genecount <-  read.csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/raw_counts.csv')
genecount <- genecount %>% column_to_rownames(var='ID')

#sample list 1 tumor/product ## 1801 = Δ133p53-P2A-CTL119 1924 = mCherry2-P2A-CTL119
#sample.list1 <-data.frame('sample'=colnames(genecount)[1:13],'groups'=c(rep('tumor',6),rep('product',6)),row.names=NULL)
sample.list2 <- data.frame('sample'=colnames(genecount)[1:6],'groups'=rep(c('mCherry2 tumor','Δ133 tumor'),3))
sample.list3 <- data.frame('sample'=colnames(genecount)[7:12],'groups'=rep(c('Δ133 product','mCherry2 product'),3))
sample.list <- rbind(sample.list2,sample.list3)
sample.list$groups <- as.factor(sample.list$groups)

tumor <- c('Δ133 tumor','mCherry2 tumor')

sample.list.g <- sample.list[sample.list$groups %in% tumor,]

c <- as.matrix(genecount[,sample.list$sample[sample.list$groups %in% tumor]])
#row.names(c) <- genecount$Gene.name[which(row.names(c) == row.names(genecount))]
dds <- DESeqDataSetFromMatrix(c, sample.list.g, ~ groups)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
c$counts <- counts(dds, normalized=TRUE)
write.csv(c$counts, file=paste(cmp, "vs", ref, ".norm_counts.csv",sep="_"))



c1 <- as.data.frame(c$counts[keep,])



cmps = c('Δ133 tumor')

refs = c('mCherry2 tumor')


fdr.cutoff=0.05			# 5% FDR cutoff
logfc.cutoff=0



for(i in 1:length(refs)){	
  i=1
  ref=refs[i]
  cmp=cmps[i]
  res <- results(dds,contrast=c("groups",cmp,ref))
  res <- res[is.na(res$padj)==FALSE,]
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
  res$gene <- genecount$Gene.name[match(row.names(res) ,row.names(genecount))]
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
  res <- cbind(res,c$counts[match(row.names(res),row.names(c$counts)),])
  DESeq2file = paste("DEGs", paste(cmp, "vs", ref, sep="_"), "csv",sep=".")	
  write.csv(res,file=DESeq2file,row.names=T)
}	



deg1 <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/Mar_27/DEGs.Δ133\ product_vs_mCherry2\ product.csv')
deg2 <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_07/DEGs.Δ133\ product_vs_mCherry2\ product.csv')
deg2 <- deg2[order(deg2$padj,decreasing=FALSE),]
deg1 <- na.omit(deg1)
deg1[deg1$padj<=0.05,]
deg2[deg2$padj<=0.05,]
length(deg1[deg1$pvalue<= 0.05,])
length(deg2[deg2$pvalue<= 0.05,])


tumor.a <- read_csv('/Users/linhuichen/Downloads/20230219_Differential_expression_analysis_table.csv')
tumor <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_07/DEGs.Δ133\ tumor_vs_mCherry2\ tumor.csv')
length(tumor.a[tumor.a$padj<= 0.05,])
length(tumor[tumor$padj<= 0.05,])

tumor.a[order(tumor.a$padj),]
tumor.a <- tumor.a[match(tumor$...1,tumor.a$ID),]


#### ### May 24 ### ####


setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24')
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

refs = c('mCherry2 product','Δ133 product')

cmps = c('mCherry2 tumor')

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


###########     GSEA ########

f.list <- list.files(pattern='.csv')
deg <- read_csv(f.list[4])

## All gene ranking
alldiff <- deg[order(deg$log2FoldChange,decreasing = T),]
id1 <- alldiff$log2FoldChange
names(id1) <- alldiff$...1
id1 <- na.omit(id1)
id1 <- sort(id1,decreasing=TRUE)



geneset <- read_excel('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/geneSets_d133.xlsx')
genelist <- c('term'=c(),'gene'=c())
for(j in colnames(geneset)){
  a <- rep(j,nrow(geneset))
  b <- geneset[,j]
  c <- cbind(a,b)
  colnames(c) <- c('term','gene')
  genelist <- rbind(c,genelist)
}

genelist <- na.omit(genelist)
gsea.re1<-clusterProfiler::GSEA(id1,TERM2GENE=genelist,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,pvalue<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

genelist.list <-genelist%>% split(.$term)  %>% lapply( "[[", 2)

gsea.re2 <- fgsea(pathways =  genelist.list,
                  stats = id1,
                  minSize=1,
                  maxSize=10000,
                  eps=0
)
colnames(gsea.re2)

#
g2 <- gsea.re2[gsea.re2$pval<0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]



dev.new()
pdf(str_replace_all(f.list[4],'csv','gsea.pdf'))
plotGseaTable(genelist.list[g2$pathway],
              id1, 
              gsea.re2,gseaParam = 0.5,
              colwidths = c(0.5,0.2,0.1,0.1,0.1)
)



library(ggsci)
library(enrichplot)
col_gsea1<-pal_simpsons()(16)

num <- 1
gseaplot2(gsea.re1,geneSetID = rownames(g1)[num],
          title = rownames(g1)[num],
          color = col_gsea1[num],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)


dev.off()



########3 Jun 1 ############
res1 <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/DEGs.Δ133\ tumor_vs_mCherry2\ tumor.csv')
res1$delabel1 <- NA
res1$delabel1[which(-log10(res1$padj)>=1.3)] <- res1$gene[which(-log10(res1$padj)>=1.3)]
res1$delabel1[which(abs(res1$log2FoldChange)<=1)] <- NA


cmp <- 'Δ133 tumor'
ref <- 'mCherry2 tumor'
ggplot(data = res1, 
       aes(x = log2FoldChange,
           y = -log10(res1$padj), 
           colour= factor(diffexpressed),
           label = delabel1)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(point.padding = 0.2, 
                  size=3.4,
                  nudge_x = .15,
                  nudge_y = .5,
                  segment.linetype = 2,
                  arrow = arrow(length = unit(0.015, "npc"))) +
 geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="grey",lwd=0.8) +
  labs(x="log2FoldChange)",
       y="-log10 (padj)",
       title=paste(cmp,' vs ',ref))  +
  theme_bw()+
  labs(color='padj<=0.05
 log2FoldChange|>=1') +
  # you should specify the color also
  #scale_color_manual(labels = c("UP 587", "DOWN 815",'NO 57333')
  #                   ,values = c("#F8766D","#00BA38","#619CFF")) +
  guides(colour = guide_legend())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.key = element_rect(fill = "snow2"))+
  ggeasy::easy_center_title()

## only for product
res <- res1
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
res$delabel <- NA
#res$gene <- genecount$Gene.name[which(row.names(res) == row.names(genecount))]
res$delabel[res$diffexpressed != "NO"] <-res$gene[which(res$diffexpressed != "NO")]
res1 <- as.data.frame(res)
ggplot(data=res1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 3, box.padding = 0.01,max.overlaps = 10)+
  geom_vline(xintercept=c(-0.01,0.01),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="grey",lwd=0.8) +
  labs(x="log2FoldChange)",
       y="-log10 (padj)",
       title=paste(cmp,' vs ',ref))  +
  theme_bw()+
  labs(color='padj<=0.05
 log2FoldChange|>=0') +
  # you should specify the color also
  #scale_color_manual(labels = c("UP 587", "DOWN 815",'NO 57333')
  #                   ,values = c("#F8766D","#00BA38","#619CFF")) +
  guides(colour = guide_legend())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.key = element_rect(fill = "snow2"))+
  ggeasy::easy_center_title()


##### heatmap of geneset ##### 

genecount <-  read.csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/raw_counts.csv')
genecount <- genecount %>% column_to_rownames(var='ID')

genecount <-  read.csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/raw_counts.csv')
genecount <- genecount %>% column_to_rownames(var='ID')

#sample list 1 tumor/product ## 1801 = Δ133p53-P2A-CTL119 1924 = mCherry2-P2A-CTL119
#sample.list1 <-data.frame('sample'=colnames(genecount)[1:13],'groups'=c(rep('tumor',6),rep('product',6)),row.names=NULL)
sample.list2 <- data.frame('sample'=colnames(genecount)[1:6],'groups'=rep(c('mCherry2 tumor','Δ133 tumor'),3))
sample.list3 <- data.frame('sample'=colnames(genecount)[7:12],'groups'=rep(c('Δ133 product','mCherry2 product'),3))
sample.list <- rbind(sample.list2,sample.list3)
sample.list$groups <- as.factor(sample.list$groups)

tumor <- c('Δ133 tumor','mCherry2 tumor')

sample.list.g <- sample.list[sample.list$groups %in% tumor,]

c <- as.matrix(genecount[,sample.list$sample[sample.list$groups %in% tumor]])
#row.names(c) <- genecount$Gene.name[which(row.names(c) == row.names(genecount))]
dds <- DESeqDataSetFromMatrix(c, sample.list.g, ~ groups)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

geneset <- read_excel('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/geneSets_d133.xlsx')
genelist <- c('term'=c(),'gene'=c())
for(j in colnames(geneset)){
  a <- rep(j,nrow(geneset))
  b <- geneset[,j]
  c <- cbind(a,b)
  colnames(c) <- c('term','gene')
  genelist <- rbind(c,genelist)
}

genelist <- na.omit(genelist)

### Building matrix 
#c <- as.matrix(genecount[,1:12])
#row.names(c) <- genecount$Gene.name[which(row.names(c) == row.names(genecount))]
#dds1 <- DESeqDataSetFromMatrix(c, sample.list, ~ groups)

#dds1 <- DESeq(dds1)

deseq2VST2  <- vst(dds)
deseq2VST2 <- assay(deseq2VST2)
deseq2VST2 <- as.data.frame(deseq2VST2)
genecount1 <- genecount[keep,]
deseq2VST2$Gene <- genecount$Gene.name[which(row.names(deseq2VST2) == row.names(genecount1))]
head(deseq2VST2)

genes <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/DEGs.Δ133\ tumor_vs_mCherry2\ tumor.csv')
genes <- genes[order(genes$padj),]

genelist1 <- genes$...1[which(genes$padj<=0.05)]
genelist2 <- genelist1[genelist1 %in% genelist$gene ]
genelist.name <- genes$gene[genes$...1 %in% genelist2]
genelist.f <- genelist.name[0:100]



deseq2VST1 <- deseq2VST2[deseq2VST2$Gene %in% genelist.name,]
row.names(deseq2VST1) <- deseq2VST1$Gene
library(reshape2)
anno <- sample.list[1:6,]
anno$test <- as.character(sample.list$groups[1:6])
anno$groups <- c()
anno$conditions <- c()
for (i in 1:nrow(anno)){
  anno$groups[i] <- as.list(strsplit(anno$test,' '))[[i]][1]
  anno$conditions[i] <- as.list(strsplit(anno$test,' '))[[i]][2]
}
anno1 <-data.frame(anno[,3],row.names=anno$sample)
colnames(anno1) <- 'group'

deseq2VST_plot <- as.data.frame(lapply(deseq2VST1,as.numeric))
row.names(deseq2VST_plot) <- deseq2VST1$Gene
pheatmap(t(deseq2VST_plot[,-7]),fontsize=6.5,main='Δ133 tumor vs mCherry2 tumor',border_color=NA
         ,clustering_distance_cols = 'euclidean'
         ,clustering_distance_rows = 'euclidean',annotation_row =anno1)



####### GOenrighc #######
path <- '/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/'
f.list <- list.files(path,pattern='.csv')
setwd(path)
deg <- read_csv(f.list[5])

## All gene ranking
deg1 <- deg[deg$padj<=0.05,]
alldiff <- deg1[order(deg1$log2FoldChange,decreasing = T),]
id1 <- alldiff$log2FoldChange
names(id1) <- alldiff$...1
id1 <- na.omit(id1)
id1 <- sort(id1,decreasing=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
ggo <- enrichGO(gene     = names(id1),
                keyType = "ENSEMBL",
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               readable = TRUE,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)

head(ggo)

setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/July_24')
write.csv(ggo,file=str_replace(f.list[5],'.csv','BP.gsea.csv'))

dev.new()
pdf(str_replace(f.list[5],'.csv','BP.gsea.pdf'),height=10)

dotplot(ggo, showCategory=30) +
  scale_size_area(max_size = 5)
dev.off()


###### New Volcano Plot ####


res1 <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/DEGs.Δ133\ product_vs_mCherry2\ product.csv')
res1$delabel1 <- NA
res1$delabel1[which(-log10(res1$padj)>=1.3)] <- res1$gene[which(-log10(res1$padj)>=1.3)]
#res1$delabel1[which(abs(res1$log2FoldChange)<=1)] <- NA
res1$diffexpressed1 <- "NO"
res1 <- res1[is.na(res1$padj)==FALSE,]

for (i in 1:nrow(res1)){
  if (res1$padj[i]<=0.05 & res1$log2FoldChange[i]<0){
    print(i)
    res1$diffexpressed1[i] <-  'DOWN'
  }else if(res1$padj[i]<=0.05 & res1$log2FoldChange[i]>0){
    print(i)
    res1$diffexpressed1[i] <-  'UP'
  }
} 



cmp <- 'Δ133 product'
ref <- 'mCherry2 product'
ggplot(data = res1, 
       aes(x = log2FoldChange,
           y = -log10(res1$padj), 
           colour= factor(diffexpressed1),
           label = delabel1)) +
  geom_point(size=1.5) + 
  theme_minimal() +
  geom_text_repel(point.padding = 0.2, 
                  size=3.0,
                  nudge_x = .15,
                  nudge_y = .5,
                  segment.linetype = 2,
                  arrow = arrow(length = unit(0.015, "npc")),
                  max.overlaps =10) +
  #geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="grey",lwd=0.8) +
  labs(x="log2FoldChange)",
       y="-log10 (padj)",
       title=paste(cmp,' vs ',ref))  +
  theme_bw()+
  labs(color='padj<=0.05') +
  # you should specify the color also
  scale_color_manual(values = c("red","grey27","blue")) +
  guides(colour = guide_legend())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.key = element_rect(fill = "snow2"))+
  ggeasy::easy_center_title()

######## Product ######## 
## only for product
res <- res1
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
res$delabel <- NA
#res$gene <- genecount$Gene.name[which(row.names(res) == row.names(genecount))]
res$delabel[res$diffexpressed != "NO"] <-res$gene[which(res$diffexpressed != "NO")]
res1 <- as.data.frame(res)
ggplot(data=res1, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 3, box.padding = 0.01,max.overlaps = 10)+
  geom_hline(yintercept = 1.301,lty=4,col="grey",lwd=0.8) +
  labs(x="log2FoldChange)",
       y="-log10 (padj)",
       title=paste(cmp,' vs ',ref))  +
  theme_bw()+
  labs(color='padj<=0.05
 log2FoldChange|>=0') +
  scale_color_manual(values = c("red","grey27","blue")) +
  # you should specify the color also
  #scale_color_manual(labels = c("UP 587", "DOWN 815",'NO 57333')
  #                   ,values = c("#F8766D","#00BA38","#619CFF")) +
  guides(colour = guide_legend())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.key = element_rect(fill = "snow2"))+
  ggeasy::easy_center_title()


#######GO term#####
path <- '/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/'
f.list <- list.files(path,pattern='.csv')
setwd(path)
deg.o <- read_csv(f.list[1])

## All gene ranking
deg <- deg.o[is.na(deg.o$log2FoldChange)==FALSE,]
n <- nrow(deg)
sdf <- sd(deg$log2FoldChange)
margin <- qt(0.975,df=n-1)*sdf/sqrt(n)
low <- mean(deg$log2FoldChange)-margin
high <- mean(deg$log2FoldChange)+margin
deg <- deg[deg$log2FoldChange<=high,]
deg <- deg[deg$log2FoldChange>=low,]

deg1 <- deg[deg$padj<=0.05,]
deg1 <- deg1[deg1$log2FoldChange>=0,]
alldiff <- deg1[order(deg1$log2FoldChange,decreasing = T),]
id1 <- alldiff$log2FoldChange
names(id1) <- alldiff$...1
id1 <- na.omit(id1)
id1 <- sort(id1,decreasing=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
ggo <- enrichGO(gene     = names(id1),
                keyType = "ENSEMBL",
                OrgDb    = org.Hs.eg.db,
                ont      = "BP",
                readable = TRUE,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

head(ggo)

setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/Aug_01')
write.csv(ggo,file=str_replace(f.list[1],'.csv','.UP.BP.gsea.csv'))

dev.new()
pdf(str_replace(f.list[1],'.csv','.UP.BP.gsea.pdf'),height=10)

dotplot(ggo, showCategory=30) +
  scale_size_area(max_size = 5)
dev.off()


path <- '/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/May_24/DEG/'
f.list <- list.files(path,pattern='.csv')
setwd(path)
deg <- read_csv(f.list[4])

## All gene ranking
deg1 <- deg[deg$padj<=0.05,]
deg1 <- deg1[deg1$log2FoldChange>=0,]
alldiff <- deg1[order(deg1$log2FoldChange,decreasing = T),]
id1 <- alldiff$log2FoldChange
names(id1) <- alldiff$...1
id1 <- na.omit(id1)
id1 <- sort(id1,decreasing=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
ggo <- enrichGO(gene     = names(id1),
                keyType = "ENSEMBL",
                OrgDb    = org.Hs.eg.db,
                ont      = "BP",
                readable = TRUE,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

head(ggo)

setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/Aug_01')
write.csv(ggo,file=str_replace(f.list[4],'.csv','.DOWN.BP.gsea.csv'))

dev.new()
pdf(str_replace(f.list[4],'.csv','.DOWN.BP.gsea.pdf'),height=10)

dotplot(ggo, showCategory=30) +
  scale_size_area(max_size = 5)
dev.off()