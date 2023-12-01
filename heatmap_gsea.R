
## Pathway Analysis 

setwd('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report')
deg <- read_csv('DEGs.mCherry2 tumor_vs_Δ133 tumor.csv')

## All gene ranking
alldiff <- deg[order(deg$log2FoldChange,decreasing = T),]
id1 <- alldiff$log2FoldChange
names(id1) <- alldiff$gene
id1 <- na.omit(id1)
id1 <- sort(id1,decreasing=TRUE)

## Load gmt file
gmtfile <- "/Users/linhuichen/Downloads/h.all.v7.4.symbols.gmt"
hallmark <- read.gmt(gmtfile)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term)  %>% lapply( "[[", 2)


gsea.re2 <- fgsea(pathways =  hallmark.list,
                  stats = id,
                  minSize=1,
                  maxSize=10000,
                  eps=0
)
colnames(gsea.re2)
#
g2 <- gsea.re2[gsea.re2$padj<0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]

dev.new()
pdf('All_gene_gsea.pdf')
plotGseaTable(hallmark.list[g2$pathway],
              id, 
              gsea.re2,gseaParam = 0.5,
              colwidths = c(0.5,0.2,0.1,0.1,0.1)
)


names(hallmark.list)
se_hall<-c(head(g2$pathway,2),tail(g2$pathway,2))

sig1<-g2[g2$pathway%in%se_hall,]
hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(9)
ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#
            ticksize=0.1,
            #colour='grey',
            linecolor=pal_line,
            linesize=2,lty=1,
            #plot number for each row
            nrow = 2,
            ncol = 2
  ) + 
  theme_bw()+
  theme(strip.text = element_text(size = 7,face = 'italic'),#
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=8,angle = -30,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=8,angle = 0,face = 'plain',vjust = 0.5))
dev.off()

deg_down <- deg[deg$log2FoldChange<0,]
down_list <- na.omit(deg_down$gene)
gsea.re2$down_gene <- c()
for (i in 1:nrow(gsea.re2)){
  print(gsea.re2$leadingEdge[[i]] [which(gsea.re2$leadingEdge[[i]] %in% down_list)])
  gsea.re2$down_gene[[i]] <- gsea.re2$leadingEdge[[i]] [which(gsea.re2$leadingEdge[[i]] %in% down_list)]
}
fwrite(gsea.re2 , file="all_gsea_mCherry2_tumor_vs_Δ133_tumor.tsv", sep="\t", sep2=c("", " ", ""))






### HEATMAP
genecount <-  read.csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/data/raw_counts.csv')
genecount <- genecount %>% column_to_rownames(var='ID')

#sample list 1 tumor/product ## 1801 = Δ133p53-P2A-CTL119 1924 = mCherry2-P2A-CTL119
#sample.list1 <-data.frame('sample'=colnames(genecount)[1:13],'groups'=c(rep('tumor',6),rep('product',6)),row.names=NULL)
sample.list2 <- data.frame('sample'=colnames(genecount)[1:6],'groups'=rep(c('mCherry2_tumor','Δ133_tumor'),3))
sample.list3 <- data.frame('sample'=colnames(genecount)[7:12],'groups'=rep(c('Δ133_product','mCherry2_product'),3))
sample.list <- rbind(sample.list2,sample.list3)
sample.list$groups <- as.factor(sample.list$groups)


## Building matrix 
c <- as.matrix(genecount[,1:12])
#row.names(c) <- genecount$Gene.name[which(row.names(c) == row.names(genecount))]
dds1 <- DESeqDataSetFromMatrix(c, sample.list, ~ groups)

dds1 <- DESeq(dds1)

deseq2VST2  <- vst(dds1)
deseq2VST2 <- assay(deseq2VST2)
deseq2VST2 <- as.data.frame(deseq2VST2)
deseq2VST2$Gene <- count$Gene.name[which(row.names(c) == row.names(count))]
head(deseq2VST2)


#mcherry vs delta tumor
genes <- read_csv('/Users/linhuichen/Desktop/projects/Regina_Young_RNASeq_Chris_2023_01/report/DEGs.Δ133\ tumor_vs_Δ133\ product.csv')

genelist <- genes$gene[which(genes$padj<=0.05)]
genelist <-  genes$gene[order(genes$padj)[0:100]]

deseq2VST1 <- deseq2VST2[deseq2VST2$Gene %in% genelist,]
row.names(deseq2VST1) <- deseq2VST1$Gene
library(reshape2)
anno <- sample.list
anno$test <- as.character(sample.list$groups)
anno$groups <- c()
anno$conditions <- c()
for (i in 1:nrow(anno)){
  anno$groups[i] <- as.list(strsplit(anno$test,'_'))[[i]][1]
  anno$conditions[i] <- as.list(strsplit(anno$test,'_'))[[i]][2]
}
anno1 <-data.frame(anno[,3:4],row.names=anno$sample)
pdf('delta133_tumor_vs_Δ_product_top_100_heatmap.pdf',width=15)
deseq2VST_plot <- as.data.frame(lapply(deseq2VST1,as.numeric))
row.names(deseq2VST_plot) <- deseq2VST1$Gene
dev.new()
pheatmap(t(deseq2VST_plot[,-13]),fontsize=6.5,main='Δ133 tumor vs Δ133 product'
         ,clustering_distance_cols = 'euclidean'
         ,clustering_distance_rows = 'euclidean',annotation_row =anno1)
dev.off()




##Down Regulate  gene ranking No Enrichment
deg_down <- deg[deg$log2FoldChange<0,]
alldiff <-deg_down[order(deg_down$log2FoldChange,decreasing = F),]
id <- alldiff$log2FoldChange
names(id) <- alldiff$gene
id <- na.omit(id)
id <- sort(id,decreasing=TRUE)

## Load gmt file
gmtfile <- "/Users/linhuichen/Downloads/h.all.v7.4.symbols.gmt"
hallmark <- read.gmt(gmtfile)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T,nPermSimple = 10000)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term)  %>% lapply( "[[", 2)


gsea.re2 <- fgsea(pathways =  hallmark.list,
                  stats = id,
                  minSize=1,
                  maxSize=10000,
                  eps=0
)
colnames(gsea.re2)
#
g2 <- gsea.re2[gsea.re2$padj<0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]

dev.new()
pdf('Down_gene_gsea.pdf')
plotGseaTable(hallmark.list[g2$pathway],
              id, 
              gsea.re2,gseaParam = 0.5,
              colwidths = c(0.5,0.2,0.1,0.1,0.1)
)

fwrite(gsea.re2 , file="Down_gene_gsea.tsv", sep="\t", sep2=c("", " ", ""))


names(hallmark.list)
se_hall<-c(head(g2$pathway,2),tail(g2$pathway,2))

sig1<-g2[g2$pathway%in%se_hall,]
hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(9)
ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#
            ticksize=0.1,
            #colour='grey',
            linecolor=pal_line,
            linesize=2,lty=1,
            #plot number for each row
            nrow = 2,
            ncol = 2
  ) + 
  theme_bw()+
  theme(strip.text = element_text(size = 7,face = 'italic'),#
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=8,angle = -30,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=8,angle = 0,face = 'plain',vjust = 0.5))
dev.off()


############