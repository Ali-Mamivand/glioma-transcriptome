##---- loading packages --------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(PoiClaClu)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReportingTools)
library(tibble)
library(dplyr)
library(tidyr)
library(fgsea)

##----loading data---------------
# Download GBM and LGG samples from TCGA
set.seed(1000)
lgg.samples <- matchedMetExp("TCGA-LGG", n = 20)
gbm.samples <- matchedMetExp("TCGA-GBM", n = 20)
samples <- c(lgg.samples,gbm.samples)
query_exp <- GDCquery(project= c("TCGA-LGG","TCGA-GBM"), 
                      data.category = "Transcriptome Profiling", 
                      workflow.type = 'HTSeq - Counts', 
                      barcode = samples)
GDCdownload(query_exp)
counts.table=GDCprepare(query_exp)
Count.Table=assay(counts.table)
table(colData(counts.table)$definition)
#remove normal samples (if any):
counts.table=counts.table[,(colData(counts.table)$definition != "Solid Tissue Normal")]
clinical_information = as.data.frame(colData(counts.table))
##------Create DEseq2 object----------
(ddsMat=DESeqDataSetFromMatrix(countData = Count.Table,
                               colData = clinical_information,
                               design = ~ project_id))
##------Quality Control-----------
dim(ddsMat)
#remove reads with 0 counts across all samples:
dds <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]
dim(dds)
dds_vst <- vst(dds)
head(assay(dds_vst), 3)

sampleDists <- dist( t( assay(dds_vst) ) ) 
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- dds_vst$project_id
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- dds_vst$project_id
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
#PCA 
plotPCA(dds_vst, intgroup = c("project_id"))
#Classical multidimensional scaling (MDS)
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(dds_vst)))
ggplot(mds, aes(X1,X2,color=project_id)) + geom_point(size=3)
##---- (DEA)------------
#DEA
dds <- DESeq(dds)
(res <- results(dds))
mcols(res, use.names=TRUE)
summary(res)
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.05)
# Visualization of gene expression
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
data <- plotCounts(dds, gene=topGene,
                   intgroup="project_id", returnData=TRUE)
ggplot(data, aes(x=project_id, y=count)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0),
             size=3)+ggtitle(topGene)
# Alternative way :
# plotCounts(dds, gene=topGene, intgroup=c("project_id"))


plotMA(res, ylim=c(-5,5))
plotMA(resLFC1, ylim=c(-5,5))
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="black", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="black")
})
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")

topVarGenes <- head(order(rowVars(assay(dds_vst)),decreasing=TRUE),20)
mat <- assay(dds_vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_vst)["project_id"])
pheatmap(mat, annotation_col=df)
##---Annotation-------------------
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
##---- report results--------
res <- res[order(res$padj),]
head(res)
res <- as.data.frame(res)
# Save DEA results in csv:
write.csv(res,"DEA.csv")
#Export result as an HTML file:
htmlRep <- HTMLReport(shortName="report", title="My report")
publish(res, htmlRep)
url <- finish(htmlRep)
browseURL(url)
##---- Pathway Enrichment -------
# Load the pathway (gene set) into a named list.
# Download mysigdb from https://www.gsea-msigdb.org/gsea/msigdb/ and put it 
# in the same directory with this code.
# Here I used h.all.v7.2.symbols.gmt:
pathway_Dir = 'F:/RNA-Seq/glioma/LGG_Normal/' #set it to the your directory that contains gmt file
pathway <- gmtPathways(paste0(pathway_Dir,"h.all.v7.2.symbols.gmt"))
# show few lines from the pathways file
head(pathway)

res_filtered <- res %>% 
  select(symbol,log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(log2FoldChange))
res_filtered
# creating  a named vector [ranked genes]
ranks <- res_filtered$stat
names(ranks) <- res_filtered$symbol

#Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathway, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

# To see what genes are in each of these pathways:
gene.in.pathway <- pathway %>% 
  enframe("pathway", "symbol") %>% 
  unnest(cols = c(symbol)) %>% 
  inner_join(res_filtered, by="symbol")
# Plot pathways based on the normalized enrichment scores(NES). 
# Color the bar indicating whether or not the pathway was significant(based on FDR):
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
fgseaResTidy=fgseaResTidy[order(fgseaResTidy$padj),]
fgseaRes=fgseaRes[order(fgseaRes$padj),]

pdf('Hallmark pathways enrichment.pdf',width = 11,height = 8.5)
ggplot(fgseaResTidy[fgseaResTidy$padj<0.0005,], aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")
dev.off()

# Enrichment plot for specific target gene with lowest FDR:
top_pathway <- fgseaRes$pathway[which.min(fgseaRes$padj)]

pdf('Enrichment_plot.pdf',width = 8,height = 6)
plotEnrichment(pathway = pathway[[top_pathway]], ranks)+
  ggtitle(top_pathway)
dev.off()

# Enrichment plot for multiple pathway:
pdf('top_ten_patway.pdf',width = 8.5,height = 11)
plotGseaTable(pathway[fgseaRes$pathway[1:10]], ranks, fgseaRes, 
              gseaParam=0.5)
dev.off()



