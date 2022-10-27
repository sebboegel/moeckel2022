library("tximport")
library("DESeq2")
library("ggplot2")
library(EnhancedVolcano)
library("pheatmap")
library(reshape2)
library("tidyr")
library(tidyverse)
dir="./samples_2022"
samples <- read.table(file.path(dir,"metadata_complete.txt"), header=TRUE)
rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")

##################################
#load biomart


#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "mmusculus_gene_ensembl",
#                         host = "ensembl.org")
#ttg <- biomaRt::getBM(
#  attributes = c("ensembl_transcript_id", "transcript_version",
#                 "ensembl_gene_id", "external_gene_name", "description",
#                 "transcript_biotype"),
#  mart = mart)

###################################
ttg=readRDS("./biomart.rds")
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))

txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)

rownames(samples) <- samples$alias
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design = ~condition)
dds=DESeq(ddsTxi)

#Figure 2 ------ PCA
vsd <- vst(dds, blind = FALSE)
nudge <- position_nudge(y = 1)

pcaData <- plotPCA(vsd, intgroup=c("condition", "strain"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=factor(strain,levels = c("BAFF-ko", "BAFF-wt" , "BAFF-R-ko",   "BAFF-R-wt"),labels=c("BAFF ko", "BAFF wt","BAFF-R ko", "BAFF-R wt")))) +
  geom_text(aes(label = name), size=3,position = nudge) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  labs(shape='strain')

##proceed with adapted metadata file as L6_wt will not be considered in the further analyses
samples <- read.table(file.path(dir,"metadata.txt"), header=TRUE)
rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")

ttg=readRDS("/Users/sebastian/rheuma_seq/rna-seq/mouse/biomart.rds")
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))

txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)

rownames(samples) <- samples$alias
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design = ~condition)
dds=DESeq(ddsTxi)
res <- results(dds)
res["Havcr1",]$padj
res["Lcn2",]$padj
res["Lyz2",]$padj
res["Fn1",]$padj
res["Kl",]$padj
res["Tnfsf13b",]$padj
res["Tnfrsf13c",]$padj
res["Tnfrsf17",]$padj
res["Tnfrsf13b",]$padj


#Figure 3 ------ Volcano Plot
EnhancedVolcano(res,lab = rownames(res), x = 'log2FoldChange',y = 'pvalue',xlim = c(-10, 15))

#Figure 6&7 -----
abundance_table=txi$abundance
colnames(abundance_table)=samples$run_accession
all_mouse <- cbind(rownames(abundance_table), data.frame(abundance_table, row.names=NULL))
colnames(all_mouse)[1]="gene"
m <- melt(all_mouse)
m <- na.omit(m)
colnames(samples)=c("variable","strain","condition")
test=merge(m,samples,by="variable")

#Figure 6----------
gene_set=c("Havcr1","Lcn2","Lyz2","Fn1","Kl")
gene_set=c("Fn1","Kl")
gene_set=c("Cd44","Il1rn")
test$facet = factor(test$gene, levels = gene_set)

test %>% 
  filter(gene %in% gene_set) %>%
  ggplot(aes(x=condition,y=value,fill=factor(strain,levels = c("BAFF-ko", "BAFF-wt" , "BAFF-R-ko",   "BAFF-R-wt"),labels=c("BAFF ko", "BAFF wt","BAFF-R ko", "BAFF-R wt")))) +
  geom_boxplot() + 
  labs(fill = "strain",y="TPM") +  
  geom_point(position=position_jitterdodge()) +
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=1,strip.position="right")

#Figure 7 --------
gene_set = c("Tnfsf13b","Tnfrsf13c","Tnfrsf17","Tnfrsf13b")

test$facet = factor(test$gene, levels = gene_set)
test %>% 
  filter(gene %in% gene_set) %>%
  ggplot(aes(x=condition,y=value,fill=factor(strain,levels = c("BAFF-ko", "BAFF-wt" , "BAFF-R-ko",   "BAFF-R-wt"),labels=c("BAFF ko", "BAFF wt","BAFF-R ko", "BAFF-R wt")))) +
  geom_boxplot() + 
  labs(fill = "strain",y="TPM") +  
  geom_point(position=position_jitterdodge()) +
  #theme_bw(base_size = 16) + factor(gene,levels=gene_set,ncol=1,strip.position="right")
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=1,strip.position="right")


#BAFF only ----------------------------

samples <- read.table(file.path(dir,"metadata_baff"), header=TRUE)
rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")

ttg=readRDS("/Users/sebastian/rheuma_seq/rna-seq/mouse/biomart.rds")
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))

txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)

rownames(samples) <- samples$alias
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design = ~condition)
dds=DESeq(ddsTxi)
res <- results(dds)
vsd <- vst(dds, blind = FALSE)
#Figure 4 ------
EnhancedVolcano(res,lab = rownames(res), x = 'log2FoldChange',y = 'pvalue',xlim = c(-10, 15))

#Figure 5 -------
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("condition")])
rownames(anno) <- colnames(mat)
colnames(anno) <- "condition"
pheatmap(mat, annotation_col = anno)


#BAFF-R only ----------------------------

samples <- read.table(file.path(dir,"metadata_baff-r"), header=TRUE)
rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")

ttg=readRDS("/Users/sebastian/rheuma_seq/rna-seq/mouse/biomart.rds")
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))

txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)

rownames(samples) <- samples$alias
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design = ~condition)
dds=DESeq(ddsTxi)
res <- results(dds)
vsd <- vst(dds, blind = FALSE)
#Figure 4 ------
EnhancedVolcano(res,lab = rownames(res), x = 'log2FoldChange',y = 'pvalue',xlim = c(-10, 15))

#Figure 5 -------
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("condition")])
rownames(anno) <- colnames(mat)
colnames(anno) <- "condition"
pheatmap(mat, annotation_col = anno)





