#### Packages ####
library("DESeq2")
library("affycoretools")
library("tidyverse")
library("pheatmap")
library("dplyr")
library("biomaRt")
library("AnnotationDbi")
library('org.Hs.eg.db')
library(readxl)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(keras)
library(tfdatasets)
library(tidyverse)
library(rsample)
library(rcompanion)
library(Rmisc)
library(nnet)
library(ccoptimalmatch)
library(MatchIt)
library(epitools)
library(optmatch)
library(gmodels)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(nnet)
library(performance)
library(abd)
library(EpiStats)
library(knitr)
library(questionr)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(tidyverse)
library(finalfit)
library(purrr)
library(clipr)
library(Seurat)
library(dplyr)
library(Matrix)
library(readr)
library(monocle3)
library(SeuratWrappers)
library(SeuratDisk)
library(celldex)
library(SingleR)
library(edgeR)
library(ComplexHeatmap)
library(ggsci)
library(clusterProfiler)
library(circlize)
library(GetoptLong)
library(EnhancedVolcano)

#### Load data ####
tpm <- read.csv("/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse SEQ/mDCM_Run1.csv")
tpm2 <- read.csv("/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse SEQ/mDCM_Run2.csv")

tpm <- merge(tpm, tpm2, by = "gene_id")
tpm <- dplyr::select(tpm, -c("Gene_name.y"))
tpm <- dplyr::select(tpm, -c("DECO1", "DECO2", "DECO3")) #"M1_12W", "M2_12W", "M3_12W", "F1_12W", "F2_12W", "F3_12W"))

write.csv(tpm, "/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse SEQ/mDCM_combined.csv")

tpm <- read.csv("/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse SEQ/Share_RNAseq/mDCM_combined.csv")
md <- read.csv("/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse SEQ/Share_RNAseq/SampleList_MouseSeq_DCM.csv")

metadata <- data.frame(md[,1:4])

#### Gene lists ####
matrisome <- read_xlsx("/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse Seq/matrisome_mouse.xlsx")
matrisome <- lapply(matrisome, as.character)
matrisome <- as.data.frame(matrisome)
matrisome_glyc <- matrisome$ECM_Glycoproteins
matrisome_Coll <- matrisome$Collagens
matrisome_Prot <- matrisome$Proteoglycans
matrisome_Aff <- matrisome$ECM_Affiliated
matrisome_Reg <- matrisome$ECM_Regulators
matrisome_all <- matrisome$All

colnames(matrisome) <- c("GP", "Collagens", "PG", "Affiliated", "Regulators", "Secreted", "All")

matrisome_long <- matrisome[,1:6]
matrisome_long <- matrisome_long %>%
  mutate(Row = row_number()) %>% 
  pivot_longer(
    cols = -Row,
    names_to = "Class",   # New column to store previous column names
    values_to = "Gene"    # New column to store values
  )

bm <- read_xlsx("/Users/noah/Desktop/ECM Paper Share/Figure 1/Mouse SEQ/BM_Genes_mouse.xlsx")
bm_genes <- bm$Symbol


#### Subset to protein-coding genes ####
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
txtypes <- genes(edb, columns=c("gene_name", "gene_biotype", "tx_biotype", "tx_id"))
protGenes <- genes(edb, filter=GeneBiotypeFilter("protein_coding"))
PG <- protGenes$gene_name
PG2 <- protGenes$gene_id

id_name <- data.frame(PG, PG2)

tpm <- tpm %>% dplyr::filter(gene_id %in% PG2)

#### DEG ####
tpm <- tpm[,-1]
genelist <- tpm[,1:2]
cts <- tpm[,3:20]
rownames(cts) <- tpm[,1]

group_coarse <- c("SHAM","SHAM","12_DCM","4_DCM","4_DCM","SHAM","12_DCM","12_DCM","4_DCM",
                  "SHAM","SHAM","12_DCM","4_DCM","4_DCM","SHAM","12_DCM","12_DCM","4_DCM")

y <- DGEList(counts = cts, group = group_coarse, genes = genelist)
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = y$samples,
  design = ~group)

dds <- DESeq(dds)
normdds <- vst(dds)
res <- results(dds, contrast=c("group","4_DCM","SHAM"))
res <- data.frame(res)
res$gene_id <- rownames(res)
res$padj <- p.adjust(res$pvalue, method = "BH")
res <- merge(res, genelist, by = "gene_id")


normalized_counts <- data.frame(assay(normdds))
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- data.frame(normalized_counts)
normalized_counts$gene_id <- rownames(normalized_counts)
colnames(genelist) <- c("gene_id", "gene_name")
normalized_counts <- merge(normalized_counts, genelist, by = "gene_id")

up <- dplyr::filter(res, padj <0.05 & log2FoldChange <0)
down <- dplyr::filter(res, padj <0.05 & log2FoldChange >0)

#### BM heatmap ####

library(GetoptLong)
x = normalized_counts %>% dplyr::filter(gene_name %in% matrisome_all)
tpm <- x[!duplicated(x$gene_name),]
rownames(tpm) <- tpm$gene_name
tpm <- na.omit(tpm)
x <- tpm[,2:19]
rownames(x) <- tpm$gene_name
x <- x[rowSums(x[])>0,]
x <- as.matrix(x)
x <- t(scale(t(x), scale = TRUE, center = TRUE))


order <- c("M1_4W", "M2_4W", "M3_4W", "F4_4W", "F5_4W", "F6_4W",
           "M1_12W", "M2_12W", "M3_12W", "F1_12W", "F2_12W", "F3_12W",
           "M1_SHAM", "M2_SHAM", "M3_SHAM", "F1_SHAM", "F2_SHAM", "F3_SHAM")

group <- factor(group_coarse, levels = c("4_DCM", "12_DCM", "SHAM"))

row_split <- matrisome_long$Class[match(rownames(x), matrisome_long$Gene)]

Heatmap(x,
        name = "Expression",
        column_order = order,
        show_row_names = FALSE, 
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        column_split = group,
        row_split = row_split,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        column_title = c("DCM-4w", "DCM-12w", "Sham"))


#### Barplots ####

## Filter and prep data
as_plot <- normalized_counts %>% dplyr::filter(gene_name %in% fibr | gene_name %in% ra_genes | gene_name %in% ba_genes)
as_plot <- data.frame(as_plot)
as_plot <- t(as_plot)
as_plot_b <- data.frame(as_plot)
as_plot_b$animals <- rownames(as_plot_b)
write_csv(as_plot_b, "toplotgp5.csv")

as_plot <- normalized_counts %>% dplyr::filter(gene_name %in% bm_genes)
as_plot <- data.frame(as_plot)
as_plot <- t(as_plot)
as_plot_b <- data.frame(as_plot)
as_plot_b$animals <- rownames(as_plot_b)
write_csv(as_plot_b, "toplotgp6.csv")

#### PCA ####
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("group"))


x <- normalized_counts[,2:19]
x<- x[rowSums(x[])>300,]
xa <- as.matrix(x)
xb <- t(xa)
xb <- xb[, colSums(xb != 0) > 0]
pca <- t(xa)
pca <- data.frame(pca)
pca$group_coarse <- c("SHAM","SHAM","12_DCM","4_DCM","4_DCM","SHAM","12_DCM","12_DCM","4_DCM",
                      "SHAM","SHAM","12_DCM","4_DCM","4_DCM","SHAM","12_DCM","12_DCM","4_DCM")
pcax <- subset(pca, select = -c(group_coarse))
pca_res <- prcomp(xb, scale. = TRUE)
library(ggfortify)
autoplot(pca_res, data = pca, colour = c("group_coarse"), size = 3) + theme_bw() + 
  scale_color_npg() + scale_fill_npg() + guides(size = "none")



#### Volcano plots ####
res <- res %>% 
  mutate(Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                                log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))
top <- 12
top_genes <- bind_rows(
  res %>% 
    dplyr::filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  res %>% 
    dplyr::filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top))
res <- na.omit(res)

ggplot(res, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = Expression), size = 2/5) + 
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"PValue")) +
  theme_pubr() + theme(legend.position="none") +
  geom_label_repel(data = top_genes, mapping = aes(log2FoldChange, -log(padj,10), 
                                                   label = Gene_name.x),size = 6,max.overlaps = 30)

#### GO ####
library(org.Mm.eg.db)
gene_list <- dplyr::filter(res, padj <0.05, log2FoldChange <0)
entrez_ids <- bitr(gene_list$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = entrez_ids$ENTREZID, 
                OrgDb = org.Mm.eg.db, 
                keyType = "ENTREZID", 
                ont = "CC", 
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05, 
                qvalueCutoff = 0.05, 
                readable = TRUE)

ego_simplified <- simplify(ego, cutoff = 0.8, by = "p.adjust", select_fun = min)

barplot(ego_simplified, showCategory = 5) + theme(legend.position="none") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
