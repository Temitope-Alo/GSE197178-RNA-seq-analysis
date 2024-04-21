# Working directory
setwd("C:/Users/temit/Desktop/R")

# load libraries
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(magick)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 1: Processing count data ----------------

counts_data <- read.csv("GSE197178_mRNA_gene.FPKM.csv")
Counts <- counts_data[, -2:-3]

rownames(Counts) <- Counts$gene_id
Counts <- Counts[, -which(names(Counts) == "gene_id")]

Counts <- round(Counts)
Counts <- Counts[which(rowSums(Counts) > 1),]

# 2: Setting condition

condition <- factor(c("Negative","Negative","Negative", "Negative", 
                      "Positive","Positive","Positive", "Positive"))

coldata <- data.frame(row.names = colnames(Counts), condition)
coldata

# 3: Constructing DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata, intgroup = "condition")
plotDispEsts(dds)

# 4: Pre-processing DGE Results I ------------------
resultsNames(dds)
res <- results(dds, name = "condition_Positive_vs_Negative")
genes <- na.omit(res)
summary(genes)
View(as.data.frame(genes))
plotMA(genes)

# HEATMAP
hmap_genes <- cbind(genes, gene_id = rownames(genes))
hmap_genes <- hmap_genes[hmap_genes$padj < 0.1,]

mat <- counts(dds,normalized = T)[rownames(hmap_genes),]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(coldata)

Heatmap(mat.z, cluster_rows = T, cluster_columns = T, 
        column_labels = colnames(mat.z), name = "Z-score", 
        row_labels = hmap_genes[rownames(mat.z),]$gene_id)

# 5: Pre-processing DGE Results II ------------------
genes$regulation <- ifelse(genes$log2FoldChange > 0, "up", ifelse(genes$log2FoldChange < 0, "down", "unchanged"))
genes <- cbind(genes, gene_id = rownames(genes))
genes <- merge(genes, counts_data[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
genes$ENTREZID <- mapIds(org.Hs.eg.db, keys = genes$gene_id, column = "ENTREZID", keytype = "ENSEMBL")
Genes_table <- as.data.frame(genes)

# Merged data frame
merged_df <- merge(genes, counts_data, by = "gene_id", all.x = TRUE)
merged_df <- merged_df[, -11]
View(as.data.frame(merged_df))

# 6: Saving file
write.csv(genes, file = "GSE197178deseq.csv")
write.csv(merged_df, file = "GSE197178merged_table.csv")

# VOLCANO 
EnhancedVolcano(genes, x = "log2FoldChange", y = "padj", lab = genes$gene_name)

# 7: GENE ONTOLOGY ---------------

genes_to_test <- rownames(res[res$log2FoldChange < 0.1,])
genes_to_test
GO_results <- enrichGO(gene = genes_to_test, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL',
                       ont = 'BP')
GO_table <- as.data.frame(GO_results)
plot(barplot(GO_results, showCategory = 12))
write.csv(GO_table, file = "GSE197178_GO_table.csv")


# 8: KEGG pathway enrichment analysis----------
kegg_pathway <- enrichKEGG(
  gene = genes$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1)
kegg_table <- as.data.frame(kegg_pathway)
plot(barplot(kegg_pathway, showCategory = 15))
write.csv(kegg_table, file = "GSE197178_KEGG_table.csv")

# 9. GSEA -------------
GSEA_genes <- res[order(-res$log2FoldChange),]
gene_list <- GSEA_genes$log2FoldChange
names(gene_list) <- rownames(GSEA_genes)
gene_list

gse <- gseGO(gene_list, ont = 'BP', keyType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db',
             eps = 1e-300)
GSEA_table <- as.data.frame(gse)
write.csv(GSEA_table, file = "GSE197178_GSEA_table.csv")
gseaplot(gse, geneSetID = 1)

