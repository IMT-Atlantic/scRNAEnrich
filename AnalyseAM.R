setwd("D:/Rproject/scRNA/Irg1KO")

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

Data <- readRDS('E:/BioBank/scRNA/GSE215195ITAscRNA/DataIndex.rds')
Data <- JoinLayers(Data)
DimPlot(Data, reduction = "umap", group.by = 'cell_type_annotation', label=TRUE)

AM <- subset(Data, subset = cell_type_annotation == "Alveolar Macrophage")
AM <- FindNeighbors(object = AM, reduction = "pca", dims = 1:7)
AM <- FindClusters(object = AM, resolution = 2)
AM <- RunUMAP(object = AM, dims = 1:7)
DimPlot(AM, reduction = "umap", group.by = 'orig.ident', label=TRUE)

# 对不同组进行差异表达分析
groups <- c("PBS_WT", "LAC_WT", "PBS_IRG1KO", "LAC_IRG1KO")
diff_genes_list <- list()
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    group1 <- groups[i]
    group2 <- groups[j]
    markers <- FindMarkers(AM, ident.1 = group1, ident.2 = group2, group.by = "orig.ident", test.use = "wilcox")
    diff_genes_list[[paste(group1, "vs", group2)]] <- markers
  }
}

process_diff_genes <- function(markers, group1, group2) {
  if (!"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)  # 将基因名从行名提取到 gene 列
  }
  markers <- markers %>%
    mutate(Change = case_when(
      avg_log2FC > 0 ~ "Up",
      avg_log2FC < 0 ~ "Down",
      TRUE ~ "No Change"
    ))
  up_genes <- markers %>% filter(Change == "Up") %>% pull(gene)
  down_genes <- markers %>% filter(Change == "Down") %>% pull(gene)
  up_genes_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  down_genes_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  up_genes_entrez_ids <- up_genes_entrez$ENTREZID
  down_genes_entrez_ids <- down_genes_entrez$ENTREZID
  kegg_up <- enrichKEGG(gene = up_genes_entrez_ids, organism = 'mmu', keyType = 'kegg')
  kegg_down <- enrichKEGG(gene = down_genes_entrez_ids, organism = 'mmu', keyType = 'kegg')
  go_up <- enrichGO(gene = up_genes_entrez_ids, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = 'ENTREZID')
  go_down <- enrichGO(gene = down_genes_entrez_ids, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = 'ENTREZID')
  list(kegg_up = kegg_up, kegg_down = kegg_down, go_up = go_up, go_down = go_down)
}

# 创建一个空列表来存储富集分析结果
enrichment_results <- list()

# 对每对组进行差异分析和富集分析
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    group1 <- groups[i]
    group2 <- groups[j]
    markers <- diff_genes_list[[paste(group1, "vs", group2)]]
    enrichment_results[[paste(group1, "vs", group2)]] <- process_diff_genes(markers, group1, group2)
    p1 <- dotplot(enrichment_results[[paste(group1, "vs", group2)]]$kegg_up)
    p2 <- dotplot(enrichment_results[[paste(group1, "vs", group2)]]$kegg_down)
    p3 <- dotplot(enrichment_results[[paste(group1, "vs", group2)]]$go_up)
    p4 <- dotplot(enrichment_results[[paste(group1, "vs", group2)]]$go_down)
    ggsave(paste0(group1, "_vs_", group2, "_kegg_up.png"), plot = p1, width = 10, height = 10, dpi = 800)
    ggsave(paste0(group1, "_vs_", group2, "_kegg_down.png"), plot = p2, width = 10, height = 10, dpi = 800)
    ggsave(paste0(group1, "_vs_", group2, "_go_up.png"), plot = p3, width = 10, height = 10, dpi = 800)
    ggsave(paste0(group1, "_vs_", group2, "_go_down.png"), plot = p4, width = 10, height = 10, dpi = 800)
  }
}
