setwd("D:/Rproject/scRNA/ZenghanAMU")

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)  
library(ggplot2)
library(readr)         
library(tidyr)  

Epithelial <- readRDS('E:/BioBank/scRNA/ZenghanAMU/Epithelial.rds')

diff_genes <- FindMarkers(Epithelial, 
                          ident.1 = "MHC2_high", 
                          ident.2 = "MHC2_low", 
                          group.by = "CellTypeMHC2",
                          min.pct = 0.25,  # Set minimum percentage of cells expressing the gene
                          logfc.threshold = 0.25)  # Set minimum log fold-change threshold

diff_genes <- diff_genes %>%
  mutate(Change = ifelse(avg_log2FC > 0, "Up", "Down"))
up_genes <- diff_genes %>% filter(Change == "Up") %>% rownames()
down_genes <- diff_genes %>% filter(Change == "Down") %>% rownames()
up_genes_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_genes_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
up_genes_entrez_ids <- up_genes_entrez$ENTREZID
down_genes_entrez_ids <- down_genes_entrez$ENTREZID
kegg_up <- enrichKEGG(gene = up_genes_entrez_ids, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_genes_entrez_ids, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05)
go_up <- enrichGO(gene = up_genes_entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = 'ENTREZID', pvalueCutoff = 0.05)
go_down <- enrichGO(gene = down_genes_entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = 'ENTREZID', pvalueCutoff = 0.05)
p1 <- dotplot(kegg_up, showCategory = 10) + ggtitle("KEGG Pathways Enriched in Upregulated Genes")
p2 <- dotplot(kegg_down, showCategory = 10) + ggtitle("KEGG Pathways Enriched in Downregulated Genes")
p3 <- dotplot(go_up, showCategory = 10) + ggtitle("GO Terms Enriched in Upregulated Genes")
p4 <- dotplot(go_down, showCategory = 10) + ggtitle("GO Terms Enriched in Downregulated Genes")
ggsave("kegg_up.png", plot = p1, width = 10, height = 10, dpi = 800)
ggsave("kegg_down.png", plot = p2, width = 10, height = 10, dpi = 800)
ggsave("go_up.png", plot = p3, width = 10, height = 10, dpi = 800)
ggsave("go_down.png", plot = p4, width = 10, height = 10, dpi = 800)

################################################################################
# Define function for analyse Basal,Luminal and Basal
find_diff_genes_MHC2 <- function(seurat_obj, celltype_bcl, 
                                 mhc2_col = "CellTypeMHC2", 
                                 bcl_col = "predicted_celltype_BCL", 
                                 group_high = "MHC2_high", 
                                 group_low = "MHC2_low", 
                                 min_pct = 0.25, 
                                 logfc_thresh = 0.25,
                                 output_dir = "Differential_Genes") {
  # 检查列名是否存在
  if(!(mhc2_col %in% colnames(seurat_obj@meta.data))){
    stop(paste("列名", mhc2_col, "在Seurat对象的元数据中未找到。"))
  }
  
  if(!(bcl_col %in% colnames(seurat_obj@meta.data))){
    stop(paste("列名", bcl_col, "在Seurat对象的元数据中未找到。"))
  }
  
  # 创建输出目录（如果不存在）
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # 筛选符合条件的细胞名称
  cells_high <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[mhc2_col]] == group_high & 
                                                 seurat_obj@meta.data[[bcl_col]] == celltype_bcl]
  cells_low <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[mhc2_col]] == group_low & 
                                                seurat_obj@meta.data[[bcl_col]] == celltype_bcl]
  
  cat(paste("正在处理细胞类型:", celltype_bcl, "\n"))
  cat("  MHC2_high 细胞数量:", length(cells_high), "\n")
  cat("  MHC2_low 细胞数量:", length(cells_low), "\n")
  
  # 检查细胞数量
  if(length(cells_high) < 3 | length(cells_low) < 3){
    warning(paste("在", celltype_bcl, "中，MHC2_high 或 MHC2_low 细胞数量少于3，可能不适合进行差异基因分析。"))
    return(NULL)
  }
  
  # 创建一个子集Seurat对象，只包含高低MHC2组的细胞
  subset_obj <- subset(seurat_obj, cells = c(cells_high, cells_low))
  
  # 添加一个新的元数据列用于分组
  subset_obj$MHC2_Group <- ifelse(subset_obj@meta.data[[mhc2_col]] == group_high, group_high, group_low)
  
  # 设置临时身份为 MHC2_Group
  Idents(subset_obj) <- "MHC2_Group"
  
  # 执行差异基因分析
  diff_genes <- FindMarkers(subset_obj, 
                            ident.1 = group_high, 
                            ident.2 = group_low, 
                            min.pct = min_pct, 
                            logfc.threshold = logfc_thresh)
  
  # 添加基因名作为一列
  diff_genes$gene <- rownames(diff_genes)
  
  # 保存结果到CSV文件
  output_file <- file.path(output_dir, paste0("Differential_Genes_MHC2_", celltype_bcl, ".csv"))
  write.csv(diff_genes, output_file, row.names = FALSE)
  
  cat(paste("差异基因已保存到:", output_file, "\n\n"))
  
  return(diff_genes)
}
################################################################################

# 使用函数对多个细胞类型进行比较
celltype_list <- c("Luminal", "Basal", "Club")
diff_genes_list <- list()
for(celltype in celltype_list){
  diff_genes <- find_diff_genes_MHC2(seurat_obj = Epithelial, 
                                     celltype_bcl = celltype, 
                                     output_dir = "Differential_Genes")
  if(!is.null(diff_genes)){
    diff_genes_list[[celltype]] <- diff_genes
  }
}
