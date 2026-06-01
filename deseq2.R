# ============================================================
# DESeq2差异表达分析 - 使用IsoQuant结果 (修复版) isoquant.py --reference /data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/genome/osa_IRGSP_1.fa --genedb /data2/users/zhouyiming/DT4_ROOT/isoquant/osa_IRGSP_1.annotation.corrected.gtf --fastq /data2/users/zhouyiming/DT4_ROOT/root_ct1/drs-ct-R-1.pass.fq /data2/users/zhouyiming/DT4_ROOT/root_ct2/drs-ct-R-2.pass.fq /data2/users/zhouyiming/DT4_ROOT/root_dt4_1/drs-dt4-R-1.pass.fq /data2/users/zhouyiming/DT4_ROOT/root_dt4-2/drs-dt4-R-2.pass.fq --data_type nanopore -o /data2/users/zhouyiming/DT4_ROOT/isoquant
# ============================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)

# ============================================================
# 第一部分：设置路径和样本信息
# ============================================================

message("============================================================")
message("DESeq2差异表达分析 - IsoQuant结果")
message(paste0("开始时间: ", Sys.time()))
message("============================================================")

# IsoQuant结果目录
isoquant_dir <- "/data1/zhouyiming/DT4/DEG/OUT"

# 输出目录
output_dir <- "/data1/zhouyiming/DT4/DEG/deseq2_results_all"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 样本信息
sample_names <- c("ctrl_1", "ctrl_2", "dt4_1", "dt4_2")
conditions <- c("ctrl", "ctrl", "dt4", "dt4")

# IsoQuant的原始列名
isoquant_colnames <- c("ctrl-1-1", "ctrl-2-1", "DT4-1-1", "DT4-2-1")

sample_info <- data.frame(
  sample = sample_names,
  condition = factor(conditions, levels = c("ctrl", "dt4")),
  row.names = sample_names
)

message("\n样本信息:")
print(sample_info)

# ============================================================
# 第二部分：读取并处理IsoQuant counts数据
# ============================================================

message("\n=== 读取IsoQuant counts数据 ===")

# 读取基因水平counts
gene_file <- file.path(isoquant_dir, "OUT.discovered_gene_grouped_counts.tsv")
message(paste0("读取基因counts: ", gene_file))

gene_counts_raw <- read.table(gene_file, header = TRUE, sep = "\t", 
                               row.names = 1, check.names = FALSE,
                               stringsAsFactors = FALSE)

message(paste0("原始维度: ", nrow(gene_counts_raw), " x ", ncol(gene_counts_raw)))
message("原始列名:")
print(colnames(gene_counts_raw))
message("\n数据类型检查:")
print(sapply(gene_counts_raw, class))

# 读取转录本水平counts
transcript_file <- file.path(isoquant_dir, "OUT.transcript_grouped_counts.tsv")
message(paste0("\n读取转录本counts: ", transcript_file))

transcript_counts_raw <- read.table(transcript_file, header = TRUE, sep = "\t", 
                                     row.names = 1, check.names = FALSE,
                                     stringsAsFactors = FALSE)

message(paste0("原始维度: ", nrow(transcript_counts_raw), " x ", ncol(transcript_counts_raw)))

# ============================================================
# 第三部分：数据清洗和格式转换（关键修复）
# ============================================================

message("\n=== 数据清洗和格式转换 ===")

# 函数：清洗和转换counts数据
clean_counts <- function(counts_df, old_names, new_names) {
  
  # 只保留需要的列
  available_cols <- intersect(old_names, colnames(counts_df))
  message(paste0("可用的样本列: ", paste(available_cols, collapse = ", ")))
  
  if(length(available_cols) == 0) {
    stop("没有找到匹配的样本列！")
  }
  
  counts_subset <- counts_df[, available_cols, drop = FALSE]
  
  # 转换为矩阵
  counts_matrix <- as.matrix(counts_subset)
  
  # 确保是数值型
  counts_matrix <- apply(counts_matrix, 2, function(x) {
    # 将 "0" 字符转换为数值0，处理可能的NA
    x <- as.numeric(as.character(x))
    x[is.na(x)] <- 0
    return(x)
  })
  
  # 恢复行名
  rownames(counts_matrix) <- rownames(counts_df)
  
  # 重命名列
  name_mapping <- setNames(new_names, old_names)
  for(i in seq_along(colnames(counts_matrix))) {
    old_name <- colnames(counts_matrix)[i]
    if(old_name %in% names(name_mapping)) {
      colnames(counts_matrix)[i] <- name_mapping[old_name]
    }
  }
  
  # 按新名称顺序排列
  available_new_names <- intersect(new_names, colnames(counts_matrix))
  counts_matrix <- counts_matrix[, available_new_names, drop = FALSE]
  
  return(counts_matrix)
}

# 处理基因counts
gene_counts_matrix <- clean_counts(gene_counts_raw, isoquant_colnames, sample_names)
message(paste0("\n处理后基因counts: ", nrow(gene_counts_matrix), " x ", ncol(gene_counts_matrix)))
message("列名: ", paste(colnames(gene_counts_matrix), collapse = ", "))
message("数据类型: ", class(gene_counts_matrix[1,1]))

# 处理转录本counts
transcript_counts_matrix <- clean_counts(transcript_counts_raw, isoquant_colnames, sample_names)
message(paste0("处理后转录本counts: ", nrow(transcript_counts_matrix), " x ", ncol(transcript_counts_matrix)))

# 显示数据预览
message("\n基因counts预览（前5行）:")
print(head(gene_counts_matrix, 5))

# ============================================================
# 第四部分：保存各样本单独的counts文件
# ============================================================

message("\n=== 保存各样本单独的counts文件 ===")

for(sample in colnames(gene_counts_matrix)) {
  sample_dir <- file.path(output_dir, sample)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 基因counts
  gene_df <- data.frame(
    gene_id = rownames(gene_counts_matrix),
    counts = gene_counts_matrix[, sample]
  )
  colnames(gene_df)[2] <- sample
  write.csv(gene_df, file.path(sample_dir, paste0(sample, "_gene_counts.csv")), row.names = FALSE)
  
  # 转录本counts
  tx_df <- data.frame(
    transcript_id = rownames(transcript_counts_matrix),
    counts = transcript_counts_matrix[, sample]
  )
  colnames(tx_df)[2] <- sample
  write.csv(tx_df, file.path(sample_dir, paste0(sample, "_transcript_counts.csv")), row.names = FALSE)
  
  message(paste0("  已保存: ", sample))
}

# 保存合并的counts矩阵
write.csv(gene_counts_matrix, file.path(output_dir, "all_samples_gene_counts.csv"), row.names = TRUE)
write.csv(transcript_counts_matrix, file.path(output_dir, "all_samples_transcript_counts.csv"), row.names = TRUE)

# ============================================================
# 第五部分：DESeq2差异表达分析 - 基因水平
# ============================================================

message("\n=== DESeq2差异表达分析 - 基因水平 ===")

# 转换为整数矩阵（关键修复）
gene_counts_int <- matrix(
  as.integer(round(gene_counts_matrix)),
  nrow = nrow(gene_counts_matrix),
  ncol = ncol(gene_counts_matrix),
  dimnames = list(rownames(gene_counts_matrix), colnames(gene_counts_matrix))
)

message(paste0("整数矩阵维度: ", nrow(gene_counts_int), " x ", ncol(gene_counts_int)))
message(paste0("数据类型: ", typeof(gene_counts_int)))

# 过滤低表达基因（至少2个样本中counts >= 10）
keep_gene <- rowSums(gene_counts_int >= 10) >= 2
gene_counts_filtered <- gene_counts_int[keep_gene, , drop = FALSE]

message(paste0("过滤前基因数量: ", nrow(gene_counts_int)))
message(paste0("过滤后基因数量: ", nrow(gene_counts_filtered)))

# 更新样本信息以匹配可用样本
sample_info_sub <- sample_info[colnames(gene_counts_filtered), , drop = FALSE]

# 创建DESeq2对象
dds_gene <- DESeqDataSetFromMatrix(
  countData = gene_counts_filtered,
  colData = sample_info_sub,
  design = ~ condition
)

# 运行DESeq2
message("运行DESeq2...")
dds_gene <- DESeq(dds_gene)

# 保存DESeq2对象
saveRDS(dds_gene, file.path(output_dir, "dds_gene.rds"))

# 提取结果 (dt4 vs ctrl)
res_gene <- results(dds_gene, contrast = c("condition", "dt4", "ctrl"))
res_gene_df <- as.data.frame(res_gene[order(res_gene$padj), ])
res_gene_df$gene_id <- rownames(res_gene_df)

# 添加显著性标记
res_gene_df$significant <- ifelse(
  !is.na(res_gene_df$padj) & res_gene_df$padj < 0.05 & abs(res_gene_df$log2FoldChange) > 1,
  ifelse(res_gene_df$log2FoldChange > 1, "Up", "Down"),
  "NS"
)

# 重排列列顺序
res_gene_df <- res_gene_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant")]

# 保存完整结果
write.csv(res_gene_df, file.path(output_dir, "DESeq2_gene_level_results.csv"), row.names = FALSE)

# 统计显著差异基因
sig_gene_up <- sum(res_gene_df$significant == "Up", na.rm = TRUE)
sig_gene_down <- sum(res_gene_df$significant == "Down", na.rm = TRUE)

message(paste0("显著上调基因 (padj<0.05, log2FC>1): ", sig_gene_up))
message(paste0("显著下调基因 (padj<0.05, log2FC<-1): ", sig_gene_down))
message(paste0("总显著差异基因: ", sig_gene_up + sig_gene_down))

# 保存显著差异基因列表
sig_genes <- res_gene_df[res_gene_df$significant != "NS", ]
write.csv(sig_genes, file.path(output_dir, "significant_genes.csv"), row.names = FALSE)

up_genes <- res_gene_df[res_gene_df$significant == "Up", ]
down_genes <- res_gene_df[res_gene_df$significant == "Down", ]
write.csv(up_genes, file.path(output_dir, "upregulated_genes.csv"), row.names = FALSE)
write.csv(down_genes, file.path(output_dir, "downregulated_genes.csv"), row.names = FALSE)

# ============================================================
# 第六部分：DESeq2差异表达分析 - 转录本水平
# ============================================================

message("\n=== DESeq2差异表达分析 - 转录本水平 ===")

# 转换为整数矩阵
transcript_counts_int <- matrix(
  as.integer(round(transcript_counts_matrix)),
  nrow = nrow(transcript_counts_matrix),
  ncol = ncol(transcript_counts_matrix),
  dimnames = list(rownames(transcript_counts_matrix), colnames(transcript_counts_matrix))
)

# 过滤低表达转录本
keep_tx <- rowSums(transcript_counts_int >= 10) >= 2
transcript_counts_filtered <- transcript_counts_int[keep_tx, , drop = FALSE]

message(paste0("过滤前转录本数量: ", nrow(transcript_counts_int)))
message(paste0("过滤后转录本数量: ", nrow(transcript_counts_filtered)))

# 创建DESeq2对象
dds_tx <- DESeqDataSetFromMatrix(
  countData = transcript_counts_filtered,
  colData = sample_info_sub,
  design = ~ condition
)

# 运行DESeq2
message("运行DESeq2...")
dds_tx <- DESeq(dds_tx)

# 保存DESeq2对象
saveRDS(dds_tx, file.path(output_dir, "dds_transcript.rds"))

# 提取结果
res_tx <- results(dds_tx, contrast = c("condition", "dt4", "ctrl"))
res_tx_df <- as.data.frame(res_tx[order(res_tx$padj), ])
res_tx_df$transcript_id <- rownames(res_tx_df)

# 添加显著性标记
res_tx_df$significant <- ifelse(
  !is.na(res_tx_df$padj) & res_tx_df$padj < 0.05 & abs(res_tx_df$log2FoldChange) > 1,
  ifelse(res_tx_df$log2FoldChange > 1, "Up", "Down"),
  "NS"
)

# 重排列列顺序
res_tx_df <- res_tx_df[, c("transcript_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant")]

# 保存完整结果
write.csv(res_tx_df, file.path(output_dir, "DESeq2_transcript_level_results.csv"), row.names = FALSE)

# 统计
sig_tx_up <- sum(res_tx_df$significant == "Up", na.rm = TRUE)
sig_tx_down <- sum(res_tx_df$significant == "Down", na.rm = TRUE)

message(paste0("显著上调转录本 (padj<0.05, log2FC>1): ", sig_tx_up))
message(paste0("显著下调转录本 (padj<0.05, log2FC<-1): ", sig_tx_down))
message(paste0("总显著差异转录本: ", sig_tx_up + sig_tx_down))

# 保存显著差异转录本
sig_transcripts <- res_tx_df[res_tx_df$significant != "NS", ]
write.csv(sig_transcripts, file.path(output_dir, "significant_transcripts.csv"), row.names = FALSE)

# ============================================================
# 第七部分：可视化
# ============================================================

message("\n=== 生成可视化图表 ===")

# 1. PCA图
message("生成PCA图...")
vsd <- vst(dds_gene, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 5) +
  geom_text(vjust = -1.5, hjust = 0.5, size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  ggtitle("PCA Plot - Gene Level (IsoQuant)") +
  scale_color_manual(values = c("ctrl" = "#2166AC", "dt4" = "#B2182B"))

ggsave(file.path(output_dir, "PCA_plot.pdf"), pca_plot, width = 8, height = 6)
ggsave(file.path(output_dir, "PCA_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# 2. 火山图 - 基因水平
message("生成火山图（基因水平）...")

volcano_gene <- ggplot(res_gene_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  theme_bw(base_size = 14) +
  labs(
    title = "Leaf Volcano Plot - Gene Level (dt4 vs ctrl)",
    subtitle = paste0("Up: ", sig_gene_up, " | Down: ", sig_gene_down),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Volcano_plot_gene.pdf"), volcano_gene, width = 8, height = 6)
ggsave(file.path(output_dir, "Volcano_plot_gene.png"), volcano_gene, width = 8, height = 6, dpi = 300)

# 3. 火山图 - 转录本水平
message("生成火山图（转录本水平）...")

volcano_tx <- ggplot(res_tx_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot - Transcript Level (dt4 vs ctrl)",
    subtitle = paste0("Up: ", sig_tx_up, " | Down: ", sig_tx_down),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Volcano_plot_transcript.pdf"), volcano_tx, width = 8, height = 6)
ggsave(file.path(output_dir, "Volcano_plot_transcript.png"), volcano_tx, width = 8, height = 6, dpi = 300)

# 4. MA图
message("生成MA图...")
pdf(file.path(output_dir, "MA_plot_gene.pdf"), width = 8, height = 6)
plotMA(res_gene, main = "MA Plot - Gene Level (dt4 vs ctrl)", ylim = c(-5, 5))
dev.off()

# 5. 热图 - Top 50 差异基因
if(nrow(sig_genes) >= 2) {
  message("生成热图（Top差异基因）...")
  
  top_n <- min(50, nrow(sig_genes))
  top_genes <- head(sig_genes$gene_id, top_n)
  
  mat <- assay(vsd)[top_genes[top_genes %in% rownames(assay(vsd))], ]
  
  if(nrow(mat) >= 2) {
    mat_scaled <- t(scale(t(mat)))
    
    annotation_col <- data.frame(
      Condition = sample_info_sub$condition,
      row.names = rownames(sample_info_sub)
    )
    
    ann_colors <- list(
      Condition = c("ctrl" = "#2166AC", "dt4" = "#B2182B")
    )
    
    pdf(file.path(output_dir, "Heatmap_top_genes.pdf"), width = 8, height = 12)
    pheatmap(mat_scaled,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             show_rownames = TRUE,
             show_colnames = TRUE,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize_row = 8,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             main = paste0("Top ", nrow(mat), " Differentially Expressed Genes"))
    dev.off()
    
    png(file.path(output_dir, "Heatmap_top_genes.png"), width = 800, height = 1200, res = 100)
    pheatmap(mat_scaled,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             show_rownames = TRUE,
             show_colnames = TRUE,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize_row = 8,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             main = paste0("Top ", nrow(mat), " Differentially Expressed Genes"))
    dev.off()
  }
}

# 6. 样本相关性热图
message("生成样本相关性热图...")
sample_cor <- cor(assay(vsd))

annotation_col <- data.frame(
  Condition = sample_info_sub$condition,
  row.names = rownames(sample_info_sub)
)

ann_colors <- list(
  Condition = c("ctrl" = "#2166AC", "dt4" = "#B2182B")
)

pdf(file.path(output_dir, "Sample_correlation_heatmap.pdf"), width = 6, height = 5)
pheatmap(sample_cor,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         display_numbers = TRUE,
         number_format = "%.3f",
         color = colorRampPalette(c("white", "#B2182B"))(100),
         main = "Sample Correlation")
dev.off()

# ============================================================
# 第八部分：生成分析报告
# ============================================================

message("\n=== 生成分析报告 ===")

summary_report <- paste0(
  "========================================\n",
  "IsoQuant + DESeq2 分析报告\n",
  "========================================\n",
  "分析完成时间: ", Sys.time(), "\n",
  "----------------------------------------\n\n",
  "【样本信息】\n",
  "对照组 (ctrl): ctrlR_1, ctrlR_2\n",
  "处理组 (dt4):  dt4R_1, dt4R_2\n\n",
  "【定量结果】\n",
  "总基因数量: ", nrow(gene_counts_matrix), "\n",
  "总转录本数量: ", nrow(transcript_counts_matrix), "\n",
  "过滤后基因数量: ", nrow(gene_counts_filtered), "\n",
  "过滤后转录本数量: ", nrow(transcript_counts_filtered), "\n\n",
  "【差异分析 - 基因水平】\n",
  "上调基因: ", sig_gene_up, "\n",
  "下调基因: ", sig_gene_down, "\n",
  "总差异基因: ", sig_gene_up + sig_gene_down, "\n\n",
  "【差异分析 - 转录本水平】\n",
  "上调转录本: ", sig_tx_up, "\n",
  "下调转录本: ", sig_tx_down, "\n",
  "总差异转录本: ", sig_tx_up + sig_tx_down, "\n",
  "========================================\n"
)

cat(summary_report)
writeLines(summary_report, file.path(output_dir, "analysis_summary.txt"))

message("\n============================================================")
message("所有分析完成!")
message(paste0("结果保存在: ", output_dir))
message("============================================================")
