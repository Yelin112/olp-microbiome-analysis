# ================================================================
# Oral 挽救基因 × 挽救菌属 基因级别相关性分析
# 目标：34个挽救基因 × 13个挽救菌属，识别协同变化的基因-菌属对
# 方法：Spearman 相关 + FDR 校正
# ================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(Hmisc)
  library(stringr)
  library(patchwork)
})

set.seed(42)

# ────────────────────────────────────────────────────────────
# 0. 路径
# ────────────────────────────────────────────────────────────

BASE <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
RNA_DIR <- file.path(BASE, "analysis_results/RNA")
INT_DIR <- file.path(RNA_DIR, "Integration")
JNT_DIR <- file.path(BASE, "analysis_results/Joint/Oral")
OUT_DIR <- file.path(BASE, "analysis_results/Joint/Oral/GeneGenus_Cor")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ────────────────────────────────────────────────────────────
# 1. 读取数据
# ────────────────────────────────────────────────────────────

message("读取 Oral 挽救基因...")
rescue_genes <- read.csv(
  file.path(RNA_DIR, "Rescue/Oral/C_rescue_cluster_genes.csv"),
  stringsAsFactors = FALSE
)
# 仅保留 is_rescue=TRUE, rescue_mag=TRUE 的簇（2和4）
rescue_cluster_summary <- read.csv(
  file.path(RNA_DIR, "Rescue/Oral/C_kmeans_cluster_summary.csv"),
  stringsAsFactors = FALSE
)
valid_clusters <- rescue_cluster_summary$cluster[
  rescue_cluster_summary$is_rescue == TRUE &
    rescue_cluster_summary$rescue_mag == TRUE
]
rescue_genes <- rescue_genes %>%
  filter(cluster %in% valid_clusters)
message(sprintf(
  "  挽救基因: %d 个（簇 %s）",
  nrow(rescue_genes),
  paste(valid_clusters, collapse = ", ")
))

message("读取 Oral 挽救菌属...")
rescue_genus_df <- read.csv(
  file.path(JNT_DIR, "rescue_genus_list.csv"),
  stringsAsFactors = FALSE
) %>%
  filter(Genus != "Others") # 排除"Others"分类
message(sprintf("  挽救菌属: %d 个", nrow(rescue_genus_df)))
print(rescue_genus_df$Genus)

message("读取 VST 表达矩阵...")
vsd <- readRDS(file.path(INT_DIR, "Oral_vsd.rds"))
vst_mat <- assay(vsd)
meta_rna <- read.csv(
  file.path(INT_DIR, "Oral_metadata.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("WT_Control", "WT_Treated", "DB_Control", "DB_Treated")
    )
  )

message("读取属水平丰度（Oral）...")
genus_all <- read.csv(
  file.path(BASE, "analysis_results/tables/g_abundance_all.csv"),
  row.names = 1,
  check.names = FALSE
)
rownames(genus_all) <- gsub('"', '', rownames(genus_all))

# ────────────────────────────────────────────────────────────
# 2. 样本对齐：RNA ↔ 16S（Oral，后缀 .O）
# ────────────────────────────────────────────────────────────

# RNA样本名 → 16S样本名
rna_to_16s_oral <- function(rna_name) {
  prefix <- str_remove(rna_name, "_Oral$")
  prefix_dot <- str_replace(prefix, "([A-Za-z]+)(\\d+)$", "\\1.\\2")
  paste0(prefix_dot, ".O")
}

meta_rna <- meta_rna %>%
  mutate(s16_sample = rna_to_16s_oral(rna_sample))

# 筛选有16S数据的样本
cols_oral <- colnames(genus_all)[str_ends(colnames(genus_all), "\\.O")]
meta_matched <- meta_rna %>%
  filter(s16_sample %in% cols_oral) %>%
  arrange(Group)

message(sprintf("RNA-16S匹配样本: %d / %d", nrow(meta_matched), nrow(meta_rna)))
if (nrow(meta_matched) < 6) {
  stop("匹配样本过少，请检查样本名对应关系")
}

# 对齐矩阵
vst_matched <- vst_mat[, meta_matched$rna_sample, drop = FALSE]
genus_site <- genus_all[, meta_matched$s16_sample, drop = FALSE]

# ────────────────────────────────────────────────────────────
# 3. 提取挽救基因 VST 矩阵 & 挽救菌属 CLR 矩阵
# ────────────────────────────────────────────────────────────

# 挽救基因 VST（基因 × 样本）
gene_ids_use <- intersect(rescue_genes$gene_id, rownames(vst_matched))
vst_rescue <- vst_matched[gene_ids_use, , drop = FALSE]

# 用 gene_name 作行名（更易读）
gn <- rescue_genes$gene_name[match(rownames(vst_rescue), rescue_genes$gene_id)]
rownames(vst_rescue) <- ifelse(is.na(gn) | gn == "", rownames(vst_rescue), gn)

message(sprintf("可用挽救基因: %d", nrow(vst_rescue)))

# 挽救菌属 CLR（菌属 × 样本）
genus_rescue_names <- rescue_genus_df$Genus
genus_names_in_mat <- intersect(genus_rescue_names, rownames(genus_site))
genus_rescue_raw <- genus_site[genus_names_in_mat, , drop = FALSE]

# CLR 转换
clr_transform <- function(mat) {
  mat_pseudo <- as.matrix(mat) + 1e-6
  log_mat <- log(mat_pseudo)
  sweep(log_mat, 2, colMeans(log_mat), "-")
}
genus_rescue_clr <- clr_transform(genus_rescue_raw)

message(sprintf("可用挽救菌属: %d", nrow(genus_rescue_clr)))

# ────────────────────────────────────────────────────────────
# 4. Spearman 相关：所有基因 × 所有菌属
# ────────────────────────────────────────────────────────────

message("计算 Spearman 相关（基因 × 菌属）...")

# 合并矩阵（样本为行）
n_gene <- nrow(vst_rescue)
n_genus <- nrow(genus_rescue_clr)

combined <- t(rbind(vst_rescue, genus_rescue_clr)) # 样本 × (基因+菌属)

cor_full <- rcorr(combined, type = "spearman")
cor_block <- cor_full$r[1:n_gene, (n_gene + 1):(n_gene + n_genus), drop = FALSE]
p_block <- cor_full$P[1:n_gene, (n_gene + 1):(n_gene + n_genus), drop = FALSE]

# FDR 校正（BH，对所有基因-菌属对）
p_adj_vec <- p.adjust(as.vector(p_block), method = "BH")
p_adj_mat <- matrix(
  p_adj_vec,
  nrow = n_gene,
  dimnames = list(rownames(cor_block), colnames(cor_block))
)

# ────────────────────────────────────────────────────────────
# 5. 汇总显著对
# ────────────────────────────────────────────────────────────

pairs_df <- expand.grid(
  Gene = rownames(cor_block),
  Genus = colnames(cor_block),
  stringsAsFactors = FALSE
) %>%
  mutate(
    r = as.vector(cor_block),
    p_raw = as.vector(p_block),
    p_adj = as.vector(p_adj_mat),
    # 基因方向（来自挽救簇的disease_dir）
    gene_dir = rescue_genes$cluster[
      match(
        rescue_genes$gene_name[match(Gene, rescue_genes$gene_name)],
        rescue_genes$gene_name
      )
    ],
    # 菌属挽救方向（delta_disease 的符号）
    genus_dir = sign(rescue_genus_df$delta_disease[
      match(Genus, rescue_genus_df$Genus)
    ])
  ) %>%
  filter(!is.na(r), !is.na(p_raw)) %>%
  arrange(p_raw)

write.csv(
  pairs_df,
  file.path(OUT_DIR, "all_gene_genus_cor.csv"),
  row.names = FALSE
)

# 显著对（放宽阈值：p_raw < 0.05，因为小样本FDR极严）
pairs_sig_raw <- pairs_df %>% filter(p_raw < 0.05)
pairs_sig_padj <- pairs_df %>% filter(p_adj < 0.2)

message(sprintf(
  "  p_raw<0.05: %d 对  |  FDR<0.2: %d 对",
  nrow(pairs_sig_raw),
  nrow(pairs_sig_padj)
))

# 写出两个阈值的结果
write.csv(
  pairs_sig_raw,
  file.path(OUT_DIR, "sig_pairs_p05.csv"),
  row.names = FALSE
)
write.csv(
  pairs_sig_padj,
  file.path(OUT_DIR, "sig_pairs_FDR02.csv"),
  row.names = FALSE
)

# ────────────────────────────────────────────────────────────
# 6. 可视化 A：相关热图（基因 × 菌属）
# ────────────────────────────────────────────────────────────

message("绘制相关热图...")

# 标记符号矩阵
sig_mark <- matrix(
  "",
  nrow = n_gene,
  ncol = n_genus,
  dimnames = list(rownames(cor_block), colnames(cor_block))
)
sig_mark[p_block < 0.1] <- "."
sig_mark[p_block < 0.05] <- "*"
sig_mark[p_adj_mat < 0.2] <- "**"

# 菌属方向注释（挽救方向：上调↑ or 下调↓）
genus_dir_label <- ifelse(
  rescue_genus_df$delta_disease[match(
    colnames(cor_block),
    rescue_genus_df$Genus
  )] >
    0,
  paste0(colnames(cor_block), " ↑DB"),
  paste0(colnames(cor_block), " ↓DB")
)
colnames_annotated <- genus_dir_label

# 基因按簇排序（先簇2后簇4，或按disease_dir）
gene_order <- rescue_genes %>%
  filter(gene_name %in% rownames(cor_block)) %>%
  arrange(cluster, desc(abs(log2FoldChange))) %>%
  pull(gene_name)
gene_order <- gene_order[gene_order %in% rownames(cor_block)]

cor_plot <- cor_block[gene_order, , drop = FALSE]
sig_plot <- sig_mark[gene_order, , drop = FALSE]
colnames(cor_plot) <- genus_dir_label
colnames(sig_plot) <- genus_dir_label

# 基因行注释（属于哪个簇 / disease_dir）
anno_row <- data.frame(
  Cluster = factor(rescue_genes$cluster[
    match(gene_order, rescue_genes$gene_name)
  ]),
  row.names = gene_order
)
anno_colors <- list(
  Cluster = setNames(
    c("#E64B35", "#4DBBD5"),
    as.character(sort(unique(anno_row$Cluster)))
  )
)

ph <- max(7, n_gene * 0.35 + 3)
pw <- max(6, n_genus * 0.9 + 4)

pdf(file.path(OUT_DIR, "heatmap_gene_genus_cor.pdf"), width = pw, height = ph)
pheatmap(
  cor_plot,
  display_numbers = sig_plot,
  number_color = "black",
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 9,
  border_color = "grey90",
  main = paste0(
    "Oral  挽救基因 × 挽救菌属  Spearman 相关\n",
    "（. p<0.1  * p<0.05  ** FDR<0.2）"
  )
)
dev.off()

png(
  file.path(OUT_DIR, "heatmap_gene_genus_cor.png"),
  width = pw * 120,
  height = ph * 120,
  res = 120
)
pheatmap(
  cor_plot,
  display_numbers = sig_plot,
  number_color = "black",
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 9,
  border_color = "grey90",
  main = paste0(
    "Oral  挽救基因 × 挽救菌属  Spearman 相关\n",
    "（. p<0.1  * p<0.05  ** FDR<0.2）"
  )
)
dev.off()
message("  热图已保存")

# ────────────────────────────────────────────────────────────
# 7. 可视化 B：显著对气泡图
# ────────────────────────────────────────────────────────────

message("绘制显著对气泡图...")

plot_pairs <- pairs_df %>%
  filter(p_raw < 0.05) %>%
  mutate(
    sig_level = case_when(
      p_adj < 0.2 ~ "FDR<0.2",
      p_raw < 0.01 ~ "p<0.01",
      TRUE ~ "p<0.05"
    ),
    direction = ifelse(r > 0, "正相关", "负相关")
  )

if (nrow(plot_pairs) > 0) {
  # 按菌属排序
  genus_order_bubble <- plot_pairs %>%
    group_by(Genus) %>%
    summarise(mean_abs_r = mean(abs(r)), .groups = "drop") %>%
    arrange(desc(mean_abs_r)) %>%
    pull(Genus)

  gene_order_bubble <- plot_pairs %>%
    group_by(Gene) %>%
    summarise(mean_abs_r = mean(abs(r)), .groups = "drop") %>%
    arrange(desc(mean_abs_r)) %>%
    pull(Gene)

  plot_pairs <- plot_pairs %>%
    mutate(
      Gene = factor(Gene, levels = gene_order_bubble),
      Genus = factor(Genus, levels = genus_order_bubble)
    )

  p_bubble <- ggplot(
    plot_pairs,
    aes(x = Genus, y = Gene, size = abs(r), color = r, shape = sig_level)
  ) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(
      colors = c("#2166AC", "white", "#B2182B"),
      limits = c(-1, 1),
      name = "Spearman r"
    ) +
    scale_size_continuous(range = c(2, 8), name = "|r|") +
    scale_shape_manual(
      values = c("FDR<0.2" = 18, "p<0.01" = 17, "p<0.05" = 16),
      name = "显著性"
    ) +
    labs(
      title = "Oral  挽救基因 × 挽救菌属  显著相关对",
      subtitle = "仅显示 p_raw < 0.05 的配对\n◆ FDR<0.2  ▲ p<0.01  ● p<0.05",
      x = "挽救菌属",
      y = "挽救基因"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.subtitle = element_text(size = 9)
    )

  ph_b <- max(5, length(gene_order_bubble) * 0.38 + 3)
  pw_b <- max(6, length(genus_order_bubble) * 0.85 + 4)
  ggsave(
    file.path(OUT_DIR, "bubble_sig_pairs.pdf"),
    p_bubble,
    width = pw_b,
    height = ph_b
  )
  ggsave(
    file.path(OUT_DIR, "bubble_sig_pairs.png"),
    p_bubble,
    width = pw_b,
    height = ph_b,
    dpi = 200
  )
  message("  气泡图已保存")
} else {
  message("  无 p<0.05 的显著对，跳过气泡图")
}

# ────────────────────────────────────────────────────────────
# 8. 可视化 C：top显著对的散点趋势图（每对一个小图）
# ────────────────────────────────────────────────────────────

message("绘制 top 显著对散点图...")

GROUP_COLORS <- c(
  "WT_Control" = "#4DBBD5",
  "WT_Treated" = "#00A087",
  "DB_Control" = "#E64B35",
  "DB_Treated" = "#F39B7F"
)

# 取 |r| 最大的前12对（或全部若 < 12）
top_pairs <- pairs_df %>%
  filter(p_raw < 0.05) %>%
  arrange(p_raw, desc(abs(r))) %>%
  head(12)

if (nrow(top_pairs) > 0) {
  plot_list <- lapply(seq_len(nrow(top_pairs)), function(i) {
    gene_i <- top_pairs$Gene[i]
    genus_i <- top_pairs$Genus[i]
    r_i <- round(top_pairs$r[i], 3)
    p_i <- signif(top_pairs$p_raw[i], 2)
    padj_i <- signif(top_pairs$p_adj[i], 2)

    df_i <- data.frame(
      x = as.numeric(genus_rescue_clr[genus_i, ]),
      y = as.numeric(vst_rescue[gene_i, ]),
      Group = meta_matched$Group,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        Group = factor(
          Group,
          levels = c("WT_Control", "WT_Treated", "DB_Control", "DB_Treated")
        )
      )

    ggplot(df_i, aes(x, y, color = Group)) +
      geom_point(size = 3, alpha = 0.85) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "grey40",
        linewidth = 0.7,
        linetype = "dashed"
      ) +
      scale_color_manual(values = GROUP_COLORS) +
      labs(
        title = paste0(gene_i, " ~ ", str_trunc(genus_i, 20)),
        subtitle = paste0("r=", r_i, "  p=", p_i, "  FDR=", padj_i),
        x = paste0("CLR (", str_trunc(genus_i, 18), ")"),
        y = paste0("VST (", gene_i, ")")
      ) +
      theme_bw(base_size = 9) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.5),
        axis.title = element_text(size = 8)
      )
  })

  # 拼图
  n_col <- min(3, nrow(top_pairs))
  n_row <- ceiling(nrow(top_pairs) / n_col)
  p_grid <- wrap_plots(plot_list, ncol = n_col) +
    plot_annotation(
      title = "Oral  Top 显著基因-菌属对  散点图",
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
      )
    )

  ggsave(
    file.path(OUT_DIR, "scatter_top_pairs.pdf"),
    p_grid,
    width = n_col * 4,
    height = n_row * 3.8
  )
  ggsave(
    file.path(OUT_DIR, "scatter_top_pairs.png"),
    p_grid,
    width = n_col * 4,
    height = n_row * 3.8,
    dpi = 200
  )
  message("  散点图已保存")
} else {
  message("  无显著对，跳过散点图")
}

# ────────────────────────────────────────────────────────────
# 9. 汇总输出
# ────────────────────────────────────────────────────────────

cat(sprintf(
  "\n%s\n   Oral 基因-菌属相关分析完成\n%s\n",
  strrep("=", 50),
  strrep("=", 50)
))
cat(sprintf(
  "   挽救基因: %d  |  挽救菌属: %d\n",
  nrow(vst_rescue),
  nrow(genus_rescue_clr)
))
cat(sprintf("   匹配样本: %d\n", nrow(meta_matched)))
cat(sprintf("   总配对数: %d\n", nrow(pairs_df)))
cat(sprintf("   p<0.05:  %d 对\n", nrow(pairs_sig_raw)))
cat(sprintf("   FDR<0.2: %d 对\n", nrow(pairs_sig_padj)))
cat(sprintf("\n   输出目录: %s\n", OUT_DIR))
cat("\n   输出文件:\n")
cat("   all_gene_genus_cor.csv      — 全部配对相关结果\n")
cat("   sig_pairs_p05.csv           — p<0.05 显著对\n")
cat("   sig_pairs_FDR02.csv         — FDR<0.2 显著对\n")
cat("   heatmap_gene_genus_cor.pdf  — 相关热图（基因×菌属）\n")
cat("   bubble_sig_pairs.pdf        — 显著对气泡图\n")
cat("   scatter_top_pairs.pdf       — Top对散点图（含趋势线）\n")
