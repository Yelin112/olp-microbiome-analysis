# ================================================================
# Oral 三状态网络图 + Eigengene × 菌群失调指数
# 在同一网络拓扑上展示 WT / DB / DB-Treated 三个状态
# 节点颜色 = 各状态下偏离WT的程度（基因：VST；菌属：CLR）
# 附图：Eigengene × MedianCLV 散点图（直接证明宿主-菌群关联）
# ================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(DESeq2)
  library(stringr)
  library(RColorBrewer)
})

set.seed(42)

# 关闭所有残留的图形设备（上次运行中断时可能遗留）
while (dev.cur() > 1) {
  dev.off()
}

# ────────────────────────────────────────────────────────────
# 0. 路径
# ────────────────────────────────────────────────────────────

BASE <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
RNA_DIR <- file.path(BASE, "analysis_results/RNA")
INT_DIR <- file.path(RNA_DIR, "Integration")
JNT_DIR <- file.path(BASE, "analysis_results/Joint/Oral")
COR_DIR <- file.path(JNT_DIR, "GeneGenus_Cor")
NET_DIR <- file.path(JNT_DIR, "Network")
OUT_DIR <- file.path(JNT_DIR, "Network_States")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ────────────────────────────────────────────────────────────
# 1. 读取数据（复用已有结果）
# ────────────────────────────────────────────────────────────

message("读取数据...")

edges_raw <- read.csv(
  file.path(COR_DIR, "sig_pairs_p05.csv"),
  stringsAsFactors = FALSE
)

rescue_genes <- read.csv(
  file.path(RNA_DIR, "Rescue/Oral/C_rescue_cluster_genes.csv"),
  stringsAsFactors = FALSE
)
cluster_sum <- read.csv(
  file.path(RNA_DIR, "Rescue/Oral/C_kmeans_cluster_summary.csv"),
  stringsAsFactors = FALSE
)
valid_clusters <- cluster_sum$cluster[
  cluster_sum$is_rescue == TRUE & cluster_sum$rescue_mag == TRUE
]
rescue_genes <- rescue_genes %>% filter(cluster %in% valid_clusters)

rescue_genus_df <- read.csv(
  file.path(JNT_DIR, "rescue_genus_list.csv"),
  stringsAsFactors = FALSE
) %>%
  filter(Genus != "Others")

# Eigengene 值（主脚本已算好）
eigengene_df <- read.csv(
  file.path(JNT_DIR, "eigengene_values.csv"),
  stringsAsFactors = FALSE
)

# 菌群失调评分
dysbiosis_df <- read.csv(
  file.path(BASE, "analysis_results/dysbiosis/dysbiosis_scores_Oral.csv"),
  stringsAsFactors = FALSE
)

# VST & 属丰度
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

genus_all <- read.csv(
  file.path(BASE, "analysis_results/tables/g_abundance_all.csv"),
  row.names = 1,
  check.names = FALSE
)
rownames(genus_all) <- gsub('"', '', rownames(genus_all))

# 样本对齐
rna_to_16s <- function(x) {
  p <- str_remove(x, "_Oral$")
  p <- str_replace(p, "([A-Za-z]+)(\\d+)$", "\\1.\\2")
  paste0(p, ".O")
}
meta_rna <- meta_rna %>% mutate(s16_sample = rna_to_16s(rna_sample))
cols_oral <- colnames(genus_all)[str_ends(colnames(genus_all), "\\.O")]
meta_matched <- meta_rna %>%
  filter(s16_sample %in% cols_oral) %>%
  arrange(Group)

vst_matched <- vst_mat[, meta_matched$rna_sample, drop = FALSE]
# 行名 ENSMUSG → gene_name
id2name <- setNames(rescue_genes$gene_name, rescue_genes$gene_id)
rownames(vst_matched) <- ifelse(
  is.na(id2name[rownames(vst_matched)]) | id2name[rownames(vst_matched)] == "",
  rownames(vst_matched),
  id2name[rownames(vst_matched)]
)

genus_site <- genus_all[, meta_matched$s16_sample, drop = FALSE]
clr_tf <- function(mat) {
  m <- as.matrix(mat) + 1e-6
  sweep(log(m), 2, colMeans(log(m)), "-")
}
genus_clr <- clr_tf(genus_site)

# ────────────────────────────────────────────────────────────
# 2. 计算三组均值 & 偏离WT程度
# ────────────────────────────────────────────────────────────

message("计算各组节点状态...")

grp3 <- c("WT_Control", "DB_Control", "DB_Treated")
states <- c("WT" = "WT_Control", "DB" = "DB_Control", "DBT" = "DB_Treated")

# 辅助：按组均值
grp_mean_gene <- function(gene, grp) {
  cols <- meta_matched$rna_sample[meta_matched$Group == grp]
  if (!gene %in% rownames(vst_matched)) {
    return(NA_real_)
  }
  mean(vst_matched[gene, cols])
}
grp_mean_genus <- function(genus, grp) {
  cols <- meta_matched$s16_sample[meta_matched$Group == grp]
  if (!genus %in% rownames(genus_clr)) {
    return(NA_real_)
  }
  mean(genus_clr[genus, cols])
}

# 所有网络节点
all_genes <- unique(edges_raw$Gene)
all_genera <- unique(edges_raw$Genus)
all_nodes <- c(all_genes, all_genera)
node_types <- c(
  rep("Gene", length(all_genes)),
  rep("Genus", length(all_genera))
)

# 每个节点 × 三组均值
node_means <- lapply(seq_along(all_nodes), function(i) {
  nd <- all_nodes[i]
  typ <- node_types[i]
  sapply(states, function(g) {
    if (typ == "Gene") {
      grp_mean_gene(nd, g)
    } else {
      grp_mean_genus(nd, g)
    }
  })
})
node_means_mat <- do.call(rbind, node_means)
rownames(node_means_mat) <- all_nodes
colnames(node_means_mat) <- names(states)

# 偏离WT的 Z-score（相对于 WT 均值归一化）
# delta = (value - WT_mean) / pooled_sd（跨节点标准化）
delta_mat <- node_means_mat - node_means_mat[, "WT"] # 各节点相对WT的偏差
# 全局标准化到 [-1, 1]（让颜色尺度统一）
max_abs <- max(abs(delta_mat), na.rm = TRUE)
delta_norm <- delta_mat / max_abs # 范围 [-1, 1]

# ────────────────────────────────────────────────────────────
# 3. 构建固定拓扑网络
# ────────────────────────────────────────────────────────────

message("构建网络拓扑...")

# 节点表（disease_dir用于形状）
dir_gene <- cluster_sum$disease_dir[
  match(
    rescue_genes$cluster[match(all_genes, rescue_genes$gene_name)],
    cluster_sum$cluster
  )
]
dir_genus <- sign(rescue_genus_df$delta_disease[
  match(all_genera, rescue_genus_df$Genus)
])

nodes_base <- data.frame(
  name = all_nodes,
  node_type = node_types,
  disease_dir = c(dir_gene, dir_genus),
  stringsAsFactors = FALSE
)

edges_base <- edges_raw %>%
  filter(Gene %in% all_nodes, Genus %in% all_nodes) %>%
  dplyr::select(from = Gene, to = Genus, r, p_raw) %>%
  mutate(
    edge_sign = ifelse(r > 0, "Positive", "Negative"),
    abs_r = abs(r)
  )

g_base <- tbl_graph(
  nodes = nodes_base,
  edges = edges_base,
  directed = FALSE
) %>%
  mutate(degree = centrality_degree())

# 固定布局（所有三张图共享）
layout_fixed <- create_layout(g_base, layout = "stress")
coords_fixed <- data.frame(
  name = layout_fixed$name,
  x = layout_fixed$x,
  y = layout_fixed$y
)

message(sprintf("  节点: %d  |  边: %d", vcount(g_base), ecount(g_base)))

# ────────────────────────────────────────────────────────────
# 4. 绘制三状态网络（WT / DB / DB-Treated）
# ────────────────────────────────────────────────────────────

message("绘制三状态网络图...")

GROUP_LABELS <- c(
  "WT" = "WT Control",
  "DB" = "DB Control",
  "DBT" = "DB Treated"
)

make_state_network <- function(state_key) {
  delta_vals <- delta_norm[all_nodes, state_key]
  nodes_s <- nodes_base %>%
    mutate(delta = delta_vals[name])

  g_s <- tbl_graph(nodes = nodes_s, edges = edges_base, directed = FALSE) %>%
    mutate(degree = centrality_degree())

  layout_s <- create_layout(
    g_s,
    layout = "manual",
    x = coords_fixed$x[match(V(g_s)$name, coords_fixed$name)],
    y = coords_fixed$y[match(V(g_s)$name, coords_fixed$name)]
  )

  ggraph(layout_s) +
    # 边（统一灰色，突出节点状态变化）
    geom_edge_link(
      aes(width = abs_r, alpha = abs_r),
      color = "grey60",
      show.legend = FALSE
    ) +
    scale_edge_width_continuous(range = c(0.4, 2.5)) +
    scale_edge_alpha_continuous(range = c(0.3, 0.8)) +

    # 节点：颜色=偏离WT程度，形状=基因/菌属
    geom_node_point(
      aes(shape = node_type, fill = delta, size = degree + 1),
      color = "grey30",
      stroke = 0.8
    ) +
    scale_shape_manual(
      values = c("Gene" = 21, "Genus" = 23),
      name = "Type",
      labels = c("Gene" = "Gene", "Genus" = "Genus")
    ) +
    scale_fill_gradientn(
      colors = c("#2166AC", "#F7F7F7", "#B2182B"),
      limits = c(-1, 1),
      name = "Deviation\nfrom WT",
      guide = guide_colorbar(barwidth = 0.8, barheight = 4)
    ) +
    scale_size_continuous(range = c(4, 10), guide = "none") +

    # 节点标签
    geom_node_text(
      aes(label = name),
      size = 3.0,
      repel = TRUE,
      max.overlaps = 25,
      bg.color = alpha("white", 0.6),
      bg.r = 0.1
    ) +

    labs(title = GROUP_LABELS[state_key]) +
    theme_graph(base_family = "", base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

p_wt <- make_state_network("WT")
p_db <- make_state_network("DB")
p_dbt <- make_state_network("DBT")

# ── 拼图：共享图例，统一色标 ────────────────────────────
p_three <- (p_wt + p_db + p_dbt) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Oral  Gene-Genus Network: Three Group States",
    subtitle = paste0(
      "Node color: deviation from WT_Control (blue = closer to WT, red = further from WT)\n",
      "Node shape: square = Gene, diamond = Genus | Edge width: |Spearman r| | p < 0.05"
    ),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave(
  file.path(OUT_DIR, "network_three_states.pdf"),
  p_three,
  width = 18,
  height = 7,
  device = "pdf"
)
ggsave(
  file.path(OUT_DIR, "network_three_states.png"),
  p_three,
  width = 18,
  height = 7,
  dpi = 200
)
message("  三状态网络图已保存")

# ────────────────────────────────────────────────────────────
# 5. 各节点偏离WT程度随组变化的轨迹热图
#    行=节点，列=WT/DB/DBT，颜色=偏离程度
# ────────────────────────────────────────────────────────────

message("绘制节点状态轨迹热图...")

library(pheatmap)

# 按节点类型和disease_dir排序
node_order <- nodes_base %>%
  arrange(node_type, disease_dir) %>%
  pull(name)
node_order <- node_order[node_order %in% rownames(delta_norm)]

delta_heat <- delta_norm[node_order, c("WT", "DB", "DBT")]
colnames(delta_heat) <- c("WT_Control", "DB_Control", "DB_Treated")

# 行注释
anno_row <- data.frame(
  Type = nodes_base$node_type[match(node_order, nodes_base$name)],
  row.names = node_order
)
anno_colors <- list(
  Type = c("Gene" = "#E64B35", "Genus" = "#4DBBD5")
)

pdf(
  file.path(OUT_DIR, "node_state_trajectory_heatmap.pdf"),
  width = 6,
  height = max(5, length(node_order) * 0.32 + 2.5)
)
pheatmap(
  delta_heat,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = "grey90",
  main = "Node State Trajectory: WT -> DB -> DB-Treated\n(Deviation from WT, normalized)"
)
dev.off()

png(
  file.path(OUT_DIR, "node_state_trajectory_heatmap.png"),
  width = 600,
  height = max(500, length(node_order) * 32 + 250),
  res = 120
)
pheatmap(
  delta_heat,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = "grey90",
  main = "Node State Trajectory: WT -> DB -> DB-Treated\n(Deviation from WT, normalized)"
)
dev.off()
message("  轨迹热图已保存")

# ────────────────────────────────────────────────────────────
# 6. Eigengene × 菌群失调指数（MedianCLV）散点图
# ────────────────────────────────────────────────────────────

message("绘制 Eigengene x 菌群失调指数...")

GROUP_COLORS <- c(
  "WT_Control" = "#4DBBD5",
  "WT_Treated" = "#00A087",
  "DB_Control" = "#E64B35",
  "DB_Treated" = "#F39B7F"
)

# 对齐 Eigengene 与失调评分
# eigengene_df 只有 rna_sample，需先转换为 16S 样本名再与 dysbiosis_df 匹配
eg_dys <- eigengene_df %>%
  mutate(s16_sample = rna_to_16s(rna_sample)) %>%
  left_join(
    dysbiosis_df %>%
      dplyr::select(SampleID_converted, MedianCLV, TreatGroup),
    by = c("s16_sample" = "SampleID_converted")
  ) %>%
  filter(!is.na(MedianCLV))

eg_cols <- colnames(eg_dys)[str_starts(colnames(eg_dys), "EG_")]
message(sprintf("  Eigengene 列: %s", paste(eg_cols, collapse = ", ")))
message(sprintf("  有效样本: %d", nrow(eg_dys)))

eg_labels <- c(
  "EG_DB_up" = "EG: DB(+) rescued",
  "EG_DB_down" = "EG: DB(-) rescued"
)

plot_eg_list <- lapply(eg_cols, function(ec) {
  df_p <- eg_dys %>%
    dplyr::select(
      sample = rna_sample,
      Group,
      Eigengene = all_of(ec),
      MedianCLV
    ) %>%
    filter(!is.na(Eigengene))

  # Spearman 相关
  cr <- cor.test(
    df_p$Eigengene,
    df_p$MedianCLV,
    method = "spearman",
    exact = FALSE
  )
  r_v <- round(cr$estimate, 3)
  p_v <- signif(cr$p.value, 2)
  sig <- ifelse(cr$p.value < 0.05, "*", ifelse(cr$p.value < 0.1, ".", "ns"))

  label_txt <- eg_labels[ec]
  if (is.na(label_txt)) {
    label_txt <- ec
  }

  ggplot(df_p, aes(Eigengene, MedianCLV, color = Group)) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "grey40",
      fill = "grey85",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_point(size = 4.5, alpha = 0.9) +
    # geom_text_repel(
    #   aes(label = sample),
    #   size = 2.5,
    #   max.overlaps = 12,
    #   segment.color = "grey60"
    # ) +
    scale_color_manual(values = GROUP_COLORS) +
    annotate(
      "text",
      x = min(df_p$Eigengene, na.rm = TRUE),
      y = max(df_p$MedianCLV, na.rm = TRUE),
      label = paste0("r = ", r_v, "\np = ", p_v, "  ", sig),
      hjust = 0,
      vjust = 1,
      size = 4,
      color = "grey30",
      fontface = "plain"
    ) +
    labs(
      title = paste0("Oral | ", label_txt, " x Dysbiosis Index"),
      subtitle = "x: Host Eigengene (PC1 of rescue gene cluster)\ny: MedianCLV (microbial dysbiosis score)",
      x = paste0(label_txt, " (VST PC1)"),
      y = "Dysbiosis Index (MedianCLV)",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "grey50"),
      legend.position = "bottom"
    )
})

if (length(plot_eg_list) > 0) {
  p_eg_combined <- wrap_plots(plot_eg_list, ncol = length(plot_eg_list)) +
    plot_annotation(
      title = "Oral: Host Rescue Eigengene x Microbial Dysbiosis Index",
      theme = theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    )
  ggsave(
    file.path(OUT_DIR, "eigengene_x_dysbiosis.pdf"),
    p_eg_combined,
    width = length(plot_eg_list) * 6,
    height = 6,
    device = "pdf"
  )
  ggsave(
    file.path(OUT_DIR, "eigengene_x_dysbiosis.png"),
    p_eg_combined,
    width = length(plot_eg_list) * 6,
    height = 6,
    dpi = 200
  )
  message("  Eigengene x 失调指数图已保存")
}

# ────────────────────────────────────────────────────────────
# 7. 汇总
# ────────────────────────────────────────────────────────────

cat(sprintf(
  "\n%s\n   三状态网络分析完成\n%s\n",
  strrep("=", 50),
  strrep("=", 50)
))
cat("\n   输出文件:\n")
cat(
  "   network_three_states.pdf          — WT/DB/DB-Treated 三状态网络并排图\n"
)
cat("   node_state_trajectory_heatmap.pdf — 各节点偏离WT程度轨迹热图\n")
cat("   eigengene_x_dysbiosis.pdf         — Eigengene x 菌群失调指数散点\n")
cat(sprintf("\n   输出目录: %s\n", OUT_DIR))
 