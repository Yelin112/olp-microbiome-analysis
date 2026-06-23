# ================================================================
# Oral 挽救基因 × 挽救菌属 网络可视化
# 展示 WT → DB → DB-Treated 的协同变化模式
#
# 设计：
#   节点 = 基因（方形）或菌属（圆形）
#   边   = 显著相关（宽度 ∝ |r|，颜色 = 正/负相关）
#   节点颜色 = WT→DB→DB-Treated 的变化方向
#   节点大小 = 在网络中的连接度（degree）
#   附图：各节点三组表达/丰度趋势折线图
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

# ────────────────────────────────────────────────────────────
# 0. 路径
# ────────────────────────────────────────────────────────────

BASE <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
RNA_DIR <- file.path(BASE, "analysis_results/RNA")
INT_DIR <- file.path(RNA_DIR, "Integration")
JNT_DIR <- file.path(BASE, "analysis_results/Joint/Oral")
COR_DIR <- file.path(JNT_DIR, "GeneGenus_Cor")
OUT_DIR <- file.path(JNT_DIR, "Network")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ────────────────────────────────────────────────────────────
# 1. 读取数据
# ────────────────────────────────────────────────────────────

message("读取相关结果与表达数据...")

# 显著对（p<0.05）
edges_raw <- read.csv(
  file.path(COR_DIR, "sig_pairs_p05.csv"),
  stringsAsFactors = FALSE
)

# 挽救基因信息
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

# 挽救菌属信息
rescue_genus_df <- read.csv(
  file.path(JNT_DIR, "rescue_genus_list.csv"),
  stringsAsFactors = FALSE
) %>%
  filter(Genus != "Others")

# VST 表达矩阵 & metadata
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

# 属水平丰度
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
meta_rna <- meta_rna %>%
  mutate(s16_sample = rna_to_16s(rna_sample))

cols_oral <- colnames(genus_all)[str_ends(colnames(genus_all), "\\.O")]
meta_matched <- meta_rna %>%
  filter(s16_sample %in% cols_oral) %>%
  arrange(Group)

vst_matched <- vst_mat[, meta_matched$rna_sample, drop = FALSE]

# 将 vst_matched 行名从 ENSMUSG ID 转为 gene_name（与 edges_raw 的 Gene 列一致）
id_to_name <- setNames(rescue_genes$gene_name, rescue_genes$gene_id)
new_rownames <- id_to_name[rownames(vst_matched)]
# 只替换能找到 gene_name 的行，其余保留原 ID
rownames(vst_matched) <- ifelse(
  is.na(new_rownames) | new_rownames == "",
  rownames(vst_matched),
  new_rownames
)

genus_site <- genus_all[, meta_matched$s16_sample, drop = FALSE]

clr_transform <- function(mat) {
  m <- as.matrix(mat) + 1e-6
  sweep(log(m), 2, colMeans(log(m)), "-")
}
genus_clr <- clr_transform(genus_site)

# ────────────────────────────────────────────────────────────
# 2. 计算节点属性（三组均值 → 方向 & 幅度）
# ────────────────────────────────────────────────────────────

grp3 <- c("WT_Control", "DB_Control", "DB_Treated")

get_grp_means <- function(mat, meta, groups = grp3) {
  sapply(groups, function(g) {
    cols <- meta$rna_sample[meta$Group == g]
    rowMeans(mat[, cols, drop = FALSE])
  })
}
get_grp_means_16s <- function(mat, meta, groups = grp3) {
  sapply(groups, function(g) {
    cols <- meta$s16_sample[meta$Group == g]
    rowMeans(mat[, cols, drop = FALSE])
  })
}

# 基因：VST 均值（Z-score）
gene_means_raw <- get_grp_means(vst_matched, meta_matched)
gene_means_z <- t(scale(t(gene_means_raw))) # 行 Z-score

# 菌属：CLR 均值
genus_means <- get_grp_means_16s(genus_clr, meta_matched)

# 节点的方向信息（disease_dir）
# 基因：disease_dir=1 → DB中上调（挽救后下调）; -1 → DB中下调（挽救后上调）
gene_node_info <- rescue_genes %>%
  dplyr::select(gene_name, cluster) %>%
  distinct() %>%
  left_join(
    cluster_sum %>% dplyr::select(cluster, disease_dir),
    by = "cluster"
  ) %>%
  # 行名已改为 gene_name，直接过滤在网络边中出现的基因即可
  filter(gene_name %in% edges_raw$Gene)

# 菌属：disease_dir = sign(delta_disease)
genus_node_info <- rescue_genus_df %>%
  mutate(disease_dir = sign(delta_disease)) %>%
  dplyr::select(Genus, disease_dir, rescue_score)

# ────────────────────────────────────────────────────────────
# 3. 构建 igraph 对象
# ────────────────────────────────────────────────────────────

message("构建网络...")

# 节点表
nodes_gene <- data.frame(
  name = gene_node_info$gene_name,
  node_type = "Gene",
  disease_dir = gene_node_info$disease_dir,
  rescue_score = NA_real_,
  stringsAsFactors = FALSE
)

nodes_genus <- data.frame(
  name = genus_node_info$Genus,
  node_type = "Genus",
  disease_dir = genus_node_info$disease_dir,
  rescue_score = genus_node_info$rescue_score,
  stringsAsFactors = FALSE
)

nodes <- bind_rows(nodes_gene, nodes_genus) %>%
  distinct(name, .keep_all = TRUE)

# 边表（只保留节点表中存在的节点）
edges <- edges_raw %>%
  filter(Gene %in% nodes$name, Genus %in% nodes$name) %>%
  dplyr::select(from = Gene, to = Genus, r, p_raw, p_adj) %>%
  mutate(
    edge_sign = ifelse(r > 0, "正相关", "负相关"),
    abs_r = abs(r),
    # 协同挽救判断：
    # 当基因disease_dir 与菌属disease_dir 相同时，正相关 = 协同挽救
    # 当方向相反时，负相关 = 协同挽救
    gene_ddir = nodes$disease_dir[match(from, nodes$name)],
    genus_ddir = nodes$disease_dir[match(to, nodes$name)],
    is_synergistic = (gene_ddir == genus_ddir & r > 0) |
      (gene_ddir != genus_ddir & r < 0)
  )

# 构建 tidygraph
g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

# 添加度数
g <- g %>%
  mutate(degree = centrality_degree())

message(sprintf("  节点: %d  |  边: %d", vcount(g), ecount(g)))

# ────────────────────────────────────────────────────────────
# 4. 主网络图
# ────────────────────────────────────────────────────────────

message("绘制主网络图...")

# 节点颜色方案：
# disease_dir=1  (DB中上调，挽救后下调) → 红色系
# disease_dir=-1 (DB中下调，挽救后上调) → 蓝色系
# Gene vs Genus用形状区分

dir_colors <- c(
  "Gene_up" = "#B2182B", # 基因在DB中上调
  "Gene_down" = "#2166AC", # 基因在DB中下调
  "Genus_up" = "#F4A582", # 菌属在DB中上调
  "Genus_down" = "#92C5DE" # 菌属在DB中下调
)

g <- g %>%
  mutate(
    color_key = paste0(node_type, "_", ifelse(disease_dir == 1, "up", "down"))
  )

p_net <- ggraph(g, layout = "stress") +
  # 边：宽度=|r|，颜色=正/负，虚线=非协同
  geom_edge_link(
    aes(
      width = abs_r,
      color = edge_sign,
      linetype = ifelse(is_synergistic, "solid", "dashed"),
      alpha = abs_r
    ),
    show.legend = TRUE
  ) +
  scale_edge_width_continuous(range = c(0.5, 3.5), name = "|r|") +
  scale_edge_color_manual(
    values = c("正相关" = "#B2182B", "负相关" = "#2166AC"),
    name = "相关方向"
  ) +
  scale_edge_alpha_continuous(range = c(0.5, 1.0), guide = "none") +
  scale_edge_linetype_identity() +

  # 节点：形状区分基因/菌属；颜色区分DB中的表达方向
  geom_node_point(
    aes(shape = node_type, color = color_key, size = degree + 1),
    stroke = 1.2
  ) +
  scale_shape_manual(values = c("Gene" = 15, "Genus" = 19), name = "节点类型") +
  scale_color_manual(
    values = dir_colors,
    labels = c(
      "Gene_up" = "基因：DB(+)（挽救后(-)）",
      "Gene_down" = "基因：DB(-)（挽救后(+)）",
      "Genus_up" = "菌属：DB(+)（挽救后(-)）",
      "Genus_down" = "菌属：DB(-)（挽救后(+)）"
    ),
    name = "变化方向"
  ) +
  scale_size_continuous(range = c(4, 10), guide = "none") +

  # 节点标签
  geom_node_label(
    aes(label = name),
    size = 3.2,
    repel = TRUE,
    label.padding = unit(0.15, "lines"),
    label.size = 0.2,
    fill = alpha("white", 0.75),
    max.overlaps = 30
  ) +

  labs(
    title = "Oral  挽救基因 — 挽救菌属  协同变化网络",
    subtitle = paste0(
      "边：实线 = 协同挽救方向，虚线 = 拮抗方向\n",
      "节点颜色：DB中的表达/丰度变化方向（挽救后均回调至WT）\n",
      "节点大小：网络连接度   边宽：|Spearman r|   阈值：p < 0.05"
    )
  ) +
  theme_graph(base_family = "", base_size = 11) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey40")
  )

ggsave(
  file.path(OUT_DIR, "network_main.pdf"),
  p_net,
  width = 12,
  height = 9,
  device = "pdf"
)
ggsave(
  file.path(OUT_DIR, "network_main.png"),
  p_net,
  width = 12,
  height = 9,
  dpi = 200
)
message("  主网络图已保存")

# ────────────────────────────────────────────────────────────
# 5. 二部图布局（基因左，菌属右，更清晰地展示连接关系）
# ────────────────────────────────────────────────────────────

message("绘制二部图...")

# 手动指定二部图坐标
gene_names <- nodes$name[nodes$node_type == "Gene"]
genus_names <- nodes$name[nodes$node_type == "Genus"]

n_gene_net <- length(gene_names)
n_genus_net <- length(genus_names)

# 按连接数排序，使布局更美观
gene_deg <- degree(g)[gene_names]
genus_deg <- degree(g)[genus_names]
gene_names <- gene_names[order(gene_deg, decreasing = TRUE)]
genus_names <- genus_names[order(genus_deg, decreasing = TRUE)]

coords <- data.frame(
  name = c(gene_names, genus_names),
  x = c(rep(0, n_gene_net), rep(3, n_genus_net)),
  y = c(
    seq(n_gene_net, 1, length.out = n_gene_net),
    seq(n_genus_net, 1, length.out = n_genus_net)
  )
)

# 用 manual layout
layout_bipart <- create_layout(
  g,
  layout = "manual",
  x = coords$x[match(V(g)$name, coords$name)],
  y = coords$y[match(V(g)$name, coords$name)]
)

p_bipart <- ggraph(layout_bipart) +
  geom_edge_arc(
    aes(width = abs_r, color = edge_sign, alpha = abs_r),
    strength = 0.15,
    show.legend = TRUE
  ) +
  scale_edge_width_continuous(range = c(0.4, 3), name = "|r|") +
  scale_edge_color_manual(
    values = c("正相关" = "#B2182B", "负相关" = "#2166AC"),
    name = "相关方向"
  ) +
  scale_edge_alpha_continuous(range = c(0.4, 0.95), guide = "none") +

  geom_node_point(
    aes(shape = node_type, color = color_key, size = degree + 1)
  ) +
  scale_shape_manual(values = c("Gene" = 15, "Genus" = 19), name = "类型") +
  scale_color_manual(values = dir_colors, name = "DB中方向") +
  scale_size_continuous(range = c(5, 12), guide = "none") +

  # 基因标签（左侧，右对齐）
  geom_node_text(
    data = ~ filter(.x, node_type == "Gene"),
    aes(label = name),
    hjust = 1.15,
    size = 3.5,
    fontface = "plain"
  ) +
  # 菌属标签（右侧，左对齐）
  geom_node_text(
    data = ~ filter(.x, node_type == "Genus"),
    aes(label = name),
    hjust = -0.15,
    size = 3.5,
    fontface = "plain"
  ) +

  # 左右标注
  annotate(
    "text",
    x = 0,
    y = max(coords$y) + 0.5,
    label = "Gene",
    size = 5,
    fontface = "plain",
    hjust = 0.5,
    color = "grey30"
  ) +
  annotate(
    "text",
    x = 3,
    y = max(coords$y) + 0.5,
    label = "Genus",
    size = 5,
    fontface = "plain",
    hjust = 0.5,
    color = "grey30"
  ) +

  xlim(-1.5, 5) +
  labs(
    title = "Oral  挽救基因 — 挽救菌属  二部网络图",
    subtitle = "边色：正相关（红）/ 负相关（蓝）  |  边宽：|Spearman r|  |  p < 0.05"
  ) +
  theme_graph(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(
  file.path(OUT_DIR, "network_bipartite.pdf"),
  p_bipart,
  width = 13,
  height = max(8, n_gene_net * 0.45 + 3),
  device = "pdf"
)
ggsave(
  file.path(OUT_DIR, "network_bipartite.png"),
  p_bipart,
  width = 13,
  height = max(8, n_gene_net * 0.45 + 3),
  dpi = 200
)
message("  二部图已保存")

# ────────────────────────────────────────────────────────────
# 6. 附图：hub节点三组趋势折线图
#    （连接度最高的基因 + 菌属各取top N）
# ────────────────────────────────────────────────────────────

message("绘制 hub 节点趋势图...")

grp_order <- c("WT_Control", "DB_Control", "DB_Treated")
grp_colors <- c(
  "WT_Control" = "#4DBBD5",
  "DB_Control" = "#E64B35",
  "DB_Treated" = "#F39B7F"
)
grp_labels <- c(
  "WT_Control" = "WT",
  "DB_Control" = "DB",
  "DB_Treated" = "DB+Sema"
)

node_degrees <- degree(g)

# Top hub 基因（最多5个）
top_hub_genes <- names(sort(
  node_degrees[names(node_degrees) %in% gene_names],
  decreasing = TRUE
))[1:min(5, n_gene_net)]

# Top hub 菌属（最多5个）
top_hub_genus <- names(sort(
  node_degrees[names(node_degrees) %in% genus_names],
  decreasing = TRUE
))[1:min(5, n_genus_net)]

make_gene_trend <- function(gene) {
  if (!gene %in% rownames(vst_matched)) {
    return(NULL)
  }
  df <- data.frame(
    sample = colnames(vst_matched),
    value = as.numeric(vst_matched[gene, ]),
    Group = meta_matched$Group,
    stringsAsFactors = FALSE
  ) %>%
    filter(Group %in% grp_order) %>%
    mutate(Group = factor(Group, levels = grp_order))

  summ <- df %>%
    group_by(Group) %>%
    summarise(
      mean_v = mean(value),
      se_v = sd(value) / sqrt(n()),
      .groups = "drop"
    )

  # 获取disease_dir
  ddir <- rescue_genes$cluster[rescue_genes$gene_name == gene][1]
  ddir_info <- cluster_sum$disease_dir[cluster_sum$cluster == ddir][1]
  subtitle_txt <- ifelse(
    ddir_info == 1,
    "DB(+) → Sema后(-)（挽救）",
    "DB(-) → Sema后(+)（挽救）"
  )

  ggplot(summ, aes(Group, mean_v, group = 1)) +
    geom_ribbon(
      aes(ymin = mean_v - se_v, ymax = mean_v + se_v),
      fill = "grey85",
      alpha = 0.6
    ) +
    geom_line(linewidth = 1.2, color = "grey40") +
    geom_point(aes(color = Group), size = 4) +
    scale_color_manual(values = grp_colors, guide = "none") +
    scale_x_discrete(labels = grp_labels) +
    labs(
      title = paste0("Gene: ", gene),
      subtitle = subtitle_txt,
      x = NULL,
      y = "VST（均值±SE）"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold.italic"),
      plot.subtitle = element_text(size = 8, color = "grey50"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
}

make_genus_trend <- function(genus) {
  if (!genus %in% rownames(genus_clr)) {
    return(NULL)
  }
  df <- data.frame(
    sample = colnames(genus_clr),
    value = as.numeric(genus_clr[genus, ]),
    Group = meta_matched$Group,
    stringsAsFactors = FALSE
  ) %>%
    filter(Group %in% grp_order) %>%
    mutate(Group = factor(Group, levels = grp_order))

  summ <- df %>%
    group_by(Group) %>%
    summarise(
      mean_v = mean(value),
      se_v = sd(value) / sqrt(n()),
      .groups = "drop"
    )

  ddir_info <- sign(rescue_genus_df$delta_disease[
    rescue_genus_df$Genus == genus
  ][1])
  subtitle_txt <- ifelse(
    ddir_info == 1,
    "DB(+) → Sema后(-)（挽救）",
    "DB(-) → Sema后(+)（挽救）"
  )

  ggplot(summ, aes(Group, mean_v, group = 1)) +
    geom_ribbon(
      aes(ymin = mean_v - se_v, ymax = mean_v + se_v),
      fill = "grey85",
      alpha = 0.6
    ) +
    geom_line(linewidth = 1.2, color = "grey40") +
    geom_point(aes(color = Group), size = 4) +
    scale_color_manual(values = grp_colors, guide = "none") +
    scale_x_discrete(labels = grp_labels) +
    labs(
      title = paste0("Genus: ", genus),
      subtitle = subtitle_txt,
      x = NULL,
      y = "CLR（均值±SE）"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold.italic"),
      plot.subtitle = element_text(size = 8, color = "grey50"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
}

# 生成趋势图列表
gene_trends <- Filter(Negate(is.null), lapply(top_hub_genes, make_gene_trend))
genus_trends <- Filter(Negate(is.null), lapply(top_hub_genus, make_genus_trend))
all_trends <- c(gene_trends, genus_trends)

if (length(all_trends) > 0) {
  p_trends <- wrap_plots(all_trends, ncol = min(5, length(all_trends))) +
    plot_annotation(
      title = "Oral  Hub 节点 WT→DB→DB-Treated 表达趋势",
      subtitle = "灰色带：±SE   点：各组均值",
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
      )
    )

  n_row_t <- ceiling(length(all_trends) / min(5, length(all_trends)))
  ggsave(
    file.path(OUT_DIR, "hub_node_trends.pdf"),
    p_trends,
    width = min(5, length(all_trends)) * 3.5,
    height = n_row_t * 3.8,
    device = "pdf"
  )
  ggsave(
    file.path(OUT_DIR, "hub_node_trends.png"),
    p_trends,
    width = min(5, length(all_trends)) * 3.5,
    height = n_row_t * 3.8,
    dpi = 200
  )
  message("  趋势图已保存")
}

# ────────────────────────────────────────────────────────────
# 7. 汇总
# ────────────────────────────────────────────────────────────

cat(sprintf("\n%s\n   网络可视化完成\n%s\n", strrep("=", 50), strrep("=", 50)))
cat(sprintf(
  "   节点: %d（基因 %d，菌属 %d）\n",
  vcount(g),
  n_gene_net,
  n_genus_net
))
cat(sprintf("   边:   %d\n", ecount(g)))
cat("\n   输出文件:\n")
cat("   network_main.pdf        — 力导向主网络图\n")
cat("   network_bipartite.pdf   — 二部图（基因左 | 菌属右）\n")
cat("   hub_node_trends.pdf     — Hub节点 WT→DB→DB-Treated 趋势\n")
cat(sprintf("\n   输出目录: %s\n", OUT_DIR))
