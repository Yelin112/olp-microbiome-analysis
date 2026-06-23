# ==============================================================================
# KEGG 通路分类注释 - 环形树状图（单坐标系版，内外完全对齐）
# 数据来源：图片手动提取
# ==============================================================================

library(ggplot2)
library(ggraph)
library(tidygraph)
library(stringr)

# ------------------------------------------------------------------------------
# 1. 数据
# ------------------------------------------------------------------------------

data <- data.frame(
  from = c(
    rep("Cellular Processes", 5),
    rep("Environmental Information Processing", 2),
    rep("Genetic Information Processing", 4),
    rep("Human Diseases", 11),
    rep("Metabolism", 12),
    rep("Organismal Systems", 8)
  ),
  to = c(
    # Cellular Processes
    "Cell growth and death",
    "Cell motility",
    "Cellular community - eukaryotes",
    "Cellular community - prokaryotes",
    "Transport and catabolism",
    # Environmental Information Processing
    "Membrane transport",
    "Signal transduction",
    # Genetic Information Processing
    "Folding, sorting and degradation",
    "Replication and repair",
    "Transcription",
    "Translation",
    # Human Diseases
    "Cancer: overview",
    "Cancer: specific types",
    "Cardiovascular disease",
    "Drug resistance: antimicrobial",
    "Drug resistance: antineoplastic",
    "Endocrine and metabolic disease",
    "Immune disease",
    "Infectious disease: bacterial",
    "Infectious disease: parasitic",
    "Infectious disease: viral",
    "Neurodegenerative disease",
    # Metabolism
    "Amino acid metabolism",
    "Biosynthesis of other secondary metabolites",
    "Carbohydrate metabolism",
    "Energy metabolism",
    "Global and overview maps",
    "Glycan biosynthesis and metabolism",
    "Lipid metabolism",
    "Metabolism of cofactors and vitamins",
    "Metabolism of other amino acids",
    "Metabolism of terpenoids and polyketides",
    "Nucleotide metabolism",
    "Xenobiotics biodegradation and metabolism",
    # Organismal Systems
    "Aging",
    "Circulatory system",
    "Digestive system",
    "Endocrine system",
    "Environmental adaptation",
    "Excretory system",
    "Immune system",
    "Nervous system"
  ),
  Count = c(
    91,
    333,
    1,
    255,
    31,
    382,
    460,
    297,
    80,
    120,
    903,
    65,
    11,
    127,
    83,
    23,
    103,
    1,
    286,
    2,
    10,
    66,
    377,
    70,
    702,
    536,
    1589,
    108,
    134,
    228,
    155,
    57,
    148,
    64,
    119,
    3,
    18,
    90,
    219,
    13,
    140,
    3
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. 配色
# ------------------------------------------------------------------------------

cats <- c(
  "Cellular Processes",
  "Environmental Information Processing",
  "Genetic Information Processing",
  "Human Diseases",
  "Metabolism",
  "Organismal Systems"
)

cols <- c(
  "Cellular Processes" = "#4472C4",
  "Environmental Information Processing" = "#70AD47",
  "Genetic Information Processing" = "#C9956C",
  "Human Diseases" = "#8080A0",
  "Metabolism" = "#17B8CE",
  "Organismal Systems" = "#F5821F"
)

cols_light <- c(
  "Cellular Processes" = "#8AAED8",
  "Environmental Information Processing" = "#A8CC88",
  "Genetic Information Processing" = "#DEC0A0",
  "Human Diseases" = "#B8B8CC",
  "Metabolism" = "#80D8E4",
  "Organismal Systems" = "#F8B878"
)

# ------------------------------------------------------------------------------
# 3. 构建图结构
# ------------------------------------------------------------------------------

edges_graph <- rbind(
  data[, c("from", "to", "Count")],
  data.frame(from = "pathway", to = cats, Count = 0, stringsAsFactors = FALSE)
)
edges_graph$colour <- ifelse(
  edges_graph$from == "pathway",
  edges_graph$to,
  edges_graph$from
)

all_nodes <- unique(c(data$from, data$to))
vertices <- data.frame(
  name = c("pathway", all_nodes),
  level = c(0L, ifelse(all_nodes %in% cats, 1L, 2L)),
  stringsAsFactors = FALSE
)
vertices$label_wrapped <- str_wrap(vertices$name, 14)

gr <- tbl_graph(nodes = vertices, edges = edges_graph, directed = TRUE)

# ------------------------------------------------------------------------------
# 4. 提取叶节点位置（关键：在 ggraph 坐标系内画条形，保证对齐）
# ------------------------------------------------------------------------------

layout_df <- create_layout(gr, layout = "dendrogram", circular = TRUE)

leaf_pos <- layout_df[layout_df$level == 2, c("name", "x", "y")]
leaf_pos <- merge(leaf_pos, data, by.x = "name", by.y = "to", all.x = FALSE)

# 叶节点到圆心的单位向量（用于向外延伸条形）
leaf_pos$r <- sqrt(leaf_pos$x^2 + leaf_pos$y^2)
leaf_pos$ux <- leaf_pos$x / leaf_pos$r
leaf_pos$uy <- leaf_pos$y / leaf_pos$r

# 条形高度：线性缩放，最长条形 = 1.5 坐标单位
bar_max <- 1.2
leaf_pos$bh <- leaf_pos$Count / max(leaf_pos$Count) * bar_max
leaf_pos$xend <- leaf_pos$x + leaf_pos$ux * leaf_pos$bh
leaf_pos$yend <- leaf_pos$y + leaf_pos$uy * leaf_pos$bh

# 标签位置（条形末端再往外一点）
gap <- 0.08
leaf_pos$tx <- leaf_pos$xend + leaf_pos$ux * gap
leaf_pos$ty <- leaf_pos$yend + leaf_pos$uy * gap

# 标签角度（左半圆旋转 180° 保证可读）
theta <- atan2(leaf_pos$y, leaf_pos$x) * 180 / pi
leaf_pos$text_angle <- ifelse(abs(theta) > 90, theta + 180, theta)
leaf_pos$text_hjust <- ifelse(abs(theta) > 90, 1, 0)

# ------------------------------------------------------------------------------
# 5. 绘图（单坐标系，内外完全对齐）
# ------------------------------------------------------------------------------

p <- ggraph(layout_df) +

  # 树状边
  geom_edge_diagonal(aes(colour = colour), alpha = 0.7, show.legend = FALSE) +

  # 条形（从叶节点向外延伸）
  geom_segment(
    data = leaf_pos,
    aes(x = x, y = y, xend = xend, yend = yend, colour = from),
    linewidth = 7,
    lineend = "butt",
    show.legend = TRUE
  ) +

  # 分类节点
  geom_node_point(
    aes(filter = level == 1, colour = name),
    size = 5,
    show.legend = FALSE
  ) +

  # 分类标签
  geom_node_text(
    aes(filter = level == 1, label = label_wrapped, colour = name),
    lineheight = 1,
    repel = TRUE,
    size = 4,
    show.legend = FALSE
  ) +

  # 通路名称标签
  geom_text(
    data = leaf_pos,
    aes(x = tx, y = ty, label = name, angle = text_angle, hjust = text_hjust),
    size = 4,
    lineheight = 1,
    show.legend = FALSE
  ) +

  coord_fixed() +
  scale_edge_colour_manual(values = cols_light) +
  scale_colour_manual(name = NULL, values = cols) +
  guides(colour = guide_legend(override.aes = list(linewidth = 4))) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.6, "cm")
  )

# ------------------------------------------------------------------------------
# 6. 保存
# ------------------------------------------------------------------------------

out_dir <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/富集重画"

ggsave(file.path(out_dir, "kegg_radial.pdf"), p, width = 9, height = 9)
ggsave(
  file.path(out_dir, "kegg_radial.png"),
  p,
  width = 12,
  height = 12,
  dpi = 150
)

cat("已保存至:", out_dir, "\n")
