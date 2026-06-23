# ==============================================================================
# KEGG 通路分类注释（宏转录组）- 环形树状图
# 数据来源：图片手动提取
# ==============================================================================

library(ggplot2)
library(ggraph)
library(tidygraph)
library(stringr)

# ------------------------------------------------------------------------------
# 1. 数据（从图片提取）
# ------------------------------------------------------------------------------

data <- data.frame(
  from = c(
    # Metabolism
    rep("Metabolism", 5),
    # Genetic Information Processing
    rep("Genetic Information Processing", 1),
    # Environmental Information Processing
    rep("Environmental Information Processing", 2),
    # Cellular Processes
    rep("Cellular Processes", 3),
    # Organismal Systems
    rep("Organismal Systems", 7),
    # Human Diseases
    rep("Human Diseases", 11)
  ),
  to = c(
    # Metabolism
    "Global and overview maps",
    "Amino acid metabolism",
    "Energy metabolism",
    "Carbohydrate metabolism",
    "Lipid metabolism",
    # Genetic Information Processing
    "Folding, sorting and degradation",
    # Environmental Information Processing
    "Signaling molecules and interaction",
    "Signal transduction",
    # Cellular Processes
    "Transport and catabolism",
    "Cell growth and death",
    "Cellular community - eukaryotes",
    # Organismal Systems
    "Immune system",
    "Development and regeneration",
    "Endocrine system",
    "Nervous system",
    "Environmental adaptation",
    "Excretory system",
    "Aging",
    # Human Diseases
    "Immune disease",
    "Infectious disease: viral",
    "Infectious disease: bacterial",
    "Infectious disease: parasitic",
    "Cancer: overview",
    "Endocrine and metabolic disease",
    "Cardiovascular disease",
    "Neurodegenerative disease",
    "Cancer: specific types",
    "Substance dependence",
    "Drug resistance: antineoplastic"
  ),
  Count = c(
    # Metabolism
    5,
    2,
    1,
    1,
    1,
    # Genetic Information Processing
    1,
    # Environmental Information Processing
    35,
    25,
    # Cellular Processes
    10,
    5,
    2,
    # Organismal Systems
    48,
    6,
    6,
    5,
    3,
    1,
    1,
    # Human Diseases
    25,
    25,
    21,
    20,
    18,
    13,
    10,
    5,
    3,
    2,
    1
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. 配色（参照原图）
# ------------------------------------------------------------------------------

cats <- c(
  "Metabolism",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Cellular Processes",
  "Organismal Systems",
  "Human Diseases"
)

cols <- c(
  "Metabolism" = "#8DC63F",
  "Genetic Information Processing" = "#9170C8",
  "Environmental Information Processing" = "#7FA0C8",
  "Cellular Processes" = "#F08078",
  "Organismal Systems" = "#E89050",
  "Human Diseases" = "#2ECC78"
)

cols_light <- c(
  "Metabolism" = "#BBDC88",
  "Genetic Information Processing" = "#C0AADE",
  "Environmental Information Processing" = "#AABEDE",
  "Cellular Processes" = "#F8B0A8",
  "Organismal Systems" = "#F0B888",
  "Human Diseases" = "#88DEB0"
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
# 4. 提取叶节点位置
# ------------------------------------------------------------------------------

layout_df <- create_layout(gr, layout = "dendrogram", circular = TRUE)

leaf_pos <- layout_df[layout_df$level == 2, c("name", "x", "y")]
leaf_pos <- merge(leaf_pos, data, by.x = "name", by.y = "to", all.x = FALSE)

leaf_pos$r <- sqrt(leaf_pos$x^2 + leaf_pos$y^2)
leaf_pos$ux <- leaf_pos$x / leaf_pos$r
leaf_pos$uy <- leaf_pos$y / leaf_pos$r

bar_max <- 1.2
leaf_pos$bh <- leaf_pos$Count / max(leaf_pos$Count) * bar_max
leaf_pos$xend <- leaf_pos$x + leaf_pos$ux * leaf_pos$bh
leaf_pos$yend <- leaf_pos$y + leaf_pos$uy * leaf_pos$bh

gap <- 0.08
leaf_pos$tx <- leaf_pos$xend + leaf_pos$ux * gap
leaf_pos$ty <- leaf_pos$yend + leaf_pos$uy * gap

theta <- atan2(leaf_pos$y, leaf_pos$x) * 180 / pi
leaf_pos$text_angle <- ifelse(abs(theta) > 90, theta + 180, theta)
leaf_pos$text_hjust <- ifelse(abs(theta) > 90, 1, 0)

# ------------------------------------------------------------------------------
# 5. 绘图
# ------------------------------------------------------------------------------

p <- ggraph(layout_df) +

  geom_edge_diagonal(
    aes(colour = colour),
    alpha = 0.7,
    width = 1,
    show.legend = FALSE
  ) +
  geom_segment(
    data = leaf_pos,
    aes(x = x, y = y, xend = xend, yend = yend, colour = from),
    linewidth = 6,
    lineend = "butt",
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(filter = level == 1, colour = name),
    size = 6,
    show.legend = FALSE
  ) +

  geom_node_text(
    aes(filter = level == 1, label = label_wrapped, colour = name),
    lineheight = 1,
    repel = TRUE,
    size = 3.5,
    show.legend = FALSE
  ) +

  geom_text(
    data = leaf_pos,
    aes(x = tx, y = ty, label = name, angle = text_angle, hjust = text_hjust),
    size = 5,
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

ggsave(file.path(out_dir, "kegg_radial_rna.pdf"), p, width = 9, height = 9)
ggsave(
  file.path(out_dir, "kegg_radial_rna.png"),
  p,
  width = 12,
  height = 12,
  dpi = 150
)

cat("已保存至:", out_dir, "\n")
