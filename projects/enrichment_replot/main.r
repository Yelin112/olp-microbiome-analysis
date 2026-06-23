library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# 1. 手动录入提取的数据 (为了演示，每个大类选取部分代表性子类)
# 条目名称
Description <- c(
  "Cell growth and death",
  "Cell motility",
  "Cellular community - eukaryotes",
  "Cellular community - prokaryotes",
  "Transport and catabolism",
  "Membrane transport",
  "Signal transduction",
  "Signaling molecules and interaction",
  "Folding, sorting and degradation",
  "Information processing in viruses",
  "Replication and repair",
  "Transcription",
  "Translation",
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
  "Substance dependence",
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
  "Aging",
  "Circulatory system",
  "Development and regeneration",
  "Digestive system",
  "Endocrine system",
  "Environmental adaptation",
  "Excretory system",
  "Immune system",
  "Nervous system",
  "Sensory system"
)

# 横杆右侧数字（从图片提取）
GeneCount <- c(
  26631,
  19124,
  460,
  71391,
  13806,
  107623,
  60986,
  482,
  44479,
  584,
  83625,
  7528,
  67325,
  15079,
  2472,
  10235,
  33026,
  8287,
  10532,
  1573,
  25284,
  1853,
  2182,
  3864,
  81,
  125002,
  31328,
  158223,
  77642,
  516325,
  85411,
  44111,
  106665,
  41112,
  24804,
  74670,
  17466,
  9344,
  296,
  249,
  5633,
  15842,
  7738,
  519,
  3172,
  2156,
  13
)

# 类别（根据颜色图例）
ONTOLOGY <- c(
  rep("Cellular Processes", 5),
  rep("Environmental Information Processing", 2),
  rep("Genetic Information Processing", 4),
  rep("Human Diseases", 7),
  rep("Metabolism", 9),
  rep("Organismal Systems", 20)
)

# 构建data.frame
data <- data.frame(Description, GeneCount, ONTOLOGY)

# 查看前几行
head(data)


#-----------------------------
# 2. 构建数据框
#-----------------------------
data <- data.frame(
  Description = Description,
  GeneCount = GeneCount,
  ONTOLOGY = ONTOLOGY,
  stringsAsFactors = FALSE
)

#-----------------------------
# 3. 数据整理
#-----------------------------
data <- data %>%
  group_by(ONTOLOGY) %>%
  arrange(GeneCount, .by_group = TRUE) %>% # 小到大排序
  mutate(
    Description_wrap = str_wrap(Description, width = 35),
    Description_wrap = factor(Description_wrap, levels = Description_wrap)
  ) %>%
  filter(GeneCount > 500) %>%
  ungroup()

#-----------------------------
# 3. 绘图
#-----------------------------
ggplot(data, aes(x = GeneCount, y = Description_wrap, fill = ONTOLOGY)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(
    aes(x = GeneCount * 0.02, label = Description_wrap),
    hjust = 0,
    size = 3,
    lineheight = 0.9
  ) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c(
      "Cellular Processes" = "#1cc7d0",
      "Environmental Information Processing" = "#2dde98",
      "Genetic Information Processing" = "#ff6c5f",
      "Human Diseases" = "#3369e7",
      "Metabolism" = "#f4b400",
      "Organismal Systems" = "#9c27b0"
    )
  ) +
  labs(x = "Gene Count", y = NULL) +
  expand_limits(x = max(data$GeneCount) * 1.15) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10)
  )


library(ggplot2)
library(dplyr)
library(stringr)

#==============================
# 1. 原始数据（使用你提供的修正版）
#==============================
data <- data.frame(
  Description = c(
    "Cell growth and death",
    "Cell motility",
    "Cellular community - eukaryotes",
    "Cellular community - prokaryotes",
    "Transport and catabolism",

    "Membrane transport",
    "Signal transduction",
    "Signaling molecules and interaction",

    "Folding, sorting and degradation",
    "Information processing in viruses",
    "Replication and repair",
    "Transcription",
    "Translation",

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
    "Substance dependence",

    "Global and overview maps",
    "Carbohydrate metabolism",
    "Amino acid metabolism",
    "Metabolism of cofactors and vitamins",
    "Glycan biosynthesis and metabolism",
    "Energy metabolism",
    "Nucleotide metabolism",
    "Lipid metabolism",
    "Metabolism of other amino acids",
    "Biosynthesis of other secondary metabolites",
    "Metabolism of terpenoids and polyketides",
    "Xenobiotics biodegradation and metabolism",

    "Endocrine system",
    "Aging",
    "Environmental adaptation",
    "Digestive system",
    "Immune system",
    "Nervous system",
    "Excretory system",
    "Circulatory system",
    "Development and regeneration",
    "Sensory system"
  ),

  GeneCount = c(
    26631,
    19124,
    460,
    71391,
    13806,
    107623,
    60986,
    482,
    44479,
    584,
    83625,
    7528,
    67325,
    15079,
    2472,
    10235,
    33026,
    8287,
    10532,
    1573,
    25284,
    1853,
    2182,
    3864,
    81,
    516325,
    158223,
    125002,
    106665,
    85411,
    77642,
    74670,
    44111,
    41112,
    31328,
    24804,
    17466,
    15842,
    9344,
    7738,
    5633,
    3172,
    2156,
    519,
    296,
    249,
    13
  ),

  ONTOLOGY = c(
    rep("Cellular Processes", 5),
    rep("Environmental Information Processing", 3),
    rep("Genetic Information Processing", 5),
    rep("Human Diseases", 12),
    rep("Metabolism", 12),
    rep("Organismal Systems", 10)
  ),

  stringsAsFactors = FALSE
)

#==============================
# 2. 按 ONTOLOGY 分组，选出每组前 9 个并排序
#==============================
data_top <- data %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(GeneCount), .by_group = TRUE) %>% # 从高到低排序
  slice_head(n = 9) %>% # 每组最多 9 条
  ungroup() %>%
  mutate(
    Description_wrap = str_wrap(Description, width = 35)
  )

#==============================
# 3. 固定 ONTOLOGY 顺序
#==============================
data_top$ONTOLOGY <- factor(
  data_top$ONTOLOGY,
  levels = c(
    "Cellular Processes",
    "Environmental Information Processing",
    "Genetic Information Processing",
    "Human Diseases",
    "Metabolism",
    "Organismal Systems"
  )
)

#==============================
# 4. 绘图
#==============================
ggplot(
  data_top,
  aes(x = GeneCount, y = reorder(Description_wrap, GeneCount), fill = ONTOLOGY)
) +
  geom_col(width = 0.75, show.legend = FALSE) +

  # 在 bar 上显示通路名称
  geom_text(
    aes(x = 0, label = Description_wrap),
    hjust = 0,
    nudge_x = max(data_top$GeneCount) * 0.01,
    size = 3,
    lineheight = 0.9
  ) +
  # 可选：显示柱子末端数值
  # geom_text(
  #   aes(x = GeneCount, label = format(GeneCount, big.mark = ",")),
  #   hjust = -0.1,
  #   size = 3
  # ) +
  facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c(
      "Cellular Processes" = "#3b5ba9",
      "Environmental Information Processing" = "#4daf4a",
      "Genetic Information Processing" = "#8c510a",
      "Human Diseases" = "#8073ac",
      "Metabolism" = "#00bfc4",
      "Organismal Systems" = "#f28e2b"
    )
  ) +

  labs(x = "Gene Count", y = NULL) +
  expand_limits(x = max(data_top$GeneCount) * 1.15) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10)
  )


# 第二个富集 ------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(patchwork)

#-----------------------------
# 1. 读取数据
#-----------------------------
df <- read_csv("富集重画//宏基因组KEGG富集通路分类注释.csv")

#-----------------------------
# 2. 数据整理
#-----------------------------
plot.data <- df %>%
  transmute(
    ID = ID,
    Description = Description,
    ONTOLOGY = level1,
    ReporterScore = ReporterScore,
    p.adjust = p.adjust
  ) %>%
  filter(
    !is.na(ONTOLOGY),
    !is.na(ReporterScore),
    !is.na(p.adjust),
    p.adjust > 0
  ) %>%
  mutate(
    sig = -log10(p.adjust),
    Change = ifelse(ReporterScore >= 0, "Up", "Down"),
    Description_wrap = str_wrap(Description, width = 35)
  ) %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(abs(ReporterScore)), .by_group = TRUE) %>% # 每个分面按绝对富集分数排序
  slice_head(n = 9) %>% # 每个分面最多9条
  ungroup()

# 固定分面顺序（可按需要修改）
plot.data$ONTOLOGY <- factor(
  plot.data$ONTOLOGY,
  levels = c(
    "Cellular Processes",
    "Environmental Information Processing",
    "Genetic Information Processing",
    "Human Diseases",
    "Metabolism",
    "Organismal Systems"
  )
)

# 全局最大值，用于统一不同分面的坐标范围
max_sig <- max(plot.data$sig, na.rm = TRUE)
max_rs <- max(abs(plot.data$ReporterScore), na.rm = TRUE)

# 颜色：按上下调显示
two_colors <- c(
  "Up" = "#E64B35",
  "Down" = "#4DBBD5"
)

#-----------------------------
# 3. 定义每个分面的左右拼图函数
#-----------------------------
make_pair_plot <- function(df_sub) {
  facet_name <- unique(df_sub$ONTOLOGY)

  df_sub <- df_sub %>%
    arrange(desc(abs(ReporterScore))) %>%
    mutate(
      Description_wrap = factor(
        Description_wrap,
        levels = rev(unique(Description_wrap))
      )
    )

  # 左图：-log10(p.adjust) 点线图
  p_left <- ggplot(df_sub, aes(x = sig, y = Description_wrap, group = 1)) +
    geom_line(color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = Change), size = 3) +
    scale_color_manual(values = two_colors) +
    scale_x_reverse(
      limits = c(max(df_sub$sig) * 1.05, 0),
      breaks = seq(0, ceiling(max(df_sub$sig) * 1.05), by = 1)
    ) +
    labs(x = expression(-log[10](p.adjust)), y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 9),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = margin(t = 5, r = 2, b = 5, l = 5)
    )

  # 右图：ReporterScore bar图
  p_right <- ggplot(
    df_sub,
    aes(x = ReporterScore, y = Description_wrap, fill = Change)
  ) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_text(
      aes(x = ReporterScore / 2, label = Description_wrap),
      size = 3,
      lineheight = 0.9
    ) +
    scale_fill_manual(values = two_colors) +
    labs(x = "ReporterScore", y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 9),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 2)
    )

  title_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = facet_name,
      fontface = "bold",
      size = 5
    ) +
    theme_void()

  title_plot /
    (p_left + p_right + plot_layout(widths = c(0.55, 1.45))) +
    plot_layout(heights = c(0.08, 1))
  # 给每个拼接图加标题
  # combined_plot <- p_left +
  #   p_right +
  #   plot_layout(widths = c(0.4, 1.3)) +
  #   plot_annotation(
  #     title = as.character(facet_name),
  #     theme = theme(
  #       plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  #     )
  #   )

  # return(combined_plot)
}

#-----------------------------
# 4. 对每个 ONTOLOGY 分别作图，再上下拼接
#-----------------------------
plot_list <- lapply(levels(plot.data$ONTOLOGY), function(x) {
  df_sub <- plot.data %>% filter(ONTOLOGY == x)
  if (nrow(df_sub) > 0) {
    make_pair_plot(df_sub)
  }
})


plot_list <- plot_list[!sapply(plot_list, is.null)]

final_plot <- wrap_plots(plotlist = plot_list, ncol = 2)
final_plot


# 宏转录组可视化 ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(patchwork)

#===========================
# 1. 读取宏转录组数据
#===========================
plot.data <- read_csv("富集重画//宏转录组KEGG富集通路分类注释.csv") %>%
  transmute(
    ID = ID,
    Description = Description,
    ONTOLOGY = level1,
    ReporterScore = ReporterScore,
    p.adjust = p.adjust
  ) %>%
  filter(
    !is.na(ONTOLOGY),
    !is.na(ReporterScore),
    !is.na(p.adjust),
    p.adjust > 0
  ) %>%
  mutate(
    sig = -log10(p.adjust),
    Change = ifelse(ReporterScore >= 0, "Up", "Down"),
    Description_wrap = str_wrap(Description, width = 35)
  ) %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(abs(ReporterScore)), .by_group = TRUE) %>% # 每个分面按绝对值排序
  slice_head(n = 9) %>% # 每个分面最多显示 9 条通路
  ungroup()

plot.data$ONTOLOGY <- factor(
  plot.data$ONTOLOGY,
  levels = unique(plot.data$ONTOLOGY)
)

max_sig <- max(plot.data$sig, na.rm = TRUE)
max_rs <- max(abs(plot.data$ReporterScore), na.rm = TRUE)

two_colors <- c("Up" = "#E64B35", "Down" = "#4DBBD5")

#===========================
# 2. 定义每个分面的左右拼图函数
#===========================
make_pair_plot <- function(df_sub) {
  facet_name <- unique(df_sub$ONTOLOGY)

  df_sub <- df_sub %>%
    arrange(desc(abs(ReporterScore))) %>%
    mutate(
      Description_wrap = factor(
        Description_wrap,
        levels = rev(unique(Description_wrap))
      ),
      label_x = ifelse(ReporterScore >= 0, -0.2, 0.2),
      hjust = ifelse(ReporterScore >= 0, 1, 0)
    )

  # 左图：-log10(p.adjust) 点线图
  p_left <- ggplot(df_sub, aes(x = sig, y = Description_wrap, group = 1)) +
    geom_line(color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = Change), size = 3) +
    scale_color_manual(values = two_colors) +
    scale_x_reverse(
      limits = c(max(df_sub$sig) * 1.05, 0),
      breaks = seq(0, ceiling(max(df_sub$sig) * 1.05), by = 1)
    ) +
    labs(x = expression(-log[10](p.adjust)), y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 9),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = margin(t = 5, r = 2, b = 5, l = 5)
    )

  # 右图：ReporterScore 柱状图
  p_right <- ggplot(
    df_sub,
    aes(x = ReporterScore, y = Description_wrap, fill = Change)
  ) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_text(
      aes(x = label_x, label = Description_wrap, hjust = hjust),
      size = 3,
      lineheight = 0.9
    ) +
    # geom_text(aes(x = ReporterScore, label = round(ReporterScore, 2),
    #               hjust = ifelse(ReporterScore >= 0, -0.15, 1.15)),
    #           size = 2.8) +
    scale_fill_manual(values = two_colors) +
    scale_x_continuous(limits = c(-max_rs * 1.25, max_rs * 1.25)) +
    labs(x = "ReporterScore", y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.text.x = element_text(size = 9),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 2)
    )

  # 添加分面标题
  title_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = facet_name,
      fontface = "bold",
      size = 5
    ) +
    theme_void()

  title_plot /
    (p_left + p_right + plot_layout(widths = c(0.55, 1.45))) +
    plot_layout(heights = c(0.08, 1))
}

#===========================
# 3. 对每个分面生成图并上下拼接
#===========================
plot_list <- lapply(levels(plot.data$ONTOLOGY), function(x) {
  df_sub <- plot.data %>% filter(ONTOLOGY == x)
  if (nrow(df_sub) > 0) {
    make_pair_plot(df_sub)
  }
})

plot_list <- plot_list[!sapply(plot_list, is.null)]

final_plot <- wrap_plots(plotlist = plot_list, ncol = 2)

#===========================
# 4. 显示最终图
#===========================
final_plot


# Red模块富集 ----------------------------------------------------------------

plot.data <- data.frame(
  Description = c(
    "T cell receptor signaling pathway",
    "Th17 cell differentiation",
    "Th1 and Th2 cell differentiation",
    "Cell adhesion molecules",
    "Chagas disease",
    "PD-L1 expression and PD-1 checkpoint pathway in cancer",
    "Antigen processing and presentation",
    "Viral protein interaction with cytokine and cytokine receptor",
    "Graft-versus-host disease",
    "Cytokine-cytokine receptor interaction",
    "Type I diabetes mellitus",
    "Inflammatory bowel disease",
    "Chemokine signaling pathway",
    "Allograft rejection",
    "Measles",
    "Primary immunodeficiency",
    "Autoimmune thyroid disease",
    "Rheumatoid arthritis",
    "Hematopoietic cell lineage",
    "Systemic lupus erythematosus",
    "NF-kappa B signaling pathway",
    "Malaria",
    "Intestinal immune network for IgA production",
    "Human immunodeficiency virus 1 infection",
    "Human T-cell leukemia virus 1 infection"
  ),
  RichFactor = c(
    0.165,
    0.145,
    0.155,
    0.102,
    0.125,
    0.123,
    0.127,
    0.104,
    0.190,
    0.052,
    0.169,
    0.109,
    0.058,
    0.079,
    0.066,
    0.071,
    0.069,
    0.059,
    0.057,
    0.048,
    0.055,
    0.105,
    0.066,
    0.044,
    0.042
  ),
  GeneNumber = c(
    18,
    15,
    15,
    15,
    12,
    12,
    10,
    12,
    6,
    15,
    6,
    6,
    12,
    9,
    9,
    9,
    9,
    9,
    9,
    12,
    12,
    6,
    8,
    12,
    12
  ),
  p.adjust = c(
    1.0e-5,
    1.0e-5,
    1.0e-5,
    1.1e-5,
    1.1e-5,
    1.1e-5,
    1.1e-5,
    1.1e-5,
    1.0e-5,
    1.2e-5,
    1.0e-5,
    1.1e-5,
    1.2e-5,
    1.3e-5,
    1.3e-5,
    1.3e-5,
    1.3e-5,
    1.2e-5,
    1.2e-5,
    1.2e-5,
    1.2e-5,
    1.1e-5,
    1.5e-5,
    2.0e-5,
    3.0e-5
  ),
  stringsAsFactors = FALSE
)


library(dplyr)
library(ggplot2)
library(stringr)

plot.data %>%
  # 按 GeneNumber 排序，取前10
  arrange(desc(GeneNumber)) %>%
  slice_head(n = 10) %>%
  mutate(
    Description_wrap = str_wrap(Description, width = 45),
    Description_wrap = factor(Description_wrap, levels = rev(Description_wrap))
  ) %>%
  ggplot(aes(x = RichFactor, y = Description_wrap, fill = -log10(p.adjust))) +
  geom_col(alpha = 0.7, width = 0.8) +

  # 通路名称写在条形左端
  geom_text(
    aes(x = 0.001, label = Description_wrap),
    hjust = 0,
    fontface = "italic",
    size = 5
  ) +

  scale_x_continuous(expand = c(0, 0)) +

  # 渐变色表示显著性
  scale_fill_gradient(
    low = "#4575b4", # 蓝色
    high = "#d73027", # 红色
    name = expression(-log[10](p.adjust))
  ) +
  labs(
    x = "RichFactor",
    y = NULL,
    title = "Top 10 of KEGG Enrichment by GeneNumber"
  ) +
  theme_classic() +
  theme(
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave("富集重画//RichFactor_KEGG_enrichment.pdf", plot = last_plot(), height = 7, width = 7)
