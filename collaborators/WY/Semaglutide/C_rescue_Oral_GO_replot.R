# ================================================================
# Oral C_rescue 聚类 GO BP 富集结果 —— 指定通路重新可视化
# 输入：analysis_results/RNA/Rescue/Oral/C_rescue_cluster_GO_BP.csv
# 输出：analysis_results/RNA/Rescue/Oral/
# ================================================================

library(ggplot2)
library(dplyr)
library(stringr)

# ────────────────────────────────────────────────────────────────
# 1. 路径
# ────────────────────────────────────────────────────────────────

INPUT_CSV <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/analysis_results/RNA/Rescue/Oral/C_rescue_cluster_GO_BP.csv"
OUT_DIR <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/analysis_results/RNA/Rescue/Oral"

# ────────────────────────────────────────────────────────────────
# 2. 选定通路
# ────────────────────────────────────────────────────────────────

KEEP_TERMS <- c(
  # 核心（必留）
  "defense response to virus",
  "response to virus",
  "response to interferon-beta",
  "cellular response to interferon-beta",
  "interferon-mediated signaling pathway",
  # 辅助
  "regulation of cytoplasmic pattern recognition receptor signaling pathway",
  "positive regulation of interferon-beta production",
  "inflammasome-mediated signaling pathway"
)

# ────────────────────────────────────────────────────────────────
# 3. 读取 & 筛选
# ────────────────────────────────────────────────────────────────

dat_raw <- read.csv(INPUT_CSV, header = TRUE, stringsAsFactors = FALSE)

dat <- dat_raw %>%
  filter(Description %in% KEEP_TERMS) %>%
  mutate(
    # GeneRatio 字符串转数值
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- as.numeric(strsplit(x, "/")[[1]])
      parts[1] / parts[2]
    }),
    neg_log10_padj = -log10(p.adjust),
    # 是否显著（p.adjust < 0.05）
    sig_label = ifelse(p.adjust < 0.05, "p.adj < 0.05", "p.adj ≥ 0.05")
  )

# 按 FoldEnrichment 排序，用于 Y 轴顺序
dat <- dat %>%
  arrange(FoldEnrichment) %>%
  mutate(Description = factor(Description, levels = Description))

cat(sprintf("筛选到 %d 条通路\n", nrow(dat)))
print(dat[, c(
  "Description",
  "Count",
  "GeneRatio",
  "FoldEnrichment",
  "p.adjust"
)])

# ────────────────────────────────────────────────────────────────
# 4. 自动换行（长标签）
# ────────────────────────────────────────────────────────────────

wrap_labels <- function(x, width = 40) {
  str_wrap(as.character(x), width = width)
}

dat$Description_wrap <- factor(
  wrap_labels(dat$Description),
  levels = wrap_labels(levels(dat$Description))
)

# ────────────────────────────────────────────────────────────────
# 5. 气泡图（Dot plot）
# ────────────────────────────────────────────────────────────────

p_dot <- ggplot(
  dat,
  aes(
    x = GeneRatio_num,
    y = Description_wrap,
    size = Count,
    fill = neg_log10_padj
  )
) +
  geom_point(
    shape = 21,
    color = "grey30",
    stroke = 0.4
  ) +
  scale_size_continuous(
    name = "Gene count",
    range = c(4, 12),
    breaks = pretty(dat$Count, n = 4)
  ) +
  scale_fill_gradient(
    name = expression(-log[10](p.adj)),
    low = "#74C5E3",
    high = "#B2182B",
    na.value = "grey80"
  ) +
  scale_x_continuous(
    name = "Gene ratio",
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = "GO Biological Process",
    subtitle = "Oral C-rescue cluster — selected pathways",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    legend.position = "right",
    legend.key.size = unit(0.9, "lines"),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor = element_blank()
  )

# ────────────────────────────────────────────────────────────────
# 6. 条形图（Bar plot，按 FoldEnrichment，-log10 p.adj 填色）
# ────────────────────────────────────────────────────────────────

p_bar <- ggplot(
  dat,
  aes(
    x = FoldEnrichment,
    y = Description_wrap,
    fill = neg_log10_padj
  )
) +
  geom_col(color = "grey30", linewidth = 0.3, width = 0.7) +
  geom_text(
    aes(label = paste0("n=", Count)),
    hjust = -0.15,
    size = 3.5,
    color = "grey30"
  ) +
  scale_fill_gradient(
    name = expression(-log[10](p.adj)),
    low = "#74C5E3",
    high = "#B2182B",
    na.value = "grey80"
  ) +
  scale_x_continuous(
    name = "Fold enrichment",
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "GO Biological Process",
    subtitle = "Oral C-rescue cluster — selected pathways",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    legend.position = "right",
    legend.key.size = unit(0.9, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ────────────────────────────────────────────────────────────────
# 7. 保存
# ────────────────────────────────────────────────────────────────

save_plot <- function(p, base_name, w = 8, h = 5) {
  pdf_path <- file.path(OUT_DIR, paste0(base_name, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(base_name, ".png"))
  ggsave(pdf_path, plot = p, width = w, height = h)
  ggsave(png_path, plot = p, width = w, height = h, dpi = 300)
  message("已保存: ", pdf_path)
  message("已保存: ", png_path)
}

save_plot(p_dot, "C_rescue_Oral_GO_dotplot", w = 8.5, h = 5.5)
save_plot(p_bar, "C_rescue_Oral_GO_barplot", w = 8.5, h = 5.5)

# ────────────────────────────────────────────────────────────────
# 8. 仿公众号风格：圆圈 + 横条 + 通路名 + 基因名
# ────────────────────────────────────────────────────────────────

# 按 FoldEnrichment 升序排列（y 轴：最高富集在顶部）
dat_s <- dat %>% arrange(FoldEnrichment)

n_s <- nrow(dat_s)
Y_STEP <- 3.2
dat_s$y <- seq(Y_STEP, n_s * Y_STEP, by = Y_STEP)

# 基因名截断显示（最多展示 10 个）
MAX_GENES <- 10
dat_s$gene_display <- sapply(dat_s$geneID, function(x) {
  g <- strsplit(x, "/")[[1]]
  if (length(g) > MAX_GENES) {
    paste(c(g[1:MAX_GENES], "..."), collapse = "/")
  } else {
    x
  }
})

# 通路名超过 42 字自动换行
dat_s$desc_wrap <- str_wrap(dat_s$Description, width = 42)

# 条形长度：FoldEnrichment 等比缩放到 BAR_MAX 宽度
BAR_MAX <- 32
bar_scale <- BAR_MAX / max(dat_s$FoldEnrichment)
dat_s$bar_end <- dat_s$FoldEnrichment * bar_scale

# ── X 布局坐标 ───────────────────────────────────────────────
X_BAR_S <- 0 # 条形起点
X_TEXT <- 0.5 # 文字起点
X_CIR <- -5.0 # 圆圈中心
X_CAT_L <- -9.0 # 分类框左边界
X_CAT_R <- -6.5 # 分类框右边界
X_CAT_M <- (X_CAT_L + X_CAT_R) / 2

# ── 颜色 ─────────────────────────────────────────────────────
CAT_COL <- "#CC79A7" # GO BP 紫粉色（与参考图 BP 配色一致）
GENE_COL <- "#7B1D5E" # 基因名文字（同系深色）

# ── Y 范围 ───────────────────────────────────────────────────
y_lo <- min(dat_s$y) - Y_STEP * 0.55
y_hi <- max(dat_s$y) + Y_STEP * 0.55

# ── X 轴刻度（还原真实 FoldEnrichment 值）───────────────────
FE_seq <- pretty(c(0, max(dat_s$FoldEnrichment)), n = 5)
FE_seq <- FE_seq[FE_seq >= 0]
x_breaks <- FE_seq * bar_scale

p_style2 <- ggplot(dat_s) +

  # 分类标签框
  annotate(
    "rect",
    xmin = X_CAT_L,
    xmax = X_CAT_R,
    ymin = y_lo,
    ymax = y_hi,
    fill = CAT_COL,
    color = NA
  ) +
  annotate(
    "text",
    x = X_CAT_M,
    y = mean(dat_s$y),
    label = "BP",
    color = "white",
    fontface = "bold",
    size = 5.5,
    angle = 90
  ) +

  # 横条（矩形）
  geom_rect(
    aes(
      xmin = X_BAR_S,
      xmax = bar_end,
      ymin = y - 0.75,
      ymax = y + 0.75
    ),
    fill = CAT_COL,
    alpha = 0.78,
    color = NA
  ) +

  # 通路名（黑色粗体，条形内左对齐）
  geom_text(
    aes(
      x = X_TEXT,
      y = y + 0.18,
      label = desc_wrap
    ),
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fontface = "bold",
    color = "black",
    lineheight = 0.9
  ) +

  # 基因名（彩色斜体小字，通路名下方）
  geom_text(
    aes(
      x = X_TEXT,
      y = y - 0.42,
      label = gene_display
    ),
    hjust = 0,
    size = 2.5,
    color = GENE_COL,
    fontface = "italic"
  ) +

  # 圆圈（大小正比于 Count）
  geom_point(
    aes(
      x = X_CIR,
      y = y,
      size = Count
    ),
    shape = 21,
    fill = CAT_COL,
    color = "white",
    stroke = 0.9
  ) +

  # 圆圈内数字
  geom_text(
    aes(
      x = X_CIR,
      y = y,
      label = Count
    ),
    size = 3.0,
    color = "white",
    fontface = "bold"
  ) +

  scale_size_continuous(
    name = "Gene count",
    range = c(5, 14),
    breaks = c(1, 3, 5, 8)
  ) +
  scale_x_continuous(
    name = "Fold enrichment",
    breaks = x_breaks,
    labels = FE_seq,
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(
    xlim = c(X_CAT_L - 0.3, max(dat_s$bar_end) + 1.5),
    ylim = c(y_lo, y_hi),
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_text(
      size = 10,
      color = "grey30",
      margin = margin(t = 3)
    ),
    axis.title.x = element_text(
      size = 11,
      color = "grey30",
      margin = margin(t = 7)
    ),
    axis.ticks.x = element_line(color = "grey60", linewidth = 0.4),
    axis.line.x = element_line(color = "grey60", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(15, 70, 30, 10, "pt"),
    plot.background = element_rect(fill = "white", color = NA)
  )

save_plot(p_style2, "C_rescue_Oral_GO_custom_style", w = 11, h = 7.5)

message("\n完成！输出目录：", OUT_DIR)
