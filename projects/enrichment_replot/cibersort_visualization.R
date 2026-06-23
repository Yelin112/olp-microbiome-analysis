library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

outdir <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/富集重画"

# ── 读取数据 ──────────────────────────────────────────────────────────────────
df <- read.csv(
  file.path(outdir, "Cibersort_data_long.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)
df$Composition <- as.numeric(df$Composition)

# ── 1. 热图：样本 × 细胞类型 ─────────────────────────────────────────────────

# 宽矩阵: 行 = 细胞类型, 列 = 样本
mat <- df %>%
  select(sample, Celltype, Composition) %>%
  pivot_wider(names_from = sample, values_from = Composition) %>%
  tibble::column_to_rownames("Celltype") %>%
  as.matrix()

# 样本注释（group 信息）
sample_anno <- df %>%
  distinct(sample, group) %>%
  arrange(group, sample) # OLP/HC 分组排列

mat <- mat[, sample_anno$sample] # 按组排列列

group_col <- c(OLP = "#E64B35", HC = "#4DBBD5")

col_anno <- HeatmapAnnotation(
  Group = sample_anno$group,
  col = list(Group = group_col),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(4, "mm")
)

# 颜色：YlGnBu 配色
col_fun <- colorRamp2(
  seq(0, max(mat, na.rm = TRUE), length.out = 9),
  brewer.pal(9, "YlGnBu")
)

# 细胞类型标签美化（去点替换为空格）
rownames(mat) <- gsub("\\.", " ", rownames(mat))

# 正方形单元格：行列数决定矩阵尺寸
cell_size <- unit(5, "mm") # 每格 5mm × 5mm

ht <- Heatmap(
  mat,
  name = "Composition",
  col = col_fun,
  width = ncol(mat) * cell_size,
  height = nrow(mat) * cell_size,
  top_annotation = col_anno,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_side = "right",
  row_dend_width = unit(15, "mm"),
  border = TRUE,
  rect_gp = gpar(col = "grey20", lwd = 0.3),
  heatmap_legend_param = list(
    title = "Composition",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  ),
  column_split = sample_anno$group,
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_gap = unit(3, "mm")
)

# PDF/PNG 尺寸根据矩阵自动计算（留出 label、legend、annotation 空间）
pdf_w <- ncol(mat) * 5 / 25.4 + 5 # 英寸
pdf_h <- nrow(mat) * 5 / 25.4 + 3

pdf(file.path(outdir, "cibersort_heatmap.pdf"), width = pdf_w, height = pdf_h)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(
  file.path(outdir, "cibersort_heatmap.png"),
  width = pdf_w,
  height = pdf_h,
  units = "in",
  res = 300
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("Heatmap saved.\n")

# ── 2. 组间比较：各细胞类型 Composition 箱线图 ───────────────────────────────
source(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/lib/compare_plot/compare_plot_optimized.R"
)

# Wilcoxon 检验（保存统计结果）
celltypes <- unique(df$Celltype)
stat_res <- lapply(celltypes, function(ct) {
  sub <- df[df$Celltype == ct, ]
  wt <- wilcox.test(Composition ~ group, data = sub, exact = FALSE)
  data.frame(Celltype = ct, p_value = wt$p.value)
}) %>%
  bind_rows()

stat_res$p_adj <- p.adjust(stat_res$p_value, method = "BH")
stat_res$sig <- ifelse(
  stat_res$p_adj < 0.001,
  "***",
  ifelse(stat_res$p_adj < 0.01, "**", ifelse(stat_res$p_adj < 0.05, "*", "ns"))
)

write.csv(
  stat_res,
  file.path(outdir, "cibersort_wilcox_stats.csv"),
  row.names = FALSE
)
cat("Wilcoxon results saved.\n")

# 标签美化，固定细胞类型顺序
df$Celltype_label <- gsub("\\.", " ", df$Celltype)
ct_order <- unique(gsub("\\.", " ", celltypes))
df$Celltype_label <- factor(df$Celltype_label, levels = ct_order)
df$group <- factor(df$group, levels = c("OLP", "HC"))

# 使用 compare_plot 绘制，自动识别样本量选择策略（≥20/组 → pure_boxplot）
p_box <- compare_plot(
  data = df,
  value.var = "Composition",
  group.by = "group",
  split.by = "Celltype_label",
  add_stat = "wilcox.test",
  comparisons = list(c("OLP", "HC")),
  hide_ns = FALSE,
  stat_label = "p.signif",
  palette = c("#E64B35", "#4DBBD5"),
  strategy = "auto",
  add_point = TRUE,
  # point_size = 1.8,
  # point_alpha = 0.55,
  # point_jitter = 0.15,
  xlab = NULL,
  ylab = "Composition",
  title = "Immune cell composition: OLP vs HC",
  theme_use = theme_cowplot
) +
  facet_wrap(~Celltype_label, scales = "free_y", ncol = 5) +
  theme(
    text = element_text(face = "bold", size = 16), # 所有字体加粗
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # X轴标签旋转45度
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5), #添加外框线
    # 有了外框后，隐藏单独的坐标轴线以保持整洁
    axis.line = element_blank()
  )

ggsave(
  file.path(outdir, "cibersort_boxplot.pdf"),
  p_box,
  width = 14,
  height = 14
)
ggsave(
  file.path(outdir, "cibersort_boxplot.png"),
  p_box,
  width = 14,
  height = 14,
  dpi = 300
)

cat("Boxplot saved.\n")
cat("\nDone. Output files:\n")
cat("  cibersort_heatmap.pdf / .png\n")
cat("  cibersort_boxplot.pdf / .png\n")
cat("  cibersort_wilcox_stats.csv\n")
