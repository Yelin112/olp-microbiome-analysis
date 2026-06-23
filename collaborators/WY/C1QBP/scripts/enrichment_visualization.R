# ==============================================================================
# C1QBP 富集分析可视化脚本
#
# 数据：CA_up 差异基因在 C1QBP high / low 组的 GO + KEGG 富集结果
# 图形（x 轴统一使用 GeneRatio）：
#   A. combined  — GO (BP/CC/MF) + KEGG 四分类综合图（high 和 low 各一张）
#   B. lollipop  — High vs Low KEGG 双向对比
#   C. bubble    — 单组 KEGG 气泡图（high 和 low 各一张）
#   D. bar       — 单组 KEGG 条形图（high 和 low 各一张）
#   E. bar dual  — High vs Low KEGG 双向条形图
#   F. point_bar — GO+KEGG 四分类分面折线点+条形
#   G. scatter   — High vs Low KEGG 散点图
#   H. radial    — GO+KEGG 环形树状图
#
# 运行方式：在 RStudio 中打开此文件，点击 Source
# ==============================================================================

library(RColorBrewer)
if (!requireNamespace("devEMF", quietly = TRUE)) {
  install.packages("devEMF")
}
if (!requireNamespace("officer", quietly = TRUE)) {
  install.packages("officer")
}
if (!requireNamespace("rvg", quietly = TRUE)) {
  install.packages("rvg")
}
library(devEMF)
library(officer)
library(rvg)

# --------------------------------------------------------------------------
# 路径设置
# --------------------------------------------------------------------------

script_dir <- dirname(normalizePath(rstudioapi::getSourceEditorContext()$path))
project_dir <- dirname(script_dir)
enrich_plot_dir <- file.path(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/enrich_plot"
)
fig_dir <- file.path(project_dir, "figures", "enrichment")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  base <- file.path(fig_dir, name)

  # PDF + PNG（静态导出）
  ggsave(paste0(base, ".pdf"), p, width = width, height = height)
  ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = 300)

  # EMF（兼容旧流程，部分环境下图形可能丢失）
  devEMF::emf(paste0(base, ".emf"), width = width, height = height)
  print(p)
  dev.off()

  # PowerPoint DrawingML 矢量图（图形+文字均可在 PPT 中直接编辑）
  # 按比例缩放至幻灯片（10 × 7.5 英寸），保持原始宽高比，居中放置
  slide_w <- 10
  slide_h <- 7.5
  scale <- min(slide_w / width, slide_h / height)
  w <- width * scale
  h <- height * scale
  pptx <<- officer::add_slide(pptx, layout = "Blank", master = "Office Theme")
  pptx <<- officer::ph_with(
    pptx,
    rvg::dml(ggobj = p, width = w, height = h),
    location = officer::ph_location(
      left = (slide_w - w) / 2,
      top = (slide_h - h) / 2,
      width = w,
      height = h
    )
  )
}

# PowerPoint 文件初始化（脚本开始时创建，每个 save_plot() 追加一张幻灯片）
pptx <- officer::read_pptx()

# --------------------------------------------------------------------------
# 加载可视化系统
# 注意：不使用 load_enrich_plot.R，因为该文件的路径检测在嵌套 source() 时不可靠
# --------------------------------------------------------------------------

source(file.path(enrich_plot_dir, "R/std_format.R"))
source(file.path(enrich_plot_dir, "R/registry.R"))
for (.f in list.files(
  file.path(enrich_plot_dir, "R/styles"),
  pattern = "\\.R$",
  full.names = TRUE
)) {
  source(.f)
}
rm(.f)
cat("enrich_plot 系统已加载，可用类型:\n")
list_enrich_plots()

# --------------------------------------------------------------------------
# 读取数据
# --------------------------------------------------------------------------

data_dir <- file.path(project_dir, "data", "enrichment")

go_high <- read.csv(
  file.path(data_dir, "high", "CA_up_GO_enrichment_sig.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)
kegg_high <- read.csv(
  file.path(data_dir, "high", "CA_up_KEGG_enrichment_sig.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)
go_low <- read.csv(
  file.path(data_dir, "low", "CA_up_GO_enrichment_sig.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)
kegg_low <- read.csv(
  file.path(data_dir, "low", "CA_up_KEGG_enrichment_sig.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)

# GO 数据按 ONTOLOGY 拆分为 BP / CC / MF 三个分类
split_go <- function(go_df) {
  lapply(c(BP = "BP", CC = "CC", MF = "MF"), function(ont) {
    go_df[go_df$ONTOLOGY == ont, ]
  })
}

go_high_list <- split_go(go_high)
go_low_list <- split_go(go_low)

# --------------------------------------------------------------------------
# 配色方案
# --------------------------------------------------------------------------

category_colors <- c(
  BP = "#437f64",
  CC = "#d88c51",
  MF = "#466277",
  KEGG = "#be5960"
)

group_colors <- c(High = "#3369e7", Low = "#ff6c5f")

# --------------------------------------------------------------------------
# 场景 A：combined 四分类综合图（GO BP/CC/MF + KEGG，x = GeneRatio）
# --------------------------------------------------------------------------

cat("\n=== 场景 A: combined 综合图 ===\n")

## A1: C1QBP High 组
p_a1 <- plot_enrich(
  c(go_high_list, list(KEGG = kegg_high)),
  type = "combined",
  top_n = 10,
  x_var = "gene_ratio",
  category_colors = category_colors,
  title = "CA Up-regulated Enrichment — C1QBP High",
  base_size = 11
)
save_plot(p_a1, "A1_combined_high", 12, 15)
cat("[PASS] A1 → A1_combined_high\n")

## A2: C1QBP Low 组
p_a2 <- plot_enrich(
  c(go_low_list, list(KEGG = kegg_low)),
  type = "combined",
  top_n = 10,
  x_var = "gene_ratio",
  category_colors = category_colors,
  title = "CA Up-regulated Enrichment — C1QBP Low",
  base_size = 11
)
save_plot(p_a2, "A2_combined_low", 12, 15)
cat("[PASS] A2 → A2_combined_low\n")

# --------------------------------------------------------------------------
# 场景 B：lollipop 双向对比（High vs Low，KEGG，x = GeneRatio）
# --------------------------------------------------------------------------

cat("\n=== 场景 B: lollipop High vs Low ===\n")

p_b <- plot_enrich(
  list(High = kegg_high, Low = kegg_low),
  type = "lollipop",
  top_n = 10,
  x_var = "gene_ratio",
  group_colors = group_colors,
  size_var = "gene_count",
  title = "KEGG Enrichment: C1QBP High vs Low",
  base_size = 11
)
save_plot(p_b, "B_lollipop_high_vs_low", 10, 7)
cat("[PASS] B → B_lollipop_high_vs_low\n")

# --------------------------------------------------------------------------
# 场景 C：bubble 单组气泡图（KEGG，x = GeneRatio）
# --------------------------------------------------------------------------

cat("\n=== 场景 C: bubble 单组气泡图 ===\n")

## C1: C1QBP High 组（两色渐变）
p_c1 <- plot_enrich(
  kegg_high,
  type = "bubble",
  top_n = 10,
  x_var = "gene_ratio",
  palette = rev(brewer.pal(5, "Spectral")),
  title = "KEGG Enrichment — C1QBP High",
  base_size = 11
)
save_plot(p_c1, "C1_bubble_kegg_high", 9, 6)
cat("[PASS] C1 → C1_bubble_kegg_high\n")

## C2: C1QBP Low 组（两色渐变）
p_c2 <- plot_enrich(
  kegg_low,
  type = "bubble",
  top_n = 10,
  x_var = "gene_ratio",
  palette = rev(brewer.pal(5, "Spectral")),
  title = "KEGG Enrichment — C1QBP Low",
  base_size = 11
)
save_plot(p_c2, "C2_bubble_kegg_low", 9, 6)
cat("[PASS] C2 → C2_bubble_kegg_low\n")

# --------------------------------------------------------------------------
# 场景 D：bar 单组（KEGG，x = GeneRatio）
# --------------------------------------------------------------------------

cat("\n=== 场景 D: bar 单组 ===\n")

## D1: C1QBP High 组（默认两色渐变）
p_d1 <- plot_enrich(
  kegg_high,
  type = "bar",
  top_n = 10,
  x_var = "gene_ratio",
  group_colors = rev(brewer.pal(5, "Spectral")),
  title = "KEGG Enrichment — C1QBP High",
  base_size = 11
)
save_plot(p_d1, "D1_bar_kegg_high", 9, 6)
cat("[PASS] D1 → D1_bar_kegg_high\n")

## D2: C1QBP Low 组（Spectral 多色板，x = GeneRatio，颜色 = p.adjust）
p_d2 <- plot_enrich(
  kegg_low,
  type = "bar",
  top_n = 10,
  x_var = "gene_ratio",
  color_var = "p_adjust",
  group_colors = rev(brewer.pal(5, "Spectral")),
  title = "KEGG Enrichment — C1QBP Low",
  base_size = 11
)
save_plot(p_d2, "D2_bar_kegg_low", 9, 6)
cat("[PASS] D2 → D2_bar_kegg_low\n")

# --------------------------------------------------------------------------
# 场景 E：bar 双向对比（High vs Low，KEGG，x = GeneRatio）
# --------------------------------------------------------------------------

cat("\n=== 场景 E: bar 双向对比 ===\n")

p_e <- plot_enrich(
  list(High = kegg_high, Low = kegg_low),
  type = "bar",
  top_n = 10,
  x_var = "gene_ratio",
  group_colors = group_colors,
  title = "KEGG Enrichment: C1QBP High vs Low",
  base_size = 11
)
save_plot(p_e, "E_bar_high_vs_low", 10, 7)
cat("[PASS] E → E_bar_high_vs_low\n")

# --------------------------------------------------------------------------
# 场景 F：point_bar 分面折线点+条形（GO+KEGG，x 左=GeneRatio，右=-log10p）
# --------------------------------------------------------------------------

cat("\n=== 场景 F: point_bar 分面图 ===\n")

p_f <- plot_enrich(
  c(go_high_list, list(KEGG = kegg_high)),
  type = "point_bar",
  top_n = 10,
  category_colors = category_colors,
  title = "GO & KEGG Enrichment — C1QBP High",
  base_size = 11
)
save_plot(p_f, "F_point_bar_high", 14, 14)
cat("[PASS] F → F_point_bar_high\n")

# --------------------------------------------------------------------------
# 场景 G：scatter 双向散点图（High vs Low，KEGG）
# --------------------------------------------------------------------------

cat("\n=== 场景 G: scatter 双向散点 ===\n")

p_g <- plot_enrich(
  list(High = kegg_high, Low = kegg_low),
  type = "scatter",
  top_n = 10,
  x_var = "gene_count",
  y_var = "p_adjust",
  size_var = "gene_ratio",
  point_shape = 21,
  stroke_color = "black",
  bg_color = "#f7f7f7",
  group_colors = group_colors,
  title = "KEGG Enrichment: C1QBP High vs Low",
  base_size = 11
)
save_plot(p_g, "G_scatter_high_vs_low", 10, 7)
cat("[PASS] G → G_scatter_high_vs_low\n")

# --------------------------------------------------------------------------
# 场景 H：radial 环形树状图（GO+KEGG）
# --------------------------------------------------------------------------

cat("\n=== 场景 H: radial 环形图 ===\n")

p_h <- plot_enrich(
  c(go_high_list, list(KEGG = kegg_high)),
  type = "radial",
  top_n = 10,
  inner_radius = -150,
  inset_left = 0.12,
  inset_right = 0.88,
  inset_bottom = 0.12,
  inset_top = 0.88,
  category_colors = category_colors,
  title = "GO & KEGG Enrichment — C1QBP High",
  base_size = 15
)
save_plot(p_h, "H_radial_high", 10, 10)
cat("[PASS] H → H_radial_high\n")


p_l <- plot_enrich(
  c(go_high_list, list(KEGG = kegg_low)),
  type = "radial",
  top_n = 10,
  inner_radius = -150,
  inset_left = 0.12,
  inset_right = 0.88,
  inset_bottom = 0.12,
  inset_top = 0.88,
  category_colors = category_colors,
  title = "GO & KEGG Enrichment — C1QBP Low",
  base_size = 15
)
save_plot(p_l, "H_radial_low", 10, 10)
cat("[PASS] H → H_radial_low\n")

# --------------------------------------------------------------------------
pptx_path <- file.path(fig_dir, "all_enrichment_plots.pptx")
print(pptx, target = pptx_path)
cat("\n所有图形已保存至：", fig_dir, "\n")
cat("PPT 矢量文件：", pptx_path, "\n")
