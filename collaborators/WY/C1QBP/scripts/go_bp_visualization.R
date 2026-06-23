# ==============================================================================
# C1QBP GO-BP 富集分析可视化脚本
#
# 重心：GO Biological Process（BP）
# 排序依据：校正后 p 值（p.adjust）升序
# 数据：CA_up 差异基因在 C1QBP high / low 组的 GO 富集结果
#
# 图形：
#   A. BP bubble   — High / Low 各一张气泡图（BP，按 p.adjust 排序）
#   B. BP bar      — High / Low 各一张条形图（BP，颜色映射 p.adjust）
#   C. BP lollipop — High vs Low 双向对比棒棒糖图（BP）
#   D. BP scatter  — High vs Low 双向散点图（BP）
#   E. GO combined — BP/CC/MF 三分类综合图（High / Low 各一张）
#   F. BP dual bar — High vs Low 双向条形图（BP）
#
# 运行方式：在 RStudio 中打开此文件，点击 Source
# ==============================================================================

library(ggplot2)
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

enrich_plot_dir <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/enrich_plot"

fig_dir <- file.path(project_dir, "figures", "go_bp")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 保存函数：PDF + PNG + EMF + PPT 矢量
# --------------------------------------------------------------------------

pptx <- officer::read_pptx()

save_plot <- function(p, name, width, height) {
  base <- file.path(fig_dir, name)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height)
  ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = 300)

  devEMF::emf(paste0(base, ".emf"), width = width, height = height)
  print(p)
  dev.off()

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

# --------------------------------------------------------------------------
# 加载可视化系统
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
go_low <- read.csv(
  file.path(data_dir, "low", "CA_up_GO_enrichment_sig.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)

# 按 p.adjust 升序确保后续 top_n 选取最显著的 term
go_high <- go_high[order(go_high$p.adjust), ]
go_low <- go_low[order(go_low$p.adjust), ]

# 拆分 BP / CC / MF
bp_high <- go_high[go_high$ONTOLOGY == "BP", ]
cc_high <- go_high[go_high$ONTOLOGY == "CC", ]
mf_high <- go_high[go_high$ONTOLOGY == "MF", ]

bp_low <- go_low[go_low$ONTOLOGY == "BP", ]
cc_low <- go_low[go_low$ONTOLOGY == "CC", ]
mf_low <- go_low[go_low$ONTOLOGY == "MF", ]

cat(sprintf(
  "\n数据概况 — High: BP=%d / CC=%d / MF=%d | Low: BP=%d / CC=%d / MF=%d\n",
  nrow(bp_high),
  nrow(cc_high),
  nrow(mf_high),
  nrow(bp_low),
  nrow(cc_low),
  nrow(mf_low)
))

# --------------------------------------------------------------------------
# 配色
# --------------------------------------------------------------------------

ont_colors <- c(BP = "#437f64", CC = "#d88c51", MF = "#466277")
group_colors <- c(High = "#3369e7", Low = "#ff6c5f")
bp_palette <- rev(brewer.pal(9, "RdYlBu")) # 蓝→黄→红渐变，低 p 值为深红

# --------------------------------------------------------------------------
# 场景 A：BP 气泡图（按 p.adjust 排序，x = GeneRatio）
# --------------------------------------------------------------------------

cat("\n=== 场景 A: BP bubble ===\n")

p_a1 <- plot_enrich(
  bp_high,
  type = "bubble",
  top_n = 10,
  x_var = "gene_ratio",
  palette = bp_palette,
  title = "GO-BP Enrichment — C1QBP High",
  base_size = 11
)
save_plot(p_a1, "A1_bp_bubble_high", 9, 6)
cat("[PASS] A1 → A1_bp_bubble_high\n")

p_a2 <- plot_enrich(
  bp_low,
  type = "bubble",
  top_n = 10,
  x_var = "gene_ratio",
  palette = bp_palette,
  title = "GO-BP Enrichment — C1QBP Low",
  base_size = 11
)
save_plot(p_a2, "A2_bp_bubble_low", 9, 6)
cat("[PASS] A2 → A2_bp_bubble_low\n")

# --------------------------------------------------------------------------
# 场景 B：BP 条形图（颜色 = p.adjust，按 p.adjust 排序）
# --------------------------------------------------------------------------

cat("\n=== 场景 B: BP bar ===\n")

p_b1 <- plot_enrich(
  bp_high,
  type = "bar",
  top_n = 10,
  x_var = "gene_ratio",
  color_var = "p_adjust",
  group_colors = bp_palette,
  title = "GO-BP Enrichment — C1QBP High",
  base_size = 11
)
save_plot(p_b1, "B1_bp_bar_high", 9, 6)
cat("[PASS] B1 → B1_bp_bar_high\n")

p_b2 <- plot_enrich(
  bp_low,
  type = "bar",
  top_n = 10,
  x_var = "gene_ratio",
  color_var = "p_adjust",
  group_colors = bp_palette,
  title = "GO-BP Enrichment — C1QBP Low",
  base_size = 11
)
save_plot(p_b2, "B2_bp_bar_low", 9, 6)
cat("[PASS] B2 → B2_bp_bar_low\n")

# --------------------------------------------------------------------------
# 场景 C：BP 双向棒棒糖（High vs Low，按 p.adjust 排序）
# --------------------------------------------------------------------------

cat("\n=== 场景 C: BP lollipop High vs Low ===\n")

p_c <- plot_enrich(
  list(High = bp_high, Low = bp_low),
  type = "lollipop",
  top_n = 10,
  x_var = "gene_ratio",
  group_colors = group_colors,
  size_var = "gene_count",
  title = "GO-BP Enrichment: C1QBP High vs Low",
  base_size = 11
)
save_plot(p_c, "C_bp_lollipop_high_vs_low", 10, 7)
cat("[PASS] C → C_bp_lollipop_high_vs_low\n")

# --------------------------------------------------------------------------
# 场景 D：BP 双向散点图（High vs Low）
# --------------------------------------------------------------------------

cat("\n=== 场景 D: BP scatter High vs Low ===\n")

p_d <- plot_enrich(
  list(High = bp_high, Low = bp_low),
  type = "scatter",
  top_n = 10,
  x_var = "gene_count",
  y_var = "p_adjust",
  size_var = "gene_ratio",
  point_shape = 21,
  stroke_color = "black",
  bg_color = "#f7f7f7",
  group_colors = group_colors,
  title = "GO-BP Enrichment: C1QBP High vs Low",
  base_size = 11
)
save_plot(p_d, "D_bp_scatter_high_vs_low", 10, 7)
cat("[PASS] D → D_bp_scatter_high_vs_low\n")

# --------------------------------------------------------------------------
# 场景 E：GO 三分类综合图（BP/CC/MF，按 p.adjust 排序，High / Low 各一张）
# --------------------------------------------------------------------------

cat("\n=== 场景 E: GO combined (BP/CC/MF) ===\n")

p_e1 <- plot_enrich(
  list(BP = bp_high, CC = cc_high, MF = mf_high),
  type = "combined",
  top_n = 10,
  x_var = "gene_ratio",
  category_colors = ont_colors,
  title = "GO Enrichment (BP/CC/MF) — C1QBP High",
  base_size = 11
)
save_plot(p_e1, "E1_go_combined_high", 12, 14)
cat("[PASS] E1 → E1_go_combined_high\n")

p_e2 <- plot_enrich(
  list(BP = bp_low, CC = cc_low, MF = mf_low),
  type = "combined",
  top_n = 10,
  x_var = "gene_ratio",
  category_colors = ont_colors,
  title = "GO Enrichment (BP/CC/MF) — C1QBP Low",
  base_size = 11
)
save_plot(p_e2, "E2_go_combined_low", 12, 14)
cat("[PASS] E2 → E2_go_combined_low\n")

# --------------------------------------------------------------------------
# 场景 F：BP 双向条形图（High vs Low）
# --------------------------------------------------------------------------

cat("\n=== 场景 F: BP bar High vs Low ===\n")

p_f <- plot_enrich(
  list(High = bp_high, Low = bp_low),
  type = "bar",
  top_n = 10,
  x_var = "gene_ratio",
  group_colors = group_colors,
  title = "GO-BP Enrichment: C1QBP High vs Low",
  base_size = 11
)
save_plot(p_f, "F_bp_bar_high_vs_low", 10, 7)
cat("[PASS] F → F_bp_bar_high_vs_low\n")

# --------------------------------------------------------------------------
# 保存 PPT
# --------------------------------------------------------------------------

pptx_path <- file.path(fig_dir, "go_bp_enrichment_plots.pptx")
print(pptx, target = pptx_path)
cat("\n所有图形已保存至：", fig_dir, "\n")
cat("PPT 矢量文件：", pptx_path, "\n")
