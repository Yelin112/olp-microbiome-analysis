# ==============================================================================
# KEGG Enrichment Bubble Plot — C1QBP+ TAMs
#
# 数据来源：high 组 KEGG 富集结果（CA_up_KEGG_enrichment_sig.csv）
# 通路：从参考图中选定的 10 条代表性 KEGG 通路
# 样式：x=GeneRatio，y=通路名（按 p.adjust 升序），大小=Count，颜色=p.adjust
# ==============================================================================

library(ggplot2)
library(scales)

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
project_dir <- dirname(dirname(script_dir))

data_dir <- file.path(project_dir, "data", "enrichment")
fig_dir <- file.path(project_dir, "figures", "go_high_fig")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 保存函数
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
  sc <- min(slide_w / width, slide_h / height)
  w <- width * sc
  h <- height * sc
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
# 读取数据，筛选目标通路
# --------------------------------------------------------------------------

kegg_high <- read.csv(
  file.path(data_dir, "high", "CA_up_KEGG_enrichment_sig.csv"),
  row.names = 1,
  stringsAsFactors = FALSE
)

target_terms <- c(
  "Proteasome",
  "Phagosome",
  "Ribosome",
  "Fc gamma R-mediated phagocytosis",
  "Regulation of actin cytoskeleton",
  "Antigen processing and presentation",
  "Apoptosis",
  "Neutrophil extracellular trap formation",
  "RNA degradation",
  "Oxidative phosphorylation"
)

df <- kegg_high[kegg_high$Description %in% target_terms, ]

# 解析 GeneRatio 为数值
ratio_parts <- strsplit(df$GeneRatio, "/")
df$GeneRatio_num <- as.numeric(sapply(ratio_parts, `[`, 1)) /
  as.numeric(sapply(ratio_parts, `[`, 2))

# y 轴排序：按 p.adjust 升序（最显著在最上方）
df <- df[order(df$p.adjust), ]
df$Term <- factor(df$Description, levels = rev(df$Description))

cat(sprintf("找到 %d / %d 条目标通路\n", nrow(df), length(target_terms)))
print(df[, c("Description", "GeneRatio", "p.adjust", "Count")])

# --------------------------------------------------------------------------
# 绘图
# --------------------------------------------------------------------------

p_bubble <- ggplot(df, aes(x = GeneRatio_num, y = Term)) +
  geom_point(
    aes(size = Count, color = p.adjust),
    shape = 16
  ) +
  scale_size_continuous(
    name = "GeneNumber",
    range = c(4, 10),
    breaks = c(4, 7, 10)
  ) +
  scale_color_gradientn(
    name = "p-value",
    colors = c("#D62828", "#F77F00", "#FCBF49", "#AED9E0", "#2B6CB0"),
    values = rescale(c(
      min(df$p.adjust),
      quantile(df$p.adjust, 0.25),
      median(df$p.adjust),
      quantile(df$p.adjust, 0.75),
      max(df$p.adjust)
    )),
    limits = c(min(df$p.adjust), max(df$p.adjust)),
    breaks = signif(quantile(df$p.adjust, c(0.25, 0.75)), 2),
    labels = function(x) formatC(x, format = "e", digits = 1)
  ) +
  scale_x_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    title = "KEGG Enrichment Analysis in C1QBP⁺ TAMs",
    x = "GeneRatio",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.key.height = unit(0.6, "cm")
  )

save_plot(p_bubble, "kegg_high_bubble", width = 6.5, height = 6)
cat("[PASS] → kegg_high_bubble (PDF/PNG/EMF/PPT)\n")

# --------------------------------------------------------------------------
pptx_path <- file.path(fig_dir, "kegg_high_bubble.pptx")
print(pptx, target = pptx_path)
cat("图形已保存至：", fig_dir, "\n")
cat("PPT 矢量文件：", pptx_path, "\n")
