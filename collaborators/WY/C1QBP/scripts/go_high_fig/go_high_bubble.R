# ==============================================================================
# GO Enrichment Bubble Plot — C1QBP+ TAMs（复现参考图）
#
# 数据来源：从参考图中提取的 10 个代表性 GO-BP 通路
# 样式：x=GeneRatio，y=通路名，大小=GeneNumber，颜色=p-value（红→蓝）
# 保存逻辑：与项目其他脚本一致（PDF + PNG + EMF + PPT 矢量）
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
project_dir <- dirname(dirname(script_dir)) # scripts/go_high_fig → scripts → project

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
# 数据（从参考图提取）
# y 轴顺序：图中从上到下排列，ggplot 中从下到上设 factor levels
# --------------------------------------------------------------------------

terms_order <- c(
  "Arp2/3 complex-mediated actin nucleation",
  "ribosome assembly",
  "regulation of actin filament depolymerization",
  "antigen processing and presentation\nof exogenous antigen",
  "phagocytosis",
  "antigen processing and presentation",
  "actin filament polymerization",
  "regulation of actin filament length",
  "regulation of actin polymerization\nor depolymerization",
  "regulation of actin filament polymerization"
)

df <- data.frame(
  Term = factor(terms_order, levels = terms_order),
  GeneRatio = c(
    0.026,
    0.030,
    0.030,
    0.041,
    0.060,
    0.060,
    0.077,
    0.077,
    0.060,
    0.077
  ),
  # 顺序同 terms_order（从下到上）：
  # Arp2/3 | ribosome | reg actin depoly | antigen exog | phago |
  # antigen | actin poly | reg actin length | reg actin poly/depoly | reg actin filament poly
  GeneNumber = c(5, 10, 10, 8, 7, 10, 7, 4, 4, 7),
  pvalue = c(
    0.002,
    0.0015,
    0.0014,
    0.0012,
    0.001,
    0.0008,
    0.0005,
    0.0004,
    0.0003,
    0.0002
  )
)

# --------------------------------------------------------------------------
# 绘图
# --------------------------------------------------------------------------

p_bubble <- ggplot(df, aes(x = GeneRatio, y = Term)) +
  geom_point(
    aes(size = GeneNumber, color = pvalue),
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
    values = rescale(c(0.0001, 0.0005, 0.001, 0.0015, 0.002)),
    limits = c(0.0001, 0.002),
    oob = squish,
    breaks = c(0.001, 0.002),
    labels = c("0.001", "0.002")
  ) +
  scale_x_continuous(
    limits = c(0, 0.10),
    breaks = seq(0, 0.10, by = 0.02),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    title = "GO Enrichment Analysis in C1QBP⁺ TAMs",
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

save_plot(p_bubble, "go_high_bubble", width = 7, height = 6)
cat("[PASS] → go_high_bubble (PDF/PNG/EMF/PPT)\n")

# --------------------------------------------------------------------------
pptx_path <- file.path(fig_dir, "go_high_bubble.pptx")
print(pptx, target = pptx_path)
cat("图形已保存至：", fig_dir, "\n")
cat("PPT 矢量文件：", pptx_path, "\n")
