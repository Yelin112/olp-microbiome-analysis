# C1QBP 单细胞数据组间比较分析
# 分析：C1QBP 表达量在 NL/LP/CA/LN 四组间的差异

# ==============================================================================
# 0. 依赖包
# ==============================================================================
library(ggplot2)
library(ggpubr)
library(dplyr)
library(paletteer)
library(officer)
library(rvg)

# 加载自定义函数
source(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/archive/compare_plot_optimized_v1_20250405.R"
)
source("E:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/my_themes.R")

# ==============================================================================
# 1. 数据读取与预处理
# ==============================================================================
meta <- read.csv(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/C1QBP/data/C1QBP_Sc/metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# 设置有序因子
group_levels <- c("NL", "LP", "CA", "LN")
meta$group <- factor(meta$group, levels = group_levels)

# 细胞类型保持原始水平（非等级）
meta$cell.type <- factor(meta$cell.type)
cell_types <- levels(meta$cell.type)

# ==============================================================================
# 2. 配色方案
# ==============================================================================
# 4组用 Shuksan2，备选 Temps
group_colors <- as.character(paletteer::paletteer_d("rcartocolor::Temps"))[2:6]
names(group_colors) <- group_levels

# ==============================================================================
# 3. 比较对设置（相邻组 + 首尾组）
# ==============================================================================
pairwise_comparisons <- list(
  c("NL", "LP"),
  c("LP", "CA"),
  c("CA", "LN"),
  c("NL", "CA"),
  c("NL", "LN"),
  c("LP", "LN")
)

# ==============================================================================
# 4. 辅助：保存函数（PDF + PNG + PPTX）
# ==============================================================================
save_plot_all <- function(p, filename_base, width = 6, height = 5) {
  dir_out <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/C1QBP/figures/C1QBP_exp_compare"

  # PDF
  ggsave(
    file.path(dir_out, paste0(filename_base, ".pdf")),
    plot = p,
    width = width,   
    height = height,
    device = cairo_pdf
  )

  # PNG
  ggsave(
    file.path(dir_out, paste0(filename_base, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )

  # PPTX (officer + rvg)
  pptx_path <- file.path(dir_out, paste0(filename_base, ".pptx"))
  pptx <- read_pptx()
  pptx <- add_slide(pptx, layout = "Blank", master = "Office Theme")
  pptx <- ph_with(
    pptx,
    value = dml(ggobj = p),
    location = ph_location(
      left = 0.5,
      top = 0.5,
      width = width,
      height = height
    )
  )
  print(pptx, target = pptx_path)

  message("已保存：", filename_base)
}

# ==============================================================================
# 5. 图1：全细胞类型汇总 —— C1QBP 各组表达量比较
# ==============================================================================
p_overall <- compare_plot(
  data = meta,
  value.var = "C1QBP",
  group.by = "group",
  plot_type = "violin",
  violin_alpha = 0.6, # 更透明
  violin_linewidth = 1.0, # 更粗的边框
  add_box = TRUE,
  add_stat = "wilcox.test",
  stat_label = "p.signif",
  comparisons = pairwise_comparisons,
  hide_ns = TRUE,
  step_increase = 0.1,
  y_expand = 0.25,
  palette = group_colors,
  xlab = "Group",
  ylab = "C1QBP Expression",
  title = "C1QBP Expression Across Groups (All Cells)",
  theme_use = theme_my_stat,
  base_size = 12
) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

save_plot_all(p_overall, "C1QBP_overall_group_compare", width = 6, height = 5.5)

# ==============================================================================
# 6. 图2：各细胞类型分面图 —— 一张图展示所有细胞类型
# ==============================================================================
p_facet <- compare_plot(
  data = meta,
  value.var = "C1QBP",
  group.by = "group",
  split.by = "cell.type",
  plot_type = "violin",
  violin_alpha = 0.6, # 更透明
  violin_linewidth = 1.0, # 更粗的边框
  add_box = TRUE,
  add_stat = "wilcox.test",
  stat_label = "p.signif",
  comparisons = pairwise_comparisons,
  hide_ns = TRUE,
  step_increase = 0.12,
  y_expand = 0.3,
  palette = group_colors,
  xlab = "Group",
  ylab = "C1QBP Expression",
  title = "C1QBP Expression by Cell Type",
  theme_use = theme_my_stat,
  base_size = 10
) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# 分面图尺寸根据细胞类型数调整
n_ct <- length(cell_types)
ncol_facet <- ceiling(sqrt(n_ct))
nrow_facet <- ceiling(n_ct / ncol_facet)
facet_w <- ncol_facet * 3.5
facet_h <- nrow_facet * 3.5

save_plot_all(
  p_facet,
  "C1QBP_celltype_facet_compare",
  width = facet_w,
  height = facet_h
)

# ==============================================================================
# 7. 图3：各细胞类型独立图（逐一保存）
# ==============================================================================
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct) # 文件名安全化

  sub_data <- meta %>% filter(cell.type == ct)

  # 跳过某组细胞数不足的情况
  group_counts <- table(sub_data$group)
  if (any(group_counts < 3)) {
    message("跳过 ", ct, "：某组细胞数不足3个")
    next
  }

  p_ct <- compare_plot(
    data = sub_data,
    value.var = "C1QBP",
    group.by = "group",
    plot_type = "violin",
    violin_alpha = 0.6, # 更透明
    violin_linewidth = 1.0, # 更粗的边框
    add_box = TRUE,
    add_stat = "wilcox.test",
    stat_label = "p.signif",
    comparisons = pairwise_comparisons,
    hide_ns = TRUE,
    step_increase = 0.1,
    y_expand = 0.25,
    palette = group_colors,
    xlab = "Group",
    ylab = "C1QBP Expression",
    title = paste0("C1QBP — ", ct),
    theme_use = theme_my_stat,
    base_size = 12
  ) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  save_plot_all(
    p_ct,
    paste0("C1QBP_", ct_safe, "_group_compare"),
    width = 5.5,
    height = 5
  )
}

message("\n全部分析完成！输出目录：")
message(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/C1QBP/figures/"
)
