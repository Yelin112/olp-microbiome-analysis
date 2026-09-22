# ============================================================================
# 使用示例：compare_plot_adaptive v3（管线式架构）
# 展示从小样本到大样本的各种典型使用场景，含 v3 新增的统计检验/分面等功能
# ============================================================================

library(ggplot2)
library(dplyr)

# 加载函数
FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/R/r_functions/dev/compare_plot_adaptive"
source(file.path(FUNC_DIR, "function.R"))

set.seed(2025)

# ============================================================================
# 示例1：基础用法 - 小样本实验（自动纯散点）
# 场景：细胞实验，每组3个生物学重复
# ============================================================================
cat("示例1: 小样本细胞实验数据\n")

data_cell <- data.frame(
  group      = rep(c("Control", "Treatment_A", "Treatment_B"), each = 3),
  expression = c(5.2, 5.5, 5.1, 7.3, 7.8, 7.5, 4.1, 4.3, 4.2)
)

result1 <- compare_plot_adaptive(
  data      = data_cell,
  group_col = "group",
  value_col = "expression",
  ylab      = "Gene Expression (log2 CPM)",
  title     = "Gene Expression - Small Experiment"
)

cat("  → 策略:", result1$strategy,
    "| 各组 n:", result1$n_per_group, "\n\n")
print(result1$plot)

# ============================================================================
# 示例2：中等样本 - 动物实验（自动柱状图+散点）+ 统计检验
# 场景：n=8，每笼8只小鼠，带 t.test
# ============================================================================
cat("示例2: 中等样本 + t.test\n")

data_mouse <- data.frame(
  group  = rep(c("WT", "KO", "Rescue"), each = 8),
  weight = c(rnorm(8, 25, 1.5), rnorm(8, 22, 2), rnorm(8, 24, 1.8))
)

result2 <- compare_plot_adaptive(
  data      = data_mouse,
  group_col = "group",
  value_col = "weight",
  palette   = "NEJM",
  add_stat  = "t.test",          # ← v3 新功能
  ylab      = "Body Weight (g)",
  title     = "Body Weight - Mouse Experiment"
)

cat("  → 策略:", result2$strategy, "\n\n")
print(result2$plot)

# ============================================================================
# 示例3：大样本 - 流式细胞仪数据（自动箱线图）+ 分面
# 场景：30个样本，按性别分面
# ============================================================================
cat("示例3: 大样本 + 分面\n")

data_flow <- data.frame(
  cell_type  = rep(rep(c("CD4+", "CD8+"), each = 15), 2),
  sex        = rep(c("Male", "Female"), each = 30),
  percentage = c(rnorm(30, 40, 8), rnorm(30, 25, 6))
)

result3 <- compare_plot_adaptive(
  data      = data_flow,
  group_col = "cell_type",
  value_col = "percentage",
  palette   = "JCO",
  split.by  = "sex",             # ← v3 新功能
  add_stat  = "wilcox.test",     # ← v3 新功能
  ylab      = "Percentage (%)",
  title     = "Cell Type Distribution by Sex"
)

cat("  → 策略:", result3$strategy,
    "| max_n:", result3$metadata$max_n, "\n\n")
print(result3$plot)

# ============================================================================
# 示例4：手动策略 + 斑马纹背景 + 自定义主题
# ============================================================================
cat("示例4: 手动策略 + 斑马纹 + theme_bw\n")

result4 <- compare_plot_adaptive(
  data      = data_mouse,
  group_col = "group",
  value_col = "weight",
  strategy  = "boxplot_points",
  add_bg    = TRUE,              # ← v3 新功能
  theme_use = theme_bw,          # ← v3 新功能
  palette   = "Lancet",
  title     = "Boxplot with Zebra Background"
)

cat("  → 策略:", result4$strategy, "(手动指定)\n\n")
print(result4$plot)

# ============================================================================
# 示例5：自定义阈值 + ANOVA
# ============================================================================
cat("示例5: 自定义阈值 + ANOVA 整体检验\n")

data_microbiome <- data.frame(
  group     = rep(c("Healthy", "Disease", "Recovery"), each = 12),
  abundance = c(rlnorm(12, 2, 0.5), rlnorm(12, 1.5, 0.6), rlnorm(12, 1.8, 0.4))
)

result5 <- compare_plot_adaptive(
  data             = data_microbiome,
  group_col        = "group",
  value_col        = "abundance",
  threshold_small  = 5,
  threshold_medium = 15,
  threshold_large  = 30,
  add_stat         = "anova",    # ← v3 新功能
  palette          = "AAAS",
  ylab             = "Relative Abundance",
  title            = "Microbiome Abundance (ANOVA)"
)

cat("  → 策略:", result5$strategy,
    "(n=12, threshold_medium=15)\n\n")
print(result5$plot)

# ============================================================================
# 示例6：在 ggplot 对象上继续定制
# ============================================================================
cat("示例6: 后处理定制\n")

result6 <- compare_plot_adaptive(
  data      = data_cell,
  group_col = "group",
  value_col = "expression",
  palette   = c("#2196F3", "#FF5722", "#4CAF50"),
  verbose   = FALSE
)

p_final <- result6$plot +
  labs(
    title    = "Gene Expression Comparison",
    subtitle = paste("Strategy:", result6$strategy,
                     "| n:", paste(result6$n_per_group, collapse = "/")),
    caption  = "Mean ± SEM"
  ) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "gray50"),
    legend.position  = "none"
  )

print(p_final)

cat("\n============================================================\n")
cat("  所有示例完成。v3 新功能演示：\n")
cat("    - add_stat (t.test / wilcox.test / anova)\n")
cat("    - split.by (facet_wrap 分面)\n")
cat("    - add_bg (斑马纹背景)\n")
cat("    - theme_use (自定义主题)\n")
cat("============================================================\n")
