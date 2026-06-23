# ============================================================================
# 使用示例：compare_plot_adaptive
# 展示从小样本到大样本的各种典型使用场景
# ============================================================================

library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 加载函数
FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/compare_plot_adaptive"
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

cat("  → 使用策略:", result1$strategy, "\n")
cat("  → 各组 n:", result1$n_per_group, "\n\n")
print(result1$plot)

# ============================================================================
# 示例2：中等样本 - 动物实验（自动柱状图+散点）
# 场景：n=8，每笼8只小鼠
# ============================================================================
cat("示例2: 中等样本动物实验\n")

data_mouse <- data.frame(
  group  = rep(c("WT", "KO", "Rescue"), each = 8),
  weight = c(
    rnorm(8, 25, 1.5),
    rnorm(8, 22, 2),
    rnorm(8, 24, 1.8)
  )
)

result2 <- compare_plot_adaptive(
  data      = data_mouse,
  group_col = "group",
  value_col = "weight",
  palette   = "NEJM",
  ylab      = "Body Weight (g)",
  title     = "Body Weight - Mouse Experiment"
)

cat("  → 使用策略:", result2$strategy, "\n\n")
print(result2$plot)

# ============================================================================
# 示例3：大样本 - 流式细胞仪数据（自动箱线图）
# 场景：30个患者样本
# ============================================================================
cat("示例3: 大样本流式细胞仪数据\n")

data_flow <- data.frame(
  cell_type  = rep(c("CD4+", "CD8+", "NK", "B cell"), each = 30),
  percentage = c(
    rnorm(30, 40, 8),
    rnorm(30, 25, 6),
    rnorm(30, 15, 4),
    rnorm(30, 20, 5)
  )
)

result3 <- compare_plot_adaptive(
  data      = data_flow,
  group_col = "cell_type",
  value_col = "percentage",
  palette   = "JCO",
  ylab      = "Percentage (%)",
  title     = "Cell Type Distribution"
)

cat("  → 使用策略:", result3$strategy, "\n\n")
print(result3$plot)

# ============================================================================
# 示例4：手动指定策略（强制使用柱状图+散点）
# 场景：虽然 n=30，但希望展示原始数据
# ============================================================================
cat("示例4: 手动指定策略 bar_points\n")

result4 <- compare_plot_adaptive(
  data      = data_flow,
  group_col = "cell_type",
  value_col = "percentage",
  strategy  = "bar_points",
  palette   = "Lancet",
  ylab      = "Percentage (%)",
  title     = "Cell Distribution (Forced Bar+Points)"
)

cat("  → 使用策略:", result4$strategy, "(手动指定)\n\n")
print(result4$plot)

# ============================================================================
# 示例5：自定义阈值
# 场景：对于微生物组数据，n=12 时仍希望用柱状图+散点
# ============================================================================
cat("示例5: 自定义策略切换阈值\n")

data_microbiome <- data.frame(
  group     = rep(c("Healthy", "Disease", "Recovery"), each = 12),
  abundance = c(
    rlnorm(12, 2, 0.5),
    rlnorm(12, 1.5, 0.6),
    rlnorm(12, 1.8, 0.4)
  )
)

result5 <- compare_plot_adaptive(
  data             = data_microbiome,
  group_col        = "group",
  value_col        = "abundance",
  threshold_small  = 5,
  threshold_medium = 15,   # 提高至15，n=12 仍用 bar_points
  threshold_large  = 30,
  palette          = "AAAS",
  ylab             = "Relative Abundance",
  title            = "Microbiome Abundance (Custom Thresholds)"
)

cat("  → 使用策略:", result5$strategy, "(n=12, threshold_medium=15)\n\n")
print(result5$plot)

# ============================================================================
# 示例6：在结果基础上继续定制 ggplot
# 场景：添加统计标注、修改主题等
# ============================================================================
cat("示例6: 在返回的 ggplot 对象上继续定制\n")

result6 <- compare_plot_adaptive(
  data      = data_cell,
  group_col = "group",
  value_col = "expression",
  palette   = c("#2196F3", "#FF5722", "#4CAF50"),  # 自定义颜色
  verbose   = FALSE
)

# 在返回的 ggplot 对象上叠加额外的主题设置
p_final <- result6$plot +
  labs(
    title    = "Gene Expression Comparison",
    subtitle = paste("Strategy:", result6$strategy,
                     "| n per group:", paste(result6$n_per_group, collapse = "/")),
    caption  = "Mean ± SEM shown"
  ) +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray50"),
    plot.caption  = element_text(size = 8, color = "gray60"),
    legend.position = "none"
  )

print(p_final)

# 保存高质量图片
# ggsave("final_expression_plot.png", p_final, width = 6, height = 5, dpi = 300)

cat("\n============================================================\n")
cat("  所有示例完成\n")
cat("  返回对象包含: $plot, $strategy, $n_per_group, $metadata\n")
cat("============================================================\n")
