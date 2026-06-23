# ============================================================================
# 测试脚本：compare_plot() 自适应策略
# 用法：在 RStudio 中直接运行，或 source("test_strategy.R")
# ============================================================================

library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

source("Scripts//compare_plot_optimized.R")

set.seed(42)

# 输出目录
if (!dir.exists("test_output")) {
  dir.create("test_output")
}

pass <- 0
fail <- 0

chk <- function(name, expr) {
  cat("▶", name, "... ")
  tryCatch(
    {
      expr
      cat("✅\n")
      pass <<- pass + 1
    },
    error = function(e) {
      cat("❌", conditionMessage(e), "\n")
      fail <<- fail + 1
    }
  )
}

cat("============================================================\n")
cat("  compare_plot() 自适应策略测试\n")
cat("============================================================\n\n")

# ── 数据准备 ──────────────────────────────────────────────────
d_tiny <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 3),
  value = rnorm(12, rep(c(30, 35, 28, 32), each = 3), 2)
)
d_small <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 7),
  value = rnorm(28, rep(c(30, 35, 28, 32), each = 7), 4)
)
d_medium <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 15),
  value = rnorm(60, rep(c(30, 35, 28, 32), each = 15), 5)
)
d_large <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 30),
  value = rnorm(120, rep(c(30, 35, 28, 32), each = 30), 5)
)

cat("--- 自动策略检测 ---\n")

chk("n=3 → pure_scatter (scatter + 散点层)", {
  p <- compare_plot(d_tiny, value.var = "value", group.by = "group")
  # scatter 策略：crossbar + errorbar + jitter = 3 层
  stopifnot(length(p$layers) >= 3)
  ggsave("test_output/s1_pure_scatter.png", p, width = 5, height = 4)
})

chk("n=7 → bar_points (柱+误差线+散点)", {
  p <- compare_plot(
    d_small,
    value.var = "value",
    group.by = "group",
    add_point = T
  )
  # bar + errorbar + jitter = 3 层
  stopifnot(length(p$layers) >= 3)
  ggsave("test_output/s2_bar_points.png", p, width = 5, height = 4)
})

chk("n=15 → boxplot_points (箱线图+散点)", {
  p <- compare_plot(d_medium, value.var = "value", group.by = "group")
  # box + mean_point + jitter = 3 层
  stopifnot(length(p$layers) >= 3)
  ggsave("test_output/s3_boxplot_points.png", p, width = 5, height = 4)
})

chk("n=30 → pure_boxplot (仅箱线图)", {
  p <- compare_plot(d_large, value.var = "value", group.by = "group")
  # box + mean_point = 2 层（无散点）
  stopifnot(length(p$layers) >= 2)
  ggsave("test_output/s4_pure_boxplot.png", p, width = 5, height = 4)
})

cat("\n--- 边界值 ---\n")

chk("n=5 → bar_points", {
  d <- data.frame(group = rep(c("A", "B"), each = 5), value = rnorm(10))
  compare_plot(d, value.var = "value", group.by = "group")
})

chk("n=10 → boxplot_points", {
  d <- data.frame(group = rep(c("A", "B"), each = 10), value = rnorm(20))
  compare_plot(d, value.var = "value", group.by = "group")
})

chk("n=20 → pure_boxplot", {
  d <- data.frame(group = rep(c("A", "B"), each = 20), value = rnorm(40))
  compare_plot(d, value.var = "value", group.by = "group")
})

cat("\n--- 手动策略覆盖 ---\n")

chk("strategy='none' + plot_type='violin' (旧行为)", {
  p <- compare_plot(
    d_large,
    value.var = "value",
    group.by = "group",
    strategy = "none",
    plot_type = "violin"
  )
  ggsave("test_output/s5_manual_violin.png", p, width = 5, height = 4)
})

chk("strategy='bar_points' 强制覆盖大样本", {
  p <- compare_plot(
    d_large,
    value.var = "value",
    group.by = "group",
    strategy = "bar_points"
  )
  ggsave("test_output/s6_forced_bar.png", p, width = 5, height = 4)
})

chk("自定义阈值 threshold_medium=20 (n=15 → bar_points)", {
  p <- compare_plot(
    d_medium,
    value.var = "value",
    group.by = "group",
    threshold_medium = 20
  )
  # 应为 bar_points 而非 boxplot_points
})

cat("\n--- 原有功能兼容性 ---\n")

chk("strategy='none' + 统计检验正常", {
  p <- compare_plot(
    d_small,
    value.var = "value",
    group.by = "group",
    strategy = "none",
    plot_type = "box",
    add_stat = "wilcox.test",
    stat_label = "p.signif"
  )
  ggsave("test_output/s7_stat_test.png", p, width = 5, height = 4)
})

chk("strategy='none' + add_bg + split.by 正常", {
  d <- cbind(d_small, batch = rep(c("X", "Y"), 14))
  p <- compare_plot(
    d,
    value.var = "value",
    group.by = "group",
    split.by = "batch",
    add_bg = TRUE,
    strategy = "none",
    plot_type = "box"
  )
  ggsave("test_output/s8_facet_bg.png", p, width = 8, height = 4)
})

chk("自动策略 + 统计检验叠加", {
  p <- compare_plot(
    d_small,
    value.var = "value",
    group.by = "group",
    add_stat = "wilcox.test"
  )
  ggsave("test_output/s9_auto_stat.png", p, width = 5, height = 4)
})

# ── 汇总 ─────────────────────────────────────────────────────
total <- pass + fail
cat("\n============================================================\n")
cat(sprintf("  结果: %d/%d 通过", pass, total))
if (fail > 0) {
  cat(sprintf("，%d 个失败", fail))
}
cat("\n  图片已保存至 test_output/\n")
cat("============================================================\n")
