# ============================================================================
# 测试文件：compare_plot_adaptive
# 创建日期：2025-04-05
# 用法：在 compare_plot_adaptive/ 目录下运行，或设置好工作目录后运行
# ============================================================================

library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 设置工作目录到函数目录（根据实际路径修改）
# setwd("e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/compare_plot_adaptive")

# 加载主函数（会自动 source 所有策略模块）
FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/compare_plot_adaptive"
source(file.path(FUNC_DIR, "strategies/detect_strategy.R"))
source(file.path(FUNC_DIR, "strategies/pure_scatter.R"))
source(file.path(FUNC_DIR, "strategies/bar_points.R"))
source(file.path(FUNC_DIR, "strategies/boxplot_points.R"))
source(file.path(FUNC_DIR, "strategies/pure_boxplot.R"))

# 手动定义辅助函数（function.R 中的内部函数）
`%||%` <- function(a, b) if (is.null(a)) b else a
.get_colors <- function(palette, n_groups) {
  sci_palettes <- list(
    "NPG" = c(
      "#E64B35",
      "#4DBBD5",
      "#00A087",
      "#3C5488",
      "#F39B7F",
      "#8491B4",
      "#91D1C2",
      "#DC0000",
      "#7E6148"
    )
  )
  if (length(palette) > 1) {
    return(palette)
  }
  if (palette %in% names(sci_palettes)) {
    base_cols <- sci_palettes[[palette]]
  } else {
    base_cols <- sci_palettes[["NPG"]]
  }
  if (n_groups > length(base_cols)) {
    colorRampPalette(base_cols)(n_groups)
  } else {
    base_cols[seq_len(n_groups)]
  }
}

# 加载主函数
source(file.path(FUNC_DIR, "function.R"))

# 创建测试输出目录
TEST_OUT <- file.path(FUNC_DIR, "test_output")
if (!dir.exists(TEST_OUT)) {
  dir.create(TEST_OUT, recursive = TRUE)
}

cat("============================================================\n")
cat("  compare_plot_adaptive 测试套件\n")
cat("============================================================\n\n")

pass_count <- 0
fail_count <- 0

run_test <- function(name, expr) {
  cat("▶ ", name, "... ")
  result <- tryCatch(
    {
      expr
      cat("✅ 通过\n")
      pass_count <<- pass_count + 1
    },
    error = function(e) {
      cat("❌ 失败:", conditionMessage(e), "\n")
      fail_count <<- fail_count + 1
    }
  )
}

set.seed(42) # 保证可重复性

# ============================================================================
# 测试1：极小样本 (n=3) -> pure_scatter
# ============================================================================
cat("\n--- 策略自动检测测试 ---\n")

data_tiny <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 3),
  value = c(rnorm(3, 35, 3), rnorm(3, 38, 4), rnorm(3, 25, 2), rnorm(3, 30, 3))
)

run_test("测试1: n=3 → pure_scatter", {
  result <- compare_plot_adaptive(data_tiny, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "pure_scatter")
  stopifnot(inherits(result$plot, "ggplot"))
  stopifnot(length(result$n_per_group) == 4)
  ggsave(
    file.path(TEST_OUT, "test1_pure_scatter.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 测试2：小样本 (n=7) -> bar_points
# ============================================================================

data_small <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 7),
  value = c(rnorm(7, 35, 5), rnorm(7, 38, 6), rnorm(7, 25, 4), rnorm(7, 30, 5))
)

run_test("测试2: n=7 → bar_points", {
  result <- compare_plot_adaptive(data_small, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "bar_points")
  ggsave(
    file.path(TEST_OUT, "test2_bar_points.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 测试3：中等样本 (n=15) -> boxplot_points
# ============================================================================

data_medium <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 15),
  value = c(
    rnorm(15, 35, 5),
    rnorm(15, 38, 6),
    rnorm(15, 25, 4),
    rnorm(15, 30, 5)
  )
)

run_test("测试3: n=15 → boxplot_points", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    verbose = FALSE
  )
  stopifnot(result$strategy == "boxplot_points")
  ggsave(
    file.path(TEST_OUT, "test3_boxplot_points.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 测试4：大样本 (n=30) -> pure_boxplot
# ============================================================================

data_large <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 30),
  value = c(
    rnorm(30, 35, 5),
    rnorm(30, 38, 6),
    rnorm(30, 25, 4),
    rnorm(30, 30, 5)
  )
)

run_test("测试4: n=30 → pure_boxplot", {
  result <- compare_plot_adaptive(data_large, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "pure_boxplot")
  ggsave(
    file.path(TEST_OUT, "test4_pure_boxplot.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 测试5：混合样本量（以 max_n 决定策略）
# ============================================================================

data_mixed <- data.frame(
  group = c(rep("A", 3), rep("B", 8), rep("C", 15), rep("D", 30)),
  value = c(
    rnorm(3, 35, 3),
    rnorm(8, 38, 4),
    rnorm(15, 25, 3),
    rnorm(30, 30, 5)
  )
)

run_test("测试5: 混合样本 max_n=30 → pure_boxplot", {
  result <- compare_plot_adaptive(data_mixed, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "pure_boxplot")
  stopifnot(result$metadata$max_n == 30)
  stopifnot(result$metadata$min_n == 3)
  ggsave(
    file.path(TEST_OUT, "test5_mixed.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 测试6：边界值（n = 5, 10, 20）
# ============================================================================
cat("\n--- 边界值测试 ---\n")

run_test("测试6a: n=5 → bar_points (边界)", {
  d <- data.frame(group = rep("A", 5), value = rnorm(5))
  result <- compare_plot_adaptive(d, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "bar_points")
})

run_test("测试6b: n=10 → boxplot_points (边界)", {
  d <- data.frame(group = rep("A", 10), value = rnorm(10))
  result <- compare_plot_adaptive(d, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "boxplot_points")
})

run_test("测试6c: n=20 → pure_boxplot (边界)", {
  d <- data.frame(group = rep("A", 20), value = rnorm(20))
  result <- compare_plot_adaptive(d, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "pure_boxplot")
})

# ============================================================================
# 测试7：手动指定策略（覆盖自动检测）
# ============================================================================
cat("\n--- 手动策略测试 ---\n")

run_test("测试7: 手动指定 bar_points（覆盖大样本）", {
  result <- compare_plot_adaptive(
    data_large,
    "group",
    "value",
    strategy = "bar_points",
    verbose = FALSE
  )
  stopifnot(result$strategy == "bar_points")
  ggsave(
    file.path(TEST_OUT, "test7_manual_override.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 测试8：缺失值处理
# ============================================================================
cat("\n--- 缺失值处理测试 ---\n")

data_na <- data.frame(
  group = rep(c("A", "B"), each = 5),
  value = c(rnorm(4), NA, rnorm(4), NA)
)

run_test("测试8: 缺失值自动移除", {
  result <- compare_plot_adaptive(data_na, "group", "value", verbose = TRUE)
  # 应该能正常运行（NA 被移除后 n=4，策略为 pure_scatter）
  stopifnot(inherits(result$plot, "ggplot"))
})

# ============================================================================
# 测试9：自定义阈值
# ============================================================================
cat("\n--- 自定义阈值测试 ---\n")

run_test("测试9: 自定义阈值 (small=8, medium=15, large=25)", {
  # n=7 在默认下是 bar_points，在 small=8 时应为 pure_scatter
  result <- compare_plot_adaptive(
    data_small,
    "group",
    "value",
    threshold_small = 8,
    threshold_medium = 15,
    threshold_large = 25,
    verbose = FALSE
  )
  stopifnot(result$strategy == "pure_scatter")
})

# ============================================================================
# 测试10：错误处理
# ============================================================================
cat("\n--- 错误处理测试 ---\n")

run_test("测试10a: 不存在的列名", {
  result <- tryCatch(
    compare_plot_adaptive(data_small, "nonexistent", "value"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

run_test("测试10b: 无效的策略名称", {
  result <- tryCatch(
    compare_plot_adaptive(
      data_small,
      "group",
      "value",
      strategy = "invalid_strategy"
    ),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

run_test("测试10c: 非 data.frame 输入", {
  result <- tryCatch(
    compare_plot_adaptive(list(a = 1), "group", "value"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

# ============================================================================
# 测试11：小提琴图模式
# ============================================================================
cat("\n--- 图形变体测试 ---\n")

run_test("测试11: boxplot_points + use_violin=TRUE", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    use_violin = TRUE,
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test11_violin_points.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试12: pure_boxplot + show_points=TRUE", {
  result <- compare_plot_adaptive(
    data_large,
    "group",
    "value",
    show_points = TRUE,
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test12_boxplot_with_points.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 汇总
# ============================================================================
total <- pass_count + fail_count
cat("\n============================================================\n")
cat(sprintf("  测试结果: %d/%d 通过", pass_count, total))
if (fail_count > 0) {
  cat(sprintf("  (%d 个失败)", fail_count))
}
cat("\n")
cat("  测试图片保存至:", TEST_OUT, "\n")
cat("============================================================\n")
