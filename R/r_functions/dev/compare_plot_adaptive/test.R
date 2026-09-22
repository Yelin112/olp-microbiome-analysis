# ============================================================================
# 测试文件：compare_plot_adaptive v3（管线式架构）
# 创建日期：2025-04-05
# 更新日期：2026-06-26（v3 重构）
#
# 用法：
#   在 R 中 source() 本文件，将自动运行全部测试并生成图片到 test_output/
#   或逐段运行以调试单个测试。
# ============================================================================

library(ggplot2)
library(dplyr)

# ---- 路径设置 ----

FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/R/r_functions/dev/compare_plot_adaptive"
UTILS_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/R/utils"

# 加载工具模块
source(file.path(UTILS_DIR, "helpers.R")) # %||%, %ni%
source(file.path(UTILS_DIR, "palette_system.R")) # get_colors()

# 加载策略模块
source(file.path(FUNC_DIR, "strategies/detect_strategy.R"))
source(file.path(FUNC_DIR, "strategies/layers_pure_scatter.R"))
source(file.path(FUNC_DIR, "strategies/layers_bar_points.R"))
source(file.path(FUNC_DIR, "strategies/layers_boxplot_points.R"))
source(file.path(FUNC_DIR, "strategies/layers_pure_boxplot.R"))

# 加载主函数
source(file.path(FUNC_DIR, "function.R"))

# 创建测试输出目录
TEST_OUT <- file.path(FUNC_DIR, "test_output")
if (!dir.exists(TEST_OUT)) {
  dir.create(TEST_OUT, recursive = TRUE)
}

# ---- 测试框架 ----

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

set.seed(42)

# ---- 测试数据 ----

data_tiny <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 3),
  value = c(rnorm(3, 35, 3), rnorm(3, 38, 4), rnorm(3, 25, 2), rnorm(3, 30, 3))
)

data_small <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 7),
  value = c(rnorm(7, 35, 5), rnorm(7, 38, 6), rnorm(7, 25, 4), rnorm(7, 30, 5))
)

data_medium <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 15),
  value = c(
    rnorm(15, 35, 5),
    rnorm(15, 38, 6),
    rnorm(15, 25, 4),
    rnorm(15, 30, 5)
  )
)

data_large <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 30),
  value = c(
    rnorm(30, 35, 5),
    rnorm(30, 38, 6),
    rnorm(30, 25, 4),
    rnorm(30, 30, 5)
  )
)

data_mixed <- data.frame(
  group = c(rep("A", 3), rep("B", 8), rep("C", 15), rep("D", 30)),
  value = c(
    rnorm(3, 35, 3),
    rnorm(8, 38, 4),
    rnorm(15, 25, 3),
    rnorm(30, 30, 5)
  )
)

# 带分面的数据
data_facet <- data.frame(
  group = rep(rep(c("A", "B"), each = 15), 2),
  sex = rep(c("Male", "Female"), each = 30),
  value = c(rnorm(30, 35, 5), rnorm(30, 32, 5))
)

cat("\n")
cat("============================================================\n")
cat("  compare_plot_adaptive v3 测试套件\n")
cat("============================================================\n")

# ============================================================================
# 一、策略自动检测
# ============================================================================
cat("\n── 一、策略自动检测 ──\n")

run_test("测试 1: n=3 → pure_scatter", {
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

run_test("测试 2: n=7 → bar_points", {
  result <- compare_plot_adaptive(data_small, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "bar_points")
  ggsave(
    file.path(TEST_OUT, "test2_bar_points.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 3: n=15 → boxplot_points", {
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

run_test("测试 4: n=30 → pure_boxplot", {
  result <- compare_plot_adaptive(data_large, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "pure_boxplot")
  ggsave(
    file.path(TEST_OUT, "test4_pure_boxplot.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 5: 混合样本量 max_n=30 → pure_boxplot", {
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
# 二、边界值
# ============================================================================
cat("\n── 二、边界值测试 ──\n")

run_test("测试 6a: n=5 → bar_points（边界）", {
  d <- data.frame(group = rep("A", 5), value = rnorm(5))
  result <- compare_plot_adaptive(d, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "bar_points")
})

run_test("测试 6b: n=10 → boxplot_points（边界）", {
  d <- data.frame(group = rep("A", 10), value = rnorm(10))
  result <- compare_plot_adaptive(d, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "boxplot_points")
})

run_test("测试 6c: n=20 → pure_boxplot（边界）", {
  d <- data.frame(group = rep("A", 20), value = rnorm(20))
  result <- compare_plot_adaptive(d, "group", "value", verbose = FALSE)
  stopifnot(result$strategy == "pure_boxplot")
})

# ============================================================================
# 三、手动策略覆盖
# ============================================================================
cat("\n── 三、手动策略覆盖 ──\n")

run_test("测试 7: 手动指定 bar_points（覆盖大样本 auto）", {
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
# 四、缺失值处理
# ============================================================================
cat("\n── 四、缺失值处理 ──\n")

run_test("测试 8: 缺失值自动移除并正常运行", {
  data_na <- data.frame(
    group = rep(c("A", "B"), each = 5),
    value = c(rnorm(4), NA, rnorm(4), NA)
  )
  result <- compare_plot_adaptive(data_na, "group", "value", verbose = FALSE)
  stopifnot(inherits(result$plot, "ggplot"))
})

# ============================================================================
# 五、自定义阈值
# ============================================================================
cat("\n── 五、自定义阈值 ──\n")

run_test("测试 9: 自定义阈值 (small=8, medium=15, large=25)", {
  # n=7 在默认下是 bar_points，改为 small=8 后应为 pure_scatter
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
# 六、错误处理
# ============================================================================
cat("\n── 六、错误处理 ──\n")

run_test("测试 10a: 不存在的列名 → error", {
  err <- tryCatch(
    compare_plot_adaptive(data_small, "nonexistent", "value"),
    error = function(e) "caught"
  )
  stopifnot(err == "caught")
})

run_test("测试 10b: 无效策略名 → error", {
  err <- tryCatch(
    compare_plot_adaptive(
      data_small,
      "group",
      "value",
      strategy = "invalid_strategy"
    ),
    error = function(e) "caught"
  )
  stopifnot(err == "caught")
})

run_test("测试 10c: 非 data.frame 输入 → error", {
  err <- tryCatch(
    compare_plot_adaptive(list(a = 1), "group", "value"),
    error = function(e) "caught"
  )
  stopifnot(err == "caught")
})

run_test("测试 10d: 全空数据 → error", {
  d <- data.frame(group = c("A", "B"), value = c(NA, NA))
  err <- tryCatch(
    compare_plot_adaptive(d, "group", "value", verbose = FALSE),
    error = function(e) "caught"
  )
  stopifnot(err == "caught")
})

# ============================================================================
# 七、图形变体（violin / show_points）
# ============================================================================
cat("\n── 七、图形变体 ──\n")

run_test("测试 11: boxplot_points + use_violin=TRUE", {
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

run_test("测试 12: pure_boxplot + show_points=TRUE", {
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
# 八、统计检验（v3 新功能）
# ============================================================================
cat("\n── 八、统计检验（v3 新增）──\n")

run_test("测试 13: t.test 两两比较", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    add_stat = "t.test",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test13_ttest.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 14: wilcox.test + 手动指定比较对", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    add_stat = "wilcox.test",
    comparisons = list(c("A", "B"), c("C", "D")),
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test14_wilcox.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 15: ANOVA 整体检验", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    add_stat = "anova",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test15_anova.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 16: kruskal.test 整体检验", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    add_stat = "kruskal.test",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test16_kruskal.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 17: 隐藏非显著 + 调整间距", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    add_stat = "t.test",
    hide_ns = TRUE,
    step_increase = 0.15,
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test17_hide_ns.png"),
    result$plot,
    width = 6,
    height = 5
  )
})

# ============================================================================
# 九、分面（v3 新功能）
# ============================================================================
cat("\n── 九、分面（v3 新增）──\n")

run_test("测试 18: split.by 分面", {
  result <- compare_plot_adaptive(
    data_facet,
    "group",
    "value",
    split.by = "sex",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test18_facet.png"),
    result$plot,
    width = 8,
    height = 4
  )
})

run_test("测试 19: 分面 + 统计检验", {
  result <- compare_plot_adaptive(
    data_facet,
    "group",
    "value",
    split.by = "sex",
    add_stat = "t.test",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test19_facet_stat.png"),
    result$plot,
    width = 8,
    height = 4
  )
})

# ============================================================================
# 十、外观定制（v3 新功能）
# ============================================================================
cat("\n── 十、外观定制（v3 新增）──\n")

run_test("测试 20: 斑马纹背景", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    add_bg = TRUE,
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test20_zebra.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 21: 自定义主题 theme_bw", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    theme_use = theme_bw,
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test21_theme_bw.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 22: 自定义配色 + 标题", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
    title = "Custom Palette + Title",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test22_custom_palette.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

run_test("测试 23: 策略=auto + 全部新功能一次性组合", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    strategy = "auto",
    add_stat = "wilcox.test",
    add_bg = TRUE,
    palette = "JCO",
    theme_use = theme_minimal,
    title = "Combined Features",
    xlab = "Group",
    ylab = "Measurement",
    verbose = FALSE
  )
  stopifnot(inherits(result$plot, "ggplot"))
  ggsave(
    file.path(TEST_OUT, "test23_combined.png"),
    result$plot,
    width = 6,
    height = 4
  )
})

# ============================================================================
# 十一、返回值的 metadata
# ============================================================================
cat("\n── 十一、返回值验证 ──\n")

run_test("测试 24: 返回列表结构完整", {
  result <- compare_plot_adaptive(
    data_medium,
    "group",
    "value",
    verbose = FALSE
  )
  stopifnot(is.list(result))
  stopifnot(c("plot", "strategy", "n_per_group", "metadata") %in% names(result))
  stopifnot(inherits(result$plot, "ggplot"))
  stopifnot(is.character(result$strategy))
  stopifnot(is.numeric(result$n_per_group))
  stopifnot(
    c("max_n", "median_n", "min_n", "thresholds") %in%
      names(result$metadata)
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
cat("  测试图片保存至: ", TEST_OUT, "\n")
cat("============================================================\n")
