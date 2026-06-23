# ============================================================================
# 测试脚本：pcoa_plot()
# 用法：在 RStudio 中 source 本文件，或设置好工作目录后直接运行
# ============================================================================

library(ggplot2)

FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/pcoa_plot"
source(file.path(FUNC_DIR, "function.R"))

TEST_OUT <- file.path(FUNC_DIR, "test_output")
if (!dir.exists(TEST_OUT)) {
  dir.create(TEST_OUT, recursive = TRUE)
}

set.seed(42)

pass <- 0
fail <- 0

chk <- function(name, expr) {
  cat("▶", name, "... ")
  tryCatch(
    {
      force(expr)
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
cat("  pcoa_plot() 测试套件\n")
cat("============================================================\n\n")

# ── 测试数据 ──────────────────────────────────────────────────────────────────

# 两组数据（类似图C，Non_RW vs RW）
d2 <- data.frame(
  PC1 = c(rnorm(25, -0.2, 0.12), rnorm(25, 0.2, 0.15)),
  PC2 = c(rnorm(25, 0.02, 0.08), rnorm(25, 0.05, 0.10)),
  group = rep(c("Non_RW", "RW"), each = 25)
)

# 三组数据（类似右图，SCNC_A / _N / _Y）
d3 <- data.frame(
  PCoA1 = c(
    rnorm(40, -0.15, 0.10),
    rnorm(40, 0.05, 0.12),
    rnorm(40, 0.18, 0.11)
  ),
  PCoA2 = c(
    rnorm(40, -0.05, 0.06),
    rnorm(40, 0.08, 0.08),
    rnorm(40, 0.01, 0.07)
  ),
  group = rep(c("SCNC_A", "SCNC_N", "SCNC_Y"), each = 40)
)

# ── Part 1：基本功能 ──────────────────────────────────────────────────────────
cat("--- 基本功能 ---\n")

chk("marginal='none' 返回 ggplot 对象", {
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none"
  )
  stopifnot(inherits(p, "ggplot"))
  ggsave(file.path(TEST_OUT, "t01_none.png"), p, width = 5, height = 4)
})

chk("marginal='none' 无椭圆", {
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none",
    add_ellipse = FALSE
  )
  stopifnot(inherits(p, "ggplot"))
  ggsave(file.path(TEST_OUT, "t02_no_ellipse.png"), p, width = 5, height = 4)
})

chk("stat_text 单行注释", {
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none",
    stat_text = "adonis R² = 0.33; P = 0.001"
  )
  stopifnot(inherits(p, "ggplot"))
  ggsave(file.path(TEST_OUT, "t03_stat_single.png"), p, width = 5, height = 4)
})

chk("stat_text 多行注释（斜体）", {
  p <- pcoa_plot(
    d3,
    x = "PCoA1",
    y = "PCoA2",
    group.by = "group",
    marginal = "none",
    stat_text = c(
      "Anosim: R = 0.477, P = 0.001",
      "Adonis: R² = 0.329, P = 0.001"
    ),
    stat_italic = TRUE
  )
  stopifnot(inherits(p, "ggplot"))
  ggsave(file.path(TEST_OUT, "t04_stat_multi.png"), p, width = 5, height = 4)
})

# ── Part 2：密度边际图 ────────────────────────────────────────────────────────
cat("\n--- 密度边际图（ggExtra） ---\n")

chk("marginal='density' 两组（复现图C风格）", {
  if (!requireNamespace("ggExtra", quietly = TRUE)) {
    stop("需要 ggExtra 包")
  }
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "density",
    stat_text = "adonis R² = 0.33; P = 0.001",
    palette = 'Paired'
  )
  ggsave(file.path(TEST_OUT, "t05_density_2grp.png"), p, width = 5, height = 5)
})

chk("marginal='density' 三组", {
  if (!requireNamespace("ggExtra", quietly = TRUE)) {
    stop("需要 ggExtra 包")
  }
  p <- pcoa_plot(
    d3,
    x = "PCoA1",
    y = "PCoA2",
    group.by = "group",
    marginal = "density"
  )
  ggsave(file.path(TEST_OUT, "t06_density_3grp.png"), p, width = 5, height = 5)
})

chk("density + add_ellipse=FALSE", {
  if (!requireNamespace("ggExtra", quietly = TRUE)) {
    stop("需要 ggExtra 包")
  }
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "density",
    add_ellipse = FALSE
  )
  ggsave(
    file.path(TEST_OUT, "t07_density_no_ellipse.png"),
    p,
    width = 5,
    height = 5
  )
})

# ── Part 3：箱线图边际图 ──────────────────────────────────────────────────────
cat("\n--- 箱线图边际图（aplot） ---\n")

chk("marginal='boxplot' 三组（复现右图风格）", {
  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("需要 aplot 包")
  }
  p <- pcoa_plot(
    d3,
    x = "PCoA1",
    y = "PCoA2",
    group.by = "group",
    marginal = "boxplot",
    stat_text = c(
      "Anosim: R = 0.477, P = 0.001",
      "Adonis: R² = 0.329, P = 0.001"
    ),
    xlab = "(PCoA1: 39.99%)",
    ylab = "(PCoA2: 17.69%)"
  )
  ggsave(file.path(TEST_OUT, "t08_boxplot_3grp.png"), p, width = 6, height = 5)
})

chk("boxplot + sig_letters 显著性字母", {
  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("需要 aplot 包")
  }
  p <- pcoa_plot(
    d3,
    x = "PCoA1",
    y = "PCoA2",
    group.by = "group",
    marginal = "boxplot",
    sig_letters = c(SCNC_A = "b", SCNC_N = "a", SCNC_Y = "b"),
    xlab = "(PCoA1: 39.99%)",
    ylab = "(PCoA2: 17.69%)"
  )
  ggsave(
    file.path(TEST_OUT, "t09_boxplot_letters.png"),
    p,
    width = 6,
    height = 5
  )
})

chk("boxplot 两组", {
  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("需要 aplot 包")
  }
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "boxplot"
  )
  ggsave(file.path(TEST_OUT, "t10_boxplot_2grp.png"), p, width = 6, height = 5)
})

# ── Part 4：外观自定义 ────────────────────────────────────────────────────────
cat("\n--- 外观自定义 ---\n")

chk("自定义颜色向量", {
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none",
    palette = c("#8B7CB3", "#E5A02A")
  )
  ggsave(file.path(TEST_OUT, "t11_custom_color.png"), p, width = 5, height = 4)
})

chk("shape=21 + stroke 参数", {
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none",
    point_shape = 21,
    point_size = 3,
    point_alpha = 0.9,
    color = "black",
    stroke = 0.5
  )
  ggsave(file.path(TEST_OUT, "t12_shape21.png"), p, width = 5, height = 4)
})

chk("自定义主题 theme_bw", {
  p <- pcoa_plot(
    d3,
    x = "PCoA1",
    y = "PCoA2",
    group.by = "group",
    marginal = "none",
    theme_use = theme_bw,
    title = "PCoA - theme_bw"
  )
  ggsave(file.path(TEST_OUT, "t13_theme_bw.png"), p, width = 5, height = 4)
})

chk("stat_text 左上角定位", {
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none",
    stat_text = "Adonis P = 0.001",
    stat_x = -Inf,
    stat_y = Inf,
    stat_hjust = -0.05,
    stat_vjust = 1.5
  )
  ggsave(file.path(TEST_OUT, "t14_stat_topleft.png"), p, width = 5, height = 4)
})

chk("marginal_size 自定义（0.4）", {
  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("需要 aplot 包")
  }
  p <- pcoa_plot(
    d3,
    x = "PCoA1",
    y = "PCoA2",
    group.by = "group",
    marginal = "boxplot",
    marginal_size = 0.4
  )
  ggsave(
    file.path(TEST_OUT, "t15_marginal_size04.png"),
    p,
    width = 6,
    height = 5
  )
})

# ── Part 5：错误处理 ──────────────────────────────────────────────────────────
cat("\n--- 错误处理 ---\n")

chk("非 data.frame 输入报错", {
  result <- tryCatch(
    pcoa_plot(list(a = 1), "PC1", "PC2", "group"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("不存在的列名报错", {
  result <- tryCatch(
    pcoa_plot(d2, "nonexist", "PC2", "group"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("非数值坐标列报错", {
  d_bad <- d2
  d_bad$PC1 <- as.character(d_bad$PC1)
  result <- tryCatch(
    pcoa_plot(d_bad, "PC1", "PC2", "group"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("含NA数据自动移除（不报错）", {
  d_na <- d2
  d_na$PC1[c(1, 5, 10)] <- NA
  p <- pcoa_plot(
    d_na,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "none"
  )
  stopifnot(inherits(p, "ggplot"))
})

chk("sig_letters 含不存在分组名（给出警告，不报错）", {
  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("需要 aplot 包")
  }
  p <- pcoa_plot(
    d2,
    x = "PC1",
    y = "PC2",
    group.by = "group",
    marginal = "boxplot",
    sig_letters = c(Non_RW = "b", RW = "a", Unknown = "c")
  )
  stopifnot(!is.null(p))
})

# ── 汇总 ──────────────────────────────────────────────────────────────────────
total <- pass + fail
cat("\n============================================================\n")
cat(sprintf("  结果：%d / %d 通过", pass, total))
if (fail > 0) {
  cat(sprintf("，%d 个失败", fail))
}
cat("\n  图片已保存至：", TEST_OUT, "\n")
cat("============================================================\n")
