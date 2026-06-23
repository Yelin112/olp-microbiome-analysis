# ============================================================================
# 测试脚本：beta_calc()
# 测试数据：vegan 内置 dune 数据集（无需外部文件）
# 用法：source("test.R")  或在 RStudio 中直接运行
# ============================================================================

library(vegan)

FUNC_DIR <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/beta_calc"
source(file.path(FUNC_DIR, "function.R"))

TEST_OUT <- file.path(FUNC_DIR, "test_output")
if (!dir.exists(TEST_OUT)) {
  dir.create(TEST_OUT, recursive = TRUE)
}

# ── 测试框架 ──────────────────────────────────────────────────────────────────
pass <- 0L
fail <- 0L

chk <- function(name, expr) {
  cat("▶", name, "... ")
  tryCatch(
    {
      force(expr)
      cat("✅\n")
      pass <<- pass + 1L
    },
    error = function(e) {
      cat("❌", conditionMessage(e), "\n")
      fail <<- fail + 1L
    }
  )
}

cat("============================================================\n")
cat("  beta_calc() 测试套件\n")
cat("============================================================\n\n")

# ── 测试数据准备 ──────────────────────────────────────────────────────────────
data(dune) # 20 sites × 30 species（样本×特征，已是正确方向）
data(dune.env) # 20 sites × 5 variables（Management, Moisture, Use, Manure, A1）
set.seed(42)

# 模拟转录组数据（40样本×100基因）
expr_mat <- matrix(
  rnorm(4000, mean = 10, sd = 3),
  nrow = 40,
  dimnames = list(paste0("S", 1:40), paste0("Gene", 1:100))
)
expr_meta <- data.frame(
  group = rep(c("Control", "TreatA", "TreatB", "TreatC"), each = 10),
  batch = rep(c("B1", "B2"), 20),
  row.names = rownames(expr_mat)
)

# ── Part 1: 基本功能 ──────────────────────────────────────────────────────────
cat("--- Part 1: 基本功能 ---\n")

chk("PCoA + Bray-Curtis（dune 数据集）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    method = "PCoA",
    dist = "bray",
    verbose = FALSE
  )
  stopifnot(is.list(res))
  stopifnot(all(
    c(
      "coords",
      "dist_mat",
      "stat_table",
      "axis_labs",
      "stat_text",
      "method",
      "dist_method"
    ) %in%
      names(res)
  ))
  stopifnot(nrow(res$coords) == nrow(dune))
  stopifnot("Management" %in% names(res$coords))
  stopifnot(res$method == "PCoA")
  cat("  axis_labs:", res$axis_labs, "\n  stat_text:", res$stat_text, "\n")
})

chk("NMDS + Bray-Curtis（返回 stress）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    method = "NMDS",
    dist = "bray",
    verbose = FALSE
  )
  stopifnot(res$method == "NMDS")
  stopifnot(!is.null(res$stress))
  stopifnot(res$stress > 0 && res$stress < 1)
  # stat_text 应包含 Stress
  stopifnot(any(grepl("Stress", res$stat_text)))
  cat("  NMDS stress =", res$stress, "\n")
})

chk("PCA（转录组数据）", {
  res <- beta_calc(
    data = expr_mat,
    meta = expr_meta,
    group = "group",
    method = "PCA",
    scale = TRUE,
    verbose = FALSE
  )
  stopifnot(res$method == "PCA")
  stopifnot(is.null(res$stress))
  stopifnot(grepl("PC1", res$axis_labs[1]))
  cat("  axis_labs:", res$axis_labs, "\n")
})

# ── Part 2: 输出结构 ──────────────────────────────────────────────────────────
cat("\n--- Part 2: 输出结构 ---\n")

chk("coords 包含分组列且行名为样本名", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    verbose = FALSE
  )
  stopifnot("Management" %in% names(res$coords))
  stopifnot("PC1" %in% names(res$coords))
  stopifnot("PC2" %in% names(res$coords))
  stopifnot(all(rownames(res$coords) == rownames(dune)))
})

chk("dist_mat 是 dist 对象且维度正确", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    verbose = FALSE
  )
  stopifnot(inherits(res$dist_mat, "dist"))
  stopifnot(attr(res$dist_mat, "Size") == nrow(dune))
})

chk("axis_labs 格式正确（含百分比）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    verbose = FALSE
  )
  stopifnot(length(res$axis_labs) == 2L)
  stopifnot(all(grepl("%", res$axis_labs)))
})

chk("stat_text 不为 NULL（默认 adonis2）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    verbose = FALSE
  )
  stopifnot(!is.null(res$stat_text))
  stopifnot(any(grepl("Adonis", res$stat_text)))
})

# ── Part 3: 统计方法 ──────────────────────────────────────────────────────────
cat("\n--- Part 3: 统计方法 ---\n")

chk("stat='anosim'", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    stat = "anosim",
    verbose = FALSE
  )
  stopifnot(res$stat_table$method == "anosim")
  stopifnot(any(grepl("ANOSIM", res$stat_text)))
})

chk("stat='mrpp'", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    stat = "mrpp",
    verbose = FALSE
  )
  stopifnot(res$stat_table$method == "mrpp")
})

chk("stat=c('adonis2','anosim')（多种方法）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    stat = c("adonis2", "anosim"),
    verbose = FALSE
  )
  stopifnot(nrow(res$stat_table) == 2L)
  stopifnot(length(res$stat_text) == 2L)
})

chk("stat='none'（跳过统计）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    stat = "none",
    verbose = FALSE
  )
  stopifnot(is.null(res$stat_table))
  stopifnot(is.null(res$stat_text))
})

# ── Part 4: 两两比较 ──────────────────────────────────────────────────────────
cat("\n--- Part 4: 两两比较 ---\n")

chk("pairwise=FALSE 时 pair_table 为 NULL", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    pairwise = FALSE,
    verbose = FALSE
  )
  stopifnot(is.null(res$pair_table))
})

chk("pairwise=TRUE 返回正确行数（C(4,2)=6 对）", {
  # Management 有 4 个水平：BF, HF, NM, SF
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    pairwise = TRUE,
    verbose = FALSE
  )
  stopifnot(!is.null(res$pair_table))
  n_pairs <- choose(nlevels(factor(dune.env$Management)), 2)
  stopifnot(nrow(res$pair_table) == n_pairs)
  stopifnot(all(c("group1", "group2", "p_value") %in% names(res$pair_table)))
  cat("  pair_table 行数：", nrow(res$pair_table), "（期望", n_pairs, "）\n")
  print(res$pair_table[, c(
    "group1",
    "group2",
    "method",
    "statistic",
    "p_value"
  )])
})

chk("pairwise + 多种方法（每对产生 2 行）", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    stat = c("adonis2", "anosim"),
    pairwise = TRUE,
    verbose = FALSE
  )
  n_pairs <- choose(nlevels(factor(dune.env$Management)), 2)
  stopifnot(nrow(res$pair_table) == n_pairs * 2L)
})

# ── Part 5: 三种输入方式 ──────────────────────────────────────────────────────
cat("\n--- Part 5: 输入方式 ---\n")

chk("otu 表输入（特征×样本，自动转置）", {
  # 模拟 OTU 表：特征×样本
  otu_tab <- t(dune) # 30 species × 20 sites
  res <- beta_calc(
    otu = otu_tab,
    meta = dune.env,
    group = "Management",
    verbose = FALSE
  )
  stopifnot(nrow(res$coords) == ncol(otu_tab)) # 20 sites
  stopifnot("Management" %in% names(res$coords))
})

chk("data 行名与 meta 行名对齐（非有序）", {
  # 打乱 meta 行顺序
  shuffled_meta <- dune.env[sample(nrow(dune.env)), ]
  res <- beta_calc(
    data = dune,
    meta = shuffled_meta,
    group = "Management",
    verbose = FALSE
  )
  # 结果应基于匹配行名，样本数不变
  stopifnot(nrow(res$coords) == nrow(dune))
})

chk("data 行名与 meta 行名均无（按行序匹配）", {
  d_noname <- dune
  m_noname <- dune.env
  rownames(d_noname) <- NULL
  rownames(m_noname) <- NULL
  res <- beta_calc(
    data = d_noname,
    meta = m_noname,
    group = "Management",
    verbose = FALSE
  )
  stopifnot(nrow(res$coords) == nrow(dune))
})

# ── Part 6: 不同距离类型 ──────────────────────────────────────────────────────
cat("\n--- Part 6: 距离类型 ---\n")

for (d in c("euclidean", "jaccard", "manhattan", "horn")) {
  local({
    dist_name <- d
    chk(paste0("dist='", dist_name, "'"), {
      res <- beta_calc(
        data = dune,
        meta = dune.env,
        group = "Management",
        dist = dist_name,
        verbose = FALSE
      )
      stopifnot(inherits(res$dist_mat, "dist"))
    })
  })
}

# ── Part 7: ndim 参数 ────────────────────────────────────────────────────────
cat("\n--- Part 7: ndim ---\n")

chk("ndim=3 返回 3 列坐标", {
  res <- beta_calc(
    data = dune,
    meta = dune.env,
    group = "Management",
    ndim = 3L,
    verbose = FALSE
  )
  coord_cols <- grep("^PC|^NMDS", names(res$coords), value = TRUE)
  stopifnot(length(coord_cols) == 3L)
  stopifnot(length(res$axis_labs) == 3L)
})

# ── Part 8: 与 pcoa_plot() 对接 ──────────────────────────────────────────────
cat("\n--- Part 8: 与 pcoa_plot() 对接 ---\n")

PCOA_FUNC <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/pcoa_plot/function.R"

if (file.exists(PCOA_FUNC)) {
  source(PCOA_FUNC)

  chk("beta_calc → pcoa_plot（density 边际）", {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
      stop("需要 ggExtra 包")
    }
    res <- beta_calc(
      data = dune,
      meta = dune.env,
      group = "Management",
      verbose = FALSE
    )
    p <- pcoa_plot(
      res$coords,
      x = "PC1",
      y = "PC2",
      group.by = "Management",
      xlab = res$axis_labs[1],
      ylab = res$axis_labs[2],
      stat_text = res$stat_text,
      marginal = "density"
    )
    ggplot2::ggsave(
      file.path(TEST_OUT, "integration_density.png"),
      p,
      width = 6,
      height = 5
    )
  })

  chk("beta_calc → pcoa_plot（boxplot 边际）", {
    if (!requireNamespace("aplot", quietly = TRUE)) {
      stop("需要 aplot 包")
    }
    res <- beta_calc(
      data = dune,
      meta = dune.env,
      group = "Management",
      stat = c("adonis2", "anosim"),
      verbose = FALSE
    )
    p <- pcoa_plot(
      res$coords,
      x = "PC1",
      y = "PC2",
      group.by = "Management",
      xlab = res$axis_labs[1],
      ylab = res$axis_labs[2],
      stat_text = res$stat_text,
      marginal = "boxplot"
    )
    ggplot2::ggsave(
      file.path(TEST_OUT, "integration_boxplot.png"),
      p,
      width = 7,
      height = 5
    )
  })

  chk("PCA → pcoa_plot（转录组场景）", {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
      stop("需要 ggExtra 包")
    }
    res <- beta_calc(
      data = expr_mat,
      meta = expr_meta,
      group = "group",
      method = "PCA",
      scale = TRUE,
      verbose = FALSE
    )
    p <- pcoa_plot(
      res$coords,
      x = "PC1",
      y = "PC2",
      group.by = "group",
      xlab = res$axis_labs[1],
      ylab = res$axis_labs[2],
      stat_text = res$stat_text,
      marginal = "density"
    )
    ggplot2::ggsave(
      file.path(TEST_OUT, "integration_pca_density.png"),
      p,
      width = 6,
      height = 5
    )
  })
} else {
  cat("  ⚠ pcoa_plot/function.R 不存在，跳过对接测试\n")
}

# ── Part 9: 错误处理 ──────────────────────────────────────────────────────────
cat("\n--- Part 9: 错误处理 ---\n")

chk("未提供任何数据输入报错", {
  result <- tryCatch(
    beta_calc(meta = dune.env, group = "Management"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("未指定 group 报错", {
  result <- tryCatch(
    beta_calc(data = dune, meta = dune.env),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("group 列不存在报错", {
  result <- tryCatch(
    beta_calc(data = dune, meta = dune.env, group = "nonexistent"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("无效 stat 选项报错", {
  result <- tryCatch(
    beta_calc(
      data = dune,
      meta = dune.env,
      group = "Management",
      stat = "invalid_stat"
    ),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

chk("data 含 NA（自动移除，不报错）", {
  d_na <- dune
  d_na[1, 1:5] <- NA
  res <- beta_calc(
    data = d_na,
    meta = dune.env,
    group = "Management",
    verbose = FALSE
  )
  # 移除 1 个样本后应有 19 行
  stopifnot(nrow(res$coords) == 19L)
})

chk("非数值 data 报错", {
  d_bad <- dune
  d_bad$V1_str <- as.character(d_bad[, 1])
  d_bad <- d_bad[, c("V1_str", names(dune)[2:5])]
  result <- tryCatch(
    beta_calc(data = d_bad, meta = dune.env, group = "Management"),
    error = function(e) "caught"
  )
  stopifnot(result == "caught")
})

# ── 汇总 ──────────────────────────────────────────────────────────────────────
total <- pass + fail
cat("\n============================================================\n")
cat(sprintf("  结果：%d / %d 通过", pass, total))
if (fail > 0L) {
  cat(sprintf("，%d 个失败", fail))
}
cat("\n  图片已保存至：", TEST_OUT, "\n")
cat("============================================================\n")
