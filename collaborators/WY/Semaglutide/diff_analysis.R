#' ============================================================================
#' diff_abundance: 通用物种丰度差异分析函数 (含效应量)
#' ============================================================================
#'
#' 独立于 microeco 包，直接对物种丰度矩阵进行差异分析。
#' 每种方法均返回其对应的效应量 (effect size) 指标。
#'
#' @param abund_table matrix 或 data.frame，行为物种/特征(features)，列为样本(samples)。
#' @param group 与 abund_table 列（样本）一一对应的分组向量。
#' @param method 差异分析方法，见下表。
#' @param p_adjust_method P 值校正方法，默认 "fdr"。
#' @param alpha 显著性阈值，默认 0.05。
#' @param filter_thres 丰度过滤阈值，默认 0。
#' @param paired 是否配对检验（仅 wilcox / t.test），默认 FALSE。
#' @param lefse_lda_cutoff LEfSe LDA 阈值，默认 2。
#' @param lefse_norm LEfSe 归一化值，默认 1e6。
#' @param boots LEfSe bootstrap 次数，默认 30。
#' @param nresam LEfSe bootstrap 采样比例，默认 2/3。
#' @param beta_pseudo betareg/glmm_beta 中 0/1 替代值，默认 .Machine$double.eps。
#' @param ... 额外参数传递给底层函数。
#'
#' @details
#' 各方法返回的效应量 (EffectSize) 指标:
#' \tabular{lll}{
#'   方法          \tab 效应量指标               \tab 说明 \cr
#'   wilcox        \tab Cliff's Delta            \tab 非参数效应量, 范围 [-1, 1] \cr
#'   t.test        \tab Cohen's d                \tab 标准化均值差 \cr
#'   KW            \tab Epsilon-squared          \tab H / (n^2-1)/(n+1), 范围 [0, 1] \cr
#'   KW_dunn       \tab r (Z / sqrt(N))          \tab 秩效应量, 范围 [0, 1] \cr
#'   anova         \tab Eta-squared              \tab SS_between / SS_total \cr
#'   lm            \tab R-squared                \tab 决定系数 \cr
#'   lefse         \tab LDA score                \tab 线性判别分析效应量 \cr
#'   betareg       \tab Pseudo R-squared         \tab Beta 回归伪决定系数 \cr
#'   DESeq2        \tab log2FoldChange           \tab 对数倍变化 \cr
#'   edgeR         \tab logFC                    \tab 对数倍变化 \cr
#'   ALDEx2_t      \tab ALDEx2 effect            \tab CLR 空间的效应量 \cr
#'   ALDEx2_kw     \tab GLM effect               \tab GLM 效应量 \cr
#'   metagenomeSeq \tab logFC                    \tab ZIG 模型对数倍变化 \cr
#'   ancombc2      \tab lfc (bias-corrected)     \tab 偏差校正后对数倍变化 \cr
#'   linda         \tab log2FoldChange           \tab LinDA 对数倍变化 \cr
#' }
#'
#' @return data.frame，统一包含以下核心列:
#'   - Feature: 物种/特征名
#'   - Comparison: 比较组
#'   - Group: 丰度最高的组
#'   - Method: 统计方法
#'   - EffectSize: 效应量数值 (各方法不同，见上表)
#'   - EffectSize_type: 效应量指标名称
#'   - P.unadj: 原始 P 值
#'   - P.adj: 校正后 P 值
#'   - Significance: 显著性标记 (***, **, *, ns)
#'   - 各方法特有的详细列
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, 10), nrow = 50, ncol = 10,
#'               dimnames = list(paste0("Sp", 1:50), paste0("S", 1:10)))
#' group <- rep(c("Ctrl", "Treat"), each = 5)
#'
#' res <- diff_abundance(mat, group, method = "wilcox")
#' head(res[, c("Feature", "EffectSize", "EffectSize_type", "P.adj")])
#'
#' @export

diff_abundance <- function(
  abund_table,
  group,
  method = c(
    "wilcox",
    "t.test",
    "KW",
    "KW_dunn",
    "anova",
    "lefse",
    "lm",
    "betareg",
    "DESeq2",
    "edgeR",
    "ALDEx2_t",
    "ALDEx2_kw",
    "metagenomeSeq",
    "ancombc2",
    "linda"
  ),
  p_adjust_method = "fdr",
  alpha = 0.05,
  filter_thres = 0,
  paired = FALSE,
  lefse_lda_cutoff = 2,
  lefse_norm = 1e6,
  boots = 30,
  nresam = 2 / 3,
  beta_pseudo = .Machine$double.eps,
  ...
) {
  method <- match.arg(method)
  abund_table <- as.data.frame(abund_table)

  if (ncol(abund_table) != length(group)) {
    stop(
      "abund_table 列数 (",
      ncol(abund_table),
      ") != group 长度 (",
      length(group),
      ")"
    )
  }
  if (is.null(rownames(abund_table))) {
    rownames(abund_table) <- paste0("Feature_", seq_len(nrow(abund_table)))
  }

  group <- as.factor(group)
  if (nlevels(group) < 2) {
    stop("group 中至少需要 2 个不同的组。")
  }

  # 丰度过滤
  if (filter_thres > 0) {
    keep <- rowMeans(abund_table) >= filter_thres
    if (sum(keep) == 0) {
      stop("filter_thres 过高，过滤后无剩余特征。")
    }
    abund_table <- abund_table[keep, , drop = FALSE]
    message("过滤后保留 ", nrow(abund_table), " 个特征")
  }

  # 分发
  res <- switch(
    method,
    "wilcox" = .diff_wilcox(abund_table, group, paired, p_adjust_method),
    "t.test" = .diff_ttest(abund_table, group, paired, p_adjust_method),
    "KW" = .diff_kw(abund_table, group, p_adjust_method),
    "KW_dunn" = .diff_kw_dunn(abund_table, group, p_adjust_method),
    "anova" = .diff_anova(abund_table, group, p_adjust_method),
    "lm" = .diff_lm(abund_table, group, p_adjust_method),
    "lefse" = .diff_lefse(
      abund_table,
      group,
      alpha,
      lefse_lda_cutoff,
      lefse_norm,
      boots,
      nresam,
      p_adjust_method
    ),
    "betareg" = .diff_betareg(abund_table, group, beta_pseudo, p_adjust_method),
    "DESeq2" = .diff_deseq2(abund_table, group, p_adjust_method, ...),
    "edgeR" = .diff_edger(abund_table, group, p_adjust_method),
    "ALDEx2_t" = .diff_aldex2(abund_table, group, "t", p_adjust_method, ...),
    "ALDEx2_kw" = .diff_aldex2(abund_table, group, "kw", p_adjust_method, ...),
    "metagenomeSeq" = .diff_metagenomeseq(
      abund_table,
      group,
      p_adjust_method,
      ...
    ),
    "ancombc2" = .diff_ancombc2(abund_table, group, p_adjust_method, ...),
    "linda" = .diff_linda(abund_table, group, p_adjust_method, ...)
  )

  # 显著性标记
  if ("P.adj" %in% colnames(res)) {
    res$Significance <- ifelse(
      is.na(res$P.adj),
      "ns",
      ifelse(
        res$P.adj < 0.001,
        "***",
        ifelse(res$P.adj < 0.01, "**", ifelse(res$P.adj < 0.05, "*", "ns"))
      )
    )
  }

  # 排序: 按效应量绝对值降序
  if ("EffectSize" %in% colnames(res)) {
    res <- res[order(abs(res$EffectSize), decreasing = TRUE, na.last = TRUE), ]
  }
  rownames(res) <- NULL
  res
}


# ============================================================================
# 内部工具函数
# ============================================================================

.get_max_group <- function(abund_table, group, use_median = TRUE) {
  fn <- if (use_median) median else mean
  sapply(seq_len(nrow(abund_table)), function(i) {
    vals <- split(as.numeric(abund_table[i, ]), group)
    stats <- sapply(vals, fn, na.rm = TRUE)
    names(which.max(stats))
  })
}

.get_pairs <- function(group_names) {
  combn(group_names, 2, simplify = FALSE)
}

# Cliff's Delta: 非参数效应量, 范围 [-1, 1]
# 解释: |d| < 0.147 可忽略, < 0.33 小, < 0.474 中, >= 0.474 大
.cliff_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  if (nx == 0 || ny == 0) {
    return(NA_real_)
  }
  more <- sum(outer(x, y, ">"))
  less <- sum(outer(x, y, "<"))
  (more - less) / (nx * ny)
}

# Cohen's d: 标准化均值差
# 解释: |d| < 0.2 可忽略, 0.2-0.5 小, 0.5-0.8 中, > 0.8 大
.cohen_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  mx <- mean(x, na.rm = TRUE)
  my <- mean(y, na.rm = TRUE)
  sx <- sd(x, na.rm = TRUE)
  sy <- sd(y, na.rm = TRUE)
  sp <- sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  if (is.na(sp) || sp == 0) {
    return(0)
  }
  (mx - my) / sp
}


# ============================================================================
# 1. Wilcoxon + Cliff's Delta
# ============================================================================
.diff_wilcox <- function(abund_table, group, paired, p_adjust_method) {
  pairs <- .get_pairs(levels(group))

  results <- lapply(pairs, function(pr) {
    idx1 <- which(group == pr[1])
    idx2 <- which(group == pr[2])
    if (paired && length(idx1) != length(idx2)) {
      stop("配对检验要求两组样本数相同")
    }

    out <- do.call(
      rbind,
      lapply(seq_len(nrow(abund_table)), function(i) {
        x <- as.numeric(abund_table[i, idx1])
        y <- as.numeric(abund_table[i, idx2])
        pval <- tryCatch(
          wilcox.test(x, y, paired = paired, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
        delta <- .cliff_delta(x, y)
        max_g <- if (median(x, na.rm = TRUE) >= median(y, na.rm = TRUE)) {
          pr[1]
        } else {
          pr[2]
        }
        data.frame(
          pval = pval,
          delta = delta,
          max_g = max_g,
          stringsAsFactors = FALSE
        )
      })
    )

    data.frame(
      Feature = rownames(abund_table),
      Comparison = paste(pr[1], "vs", pr[2]),
      Group = out$max_g,
      Method = "Wilcoxon",
      EffectSize = round(out$delta, 4),
      EffectSize_type = "Cliff's Delta",
      P.unadj = out$pval,
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, results)
  res$P.adj <- p.adjust(res$P.unadj, method = p_adjust_method)
  res
}


# ============================================================================
# 2. t.test + Cohen's d
# ============================================================================
.diff_ttest <- function(abund_table, group, paired, p_adjust_method) {
  pairs <- .get_pairs(levels(group))

  results <- lapply(pairs, function(pr) {
    idx1 <- which(group == pr[1])
    idx2 <- which(group == pr[2])
    if (paired && length(idx1) != length(idx2)) {
      stop("配对检验要求两组样本数相同")
    }

    out <- do.call(
      rbind,
      lapply(seq_len(nrow(abund_table)), function(i) {
        x <- as.numeric(abund_table[i, idx1])
        y <- as.numeric(abund_table[i, idx2])
        pval <- tryCatch(
          t.test(x, y, paired = paired)$p.value,
          error = function(e) NA_real_
        )
        d <- .cohen_d(x, y)
        max_g <- if (mean(x, na.rm = TRUE) >= mean(y, na.rm = TRUE)) {
          pr[1]
        } else {
          pr[2]
        }
        data.frame(pval = pval, d = d, max_g = max_g, stringsAsFactors = FALSE)
      })
    )

    data.frame(
      Feature = rownames(abund_table),
      Comparison = paste(pr[1], "vs", pr[2]),
      Group = out$max_g,
      Method = "t.test",
      EffectSize = round(out$d, 4),
      EffectSize_type = "Cohen's d",
      P.unadj = out$pval,
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, results)
  res$P.adj <- p.adjust(res$P.unadj, method = p_adjust_method)
  res
}


# ============================================================================
# 3. Kruskal-Wallis + Epsilon-squared
# ============================================================================
.diff_kw <- function(abund_table, group, p_adjust_method) {
  n <- length(group)
  max_grp <- .get_max_group(abund_table, group, use_median = TRUE)

  out <- do.call(
    rbind,
    lapply(seq_len(nrow(abund_table)), function(i) {
      vals <- as.numeric(abund_table[i, ])
      tryCatch(
        {
          kt <- kruskal.test(vals ~ group)
          H <- as.numeric(kt$statistic)
          # Epsilon-squared = H / ((n^2 - 1) / (n + 1))
          eps2 <- H / ((n^2 - 1) / (n + 1))
          data.frame(H = H, pval = kt$p.value, eps2 = eps2)
        },
        error = function(e) {
          data.frame(H = NA_real_, pval = NA_real_, eps2 = NA_real_)
        }
      )
    })
  )

  data.frame(
    Feature = rownames(abund_table),
    Comparison = paste(levels(group), collapse = " | "),
    Group = max_grp,
    Method = "Kruskal-Wallis",
    EffectSize = round(out$eps2, 4),
    EffectSize_type = "Epsilon-squared",
    H_statistic = round(out$H, 4),
    P.unadj = out$pval,
    P.adj = p.adjust(out$pval, method = p_adjust_method),
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 4. Dunn's Test + rank effect size r = |Z| / sqrt(N)
# ============================================================================
.diff_kw_dunn <- function(abund_table, group, p_adjust_method) {
  if (!requireNamespace("FSA", quietly = TRUE)) {
    stop("需要安装 FSA 包: install.packages('FSA')")
  }

  n_total <- length(group)

  results <- lapply(seq_len(nrow(abund_table)), function(i) {
    vals <- as.numeric(abund_table[i, ])
    tryCatch(
      {
        dt <- FSA::dunnTest(vals ~ group, method = p_adjust_method)
        rr <- dt$res
        rr$r <- abs(rr$Z) / sqrt(n_total)
        grp_max <- sapply(seq_len(nrow(rr)), function(j) {
          pair <- trimws(unlist(strsplit(rr$Comparison[j], " - ")))
          m1 <- median(vals[group == pair[1]], na.rm = TRUE)
          m2 <- median(vals[group == pair[2]], na.rm = TRUE)
          if (m1 >= m2) pair[1] else pair[2]
        })
        data.frame(
          Feature = rownames(abund_table)[i],
          Comparison = rr$Comparison,
          Group = grp_max,
          Method = "Dunn",
          EffectSize = round(rr$r, 4),
          EffectSize_type = "r (|Z|/sqrt(N))",
          Z_statistic = round(rr$Z, 4),
          P.unadj = rr$P.unadj,
          P.adj = rr$P.adj,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) NULL
    )
  })
  do.call(rbind, results)
}


# ============================================================================
# 5. ANOVA + Eta-squared
# ============================================================================
.diff_anova <- function(abund_table, group, p_adjust_method) {
  max_grp <- .get_max_group(abund_table, group, use_median = FALSE)

  out <- do.call(
    rbind,
    lapply(seq_len(nrow(abund_table)), function(i) {
      vals <- as.numeric(abund_table[i, ])
      tryCatch(
        {
          fit <- aov(vals ~ group)
          ss <- summary(fit)[[1]]
          ss_between <- ss[["Sum Sq"]][1]
          ss_total <- sum(ss[["Sum Sq"]])
          eta2 <- ss_between / ss_total
          data.frame(
            Fval = ss[["F value"]][1],
            pval = ss[["Pr(>F)"]][1],
            eta2 = eta2
          )
        },
        error = function(e) {
          data.frame(Fval = NA_real_, pval = NA_real_, eta2 = NA_real_)
        }
      )
    })
  )

  data.frame(
    Feature = rownames(abund_table),
    Comparison = paste(levels(group), collapse = " | "),
    Group = max_grp,
    Method = "ANOVA",
    EffectSize = round(out$eta2, 4),
    EffectSize_type = "Eta-squared",
    F_statistic = round(out$Fval, 4),
    P.unadj = out$pval,
    P.adj = p.adjust(out$pval, method = p_adjust_method),
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 6. 线性模型 + R-squared
# ============================================================================
.diff_lm <- function(abund_table, group, p_adjust_method) {
  max_grp <- .get_max_group(abund_table, group, use_median = FALSE)

  out <- do.call(
    rbind,
    lapply(seq_len(nrow(abund_table)), function(i) {
      vals <- as.numeric(abund_table[i, ])
      tryCatch(
        {
          fit <- lm(vals ~ group)
          sf <- summary(fit)
          fstat <- sf$fstatistic
          pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
          data.frame(
            R2 = sf$r.squared,
            adj_R2 = sf$adj.r.squared,
            pval = as.numeric(pval)
          )
        },
        error = function(e) {
          data.frame(R2 = NA_real_, adj_R2 = NA_real_, pval = NA_real_)
        }
      )
    })
  )

  data.frame(
    Feature = rownames(abund_table),
    Comparison = paste(levels(group), collapse = " | "),
    Group = max_grp,
    Method = "lm",
    EffectSize = round(out$R2, 4),
    EffectSize_type = "R-squared",
    Adj_R2 = round(out$adj_R2, 4),
    P.unadj = out$pval,
    P.adj = p.adjust(out$pval, method = p_adjust_method),
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 7. LEfSe + LDA score
# ============================================================================
.diff_lefse <- function(
  abund_table,
  group,
  alpha,
  lda_cutoff,
  lefse_norm,
  boots,
  nresam,
  p_adjust_method
) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("需要安装 MASS 包: install.packages('MASS')")
  }

  group_names <- levels(group)
  tab <- abund_table

  # 归一化
  if (lefse_norm > 0) {
    cs <- colSums(tab)
    cs[cs == 0] <- 1
    tab <- sweep(tab, 2, cs, "/") * lefse_norm
  }

  # Step 1: KW 筛选
  kw_pvals <- sapply(seq_len(nrow(tab)), function(i) {
    tryCatch(
      kruskal.test(as.numeric(tab[i, ]) ~ group)$p.value,
      error = function(e) NA_real_
    )
  })
  kw_padj <- p.adjust(kw_pvals, method = p_adjust_method)

  # Step 2: 配对 Wilcoxon 一致性检验
  pairs <- .get_pairs(group_names)
  wilcox_pass <- rep(TRUE, nrow(tab))
  for (pr in pairs) {
    i1 <- which(group == pr[1])
    i2 <- which(group == pr[2])
    wp <- sapply(seq_len(nrow(tab)), function(i) {
      tryCatch(
        wilcox.test(
          as.numeric(tab[i, i1]),
          as.numeric(tab[i, i2]),
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )
    })
    wilcox_pass <- wilcox_pass & (wp < alpha | is.na(wp))
  }

  sig_idx <- which(kw_padj < alpha & wilcox_pass)
  if (length(sig_idx) == 0) {
    message("LEfSe: 无显著差异特征")
    return(data.frame(
      Feature = character(),
      Comparison = character(),
      Group = character(),
      Method = character(),
      EffectSize = numeric(),
      EffectSize_type = character(),
      P.unadj = numeric(),
      P.adj = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  # Step 3: LDA 效应量
  lda_scores <- sapply(sig_idx, function(i) {
    fv <- as.numeric(tab[i, ])
    bs <- replicate(boots, {
      si <- unlist(lapply(group_names, function(g) {
        gi <- which(group == g)
        sample(gi, max(1, round(length(gi) * nresam)), replace = TRUE)
      }))
      tryCatch(
        {
          df_b <- data.frame(val = fv[si], grp = group[si])
          fit <- MASS::lda(grp ~ val, data = df_b)
          pred <- predict(fit, df_b)
          means <- tapply(pred$x[, 1], df_b$grp, mean)
          abs(max(means) - min(means))
        },
        error = function(e) 0
      )
    })
    log10(1 + median(bs, na.rm = TRUE))
  })

  max_grp <- .get_max_group(
    tab[sig_idx, , drop = FALSE],
    group,
    use_median = TRUE
  )

  res <- data.frame(
    Feature = rownames(tab)[sig_idx],
    Comparison = paste(group_names, collapse = " | "),
    Group = max_grp,
    Method = "LEfSe",
    EffectSize = round(lda_scores, 4),
    EffectSize_type = "LDA score",
    P.unadj = kw_pvals[sig_idx],
    P.adj = kw_padj[sig_idx],
    stringsAsFactors = FALSE
  )
  res[res$EffectSize >= lda_cutoff, ]
}


# ============================================================================
# 8. Beta 回归 + Pseudo R-squared
# ============================================================================
.diff_betareg <- function(abund_table, group, beta_pseudo, p_adjust_method) {
  if (!requireNamespace("betareg", quietly = TRUE)) {
    stop("需要安装 betareg 包: install.packages('betareg')")
  }

  max_grp <- .get_max_group(abund_table, group, use_median = FALSE)

  res_list <- lapply(seq_len(nrow(abund_table)), function(i) {
    vals <- as.numeric(abund_table[i, ])
    vals[vals <= 0] <- beta_pseudo
    vals[vals >= 1] <- 1 / (1 + beta_pseudo)
    tryCatch(
      {
        df <- data.frame(y = vals, grp = group)
        fit <- betareg::betareg(y ~ grp, data = df)
        sf <- summary(fit)
        pseudo_r2 <- sf$pseudo.r.squared
        coefs <- sf$coefficients$mean
        pval <- if (nrow(coefs) > 1) {
          min(coefs[-1, "Pr(>|z|)"], na.rm = TRUE)
        } else {
          NA_real_
        }
        est <- if (nrow(coefs) > 1) coefs[2, "Estimate"] else NA_real_
        se <- if (nrow(coefs) > 1) coefs[2, "Std. Error"] else NA_real_
        data.frame(
          pseudo_r2 = pseudo_r2,
          Estimate = est,
          Std.Error = se,
          P.unadj = pval
        )
      },
      error = function(e) {
        data.frame(
          pseudo_r2 = NA_real_,
          Estimate = NA_real_,
          Std.Error = NA_real_,
          P.unadj = NA_real_
        )
      }
    )
  })

  res <- do.call(rbind, res_list)
  data.frame(
    Feature = rownames(abund_table),
    Comparison = paste(levels(group), collapse = " | "),
    Group = max_grp,
    Method = "betareg",
    EffectSize = round(res$pseudo_r2, 4),
    EffectSize_type = "Pseudo R-squared",
    Estimate = round(res$Estimate, 4),
    Std.Error = round(res$Std.Error, 4),
    P.unadj = res$P.unadj,
    P.adj = p.adjust(res$P.unadj, method = p_adjust_method),
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 9. DESeq2 + log2FoldChange
# ============================================================================
.diff_deseq2 <- function(abund_table, group, p_adjust_method, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("需要安装 DESeq2: BiocManager::install('DESeq2')")
  }

  ct <- round(as.matrix(abund_table))
  ct[ct < 0] <- 0
  pairs <- .get_pairs(levels(group))

  results <- lapply(pairs, function(pr) {
    idx <- which(group %in% pr)
    sc <- ct[, idx, drop = FALSE]
    sg <- droplevels(group[idx])
    keep <- rowSums(sc) > 0
    sc <- sc[keep, , drop = FALSE]
    tryCatch(
      {
        cd <- data.frame(group = sg, row.names = colnames(sc))
        dds <- DESeq2::DESeqDataSetFromMatrix(sc, cd, ~group)
        dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)
        dr <- as.data.frame(DESeq2::results(
          dds,
          pAdjustMethod = p_adjust_method
        ))
        data.frame(
          Feature = rownames(dr),
          Comparison = paste(pr[1], "vs", pr[2]),
          Group = ifelse(dr$log2FoldChange > 0, pr[2], pr[1]),
          Method = "DESeq2",
          EffectSize = round(dr$log2FoldChange, 4),
          EffectSize_type = "log2FoldChange",
          baseMean = round(dr$baseMean, 2),
          lfcSE = round(dr$lfcSE, 4),
          stat = round(dr$stat, 4),
          P.unadj = dr$pvalue,
          P.adj = dr$padj,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        warning("DESeq2 ", pr[1], " vs ", pr[2], " 失败: ", e$message)
        NULL
      }
    )
  })
  do.call(rbind, results)
}


# ============================================================================
# 10. edgeR + logFC
# ============================================================================
.diff_edger <- function(abund_table, group, p_adjust_method) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("需要安装 edgeR: BiocManager::install('edgeR')")
  }

  ct <- round(as.matrix(abund_table))
  ct[ct < 0] <- 0
  pairs <- .get_pairs(levels(group))

  results <- lapply(pairs, function(pr) {
    idx <- which(group %in% pr)
    sc <- ct[, idx, drop = FALSE]
    sg <- droplevels(group[idx])
    keep <- rowSums(sc) > 0
    sc <- sc[keep, , drop = FALSE]
    tryCatch(
      {
        dge <- edgeR::DGEList(counts = sc, group = sg)
        dge <- edgeR::calcNormFactors(dge)
        dge <- edgeR::estimateDisp(dge)
        et <- edgeR::exactTest(dge)
        tt <- edgeR::topTags(et, n = Inf, adjust.method = p_adjust_method)$table
        data.frame(
          Feature = rownames(tt),
          Comparison = paste(pr[1], "vs", pr[2]),
          Group = ifelse(tt$logFC > 0, pr[2], pr[1]),
          Method = "edgeR",
          EffectSize = round(tt$logFC, 4),
          EffectSize_type = "logFC",
          logCPM = round(tt$logCPM, 4),
          P.unadj = tt$PValue,
          P.adj = tt$FDR,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        warning("edgeR ", pr[1], " vs ", pr[2], " 失败: ", e$message)
        NULL
      }
    )
  })
  do.call(rbind, results)
}


# ============================================================================
# 11. ALDEx2 + effect / glm effect
# ============================================================================
.diff_aldex2 <- function(abund_table, group, test, p_adjust_method, ...) {
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("需要安装 ALDEx2: BiocManager::install('ALDEx2')")
  }

  ct <- round(as.matrix(abund_table))
  ct[ct < 0] <- 0
  keep <- rowSums(ct) > 0
  ct <- ct[keep, , drop = FALSE]

  tryCatch(
    {
      if (test == "t") {
        ar <- ALDEx2::aldex(
          ct,
          conditions = as.character(group),
          test = "t",
          ...
        )
        data.frame(
          Feature = rownames(ar),
          Comparison = paste(levels(group), collapse = " vs "),
          Group = ifelse(ar$diff.btw > 0, levels(group)[2], levels(group)[1]),
          Method = "ALDEx2_t",
          EffectSize = round(ar$effect, 4),
          EffectSize_type = "ALDEx2 effect",
          diff.btw = round(ar$diff.btw, 4),
          diff.win = round(ar$diff.win, 4),
          P.unadj = ar$wi.ep,
          P.adj = ar$wi.eBH,
          stringsAsFactors = FALSE
        )
      } else {
        ar <- ALDEx2::aldex(
          ct,
          conditions = as.character(group),
          test = "kw",
          ...
        )
        data.frame(
          Feature = rownames(ar),
          Comparison = paste(levels(group), collapse = " | "),
          Group = NA_character_,
          Method = "ALDEx2_kw",
          EffectSize = round(ar$glm.effect, 4),
          EffectSize_type = "GLM effect",
          P.unadj = ar$kw.ep,
          P.adj = ar$kw.eBH,
          stringsAsFactors = FALSE
        )
      }
    },
    error = function(e) stop("ALDEx2 运行失败: ", e$message)
  )
}


# ============================================================================
# 12. metagenomeSeq + logFC (ZIG)
# ============================================================================
.diff_metagenomeseq <- function(abund_table, group, p_adjust_method, ...) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("需要安装 metagenomeSeq: BiocManager::install('metagenomeSeq')")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("需要安装 Biobase: BiocManager::install('Biobase')")
  }

  ct <- round(as.matrix(abund_table))
  ct[ct < 0] <- 0
  pairs <- .get_pairs(levels(group))

  results <- lapply(pairs, function(pr) {
    idx <- which(group %in% pr)
    sc <- ct[, idx, drop = FALSE]
    sg <- droplevels(group[idx])
    tryCatch(
      {
        pheno <- Biobase::AnnotatedDataFrame(data.frame(
          group = sg,
          row.names = colnames(sc)
        ))
        mgs <- metagenomeSeq::newMRexperiment(sc, phenoData = pheno)
        p <- metagenomeSeq::cumNormStat(mgs, pFlag = FALSE)
        mgs <- metagenomeSeq::cumNorm(mgs, p = p)
        mod <- model.matrix(~group, data = Biobase::pData(mgs))
        fit <- metagenomeSeq::fitZig(mgs, mod, verbose = FALSE)
        mr <- metagenomeSeq::MRcoefs(
          fit,
          number = nrow(sc),
          adjustMethod = p_adjust_method,
          ...
        )
        data.frame(
          Feature = rownames(mr),
          Comparison = paste(pr[1], "vs", pr[2]),
          Group = ifelse(mr[, 1] > 0, pr[2], pr[1]),
          Method = "metagenomeSeq",
          EffectSize = round(mr[, 1], 4),
          EffectSize_type = "logFC (ZIG)",
          P.unadj = mr$pvalues,
          P.adj = mr$adjPvalues,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        warning("metagenomeSeq ", pr[1], " vs ", pr[2], " 失败: ", e$message)
        NULL
      }
    )
  })
  do.call(rbind, results)
}


# ============================================================================
# 13. ANCOM-BC2 + lfc (bias-corrected)
# ============================================================================
.diff_ancombc2 <- function(abund_table, group, p_adjust_method, ...) {
  if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
    stop("需要安装 ANCOMBC: BiocManager::install('ANCOMBC')")
  }
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    stop(
      "需要安装 TreeSummarizedExperiment: BiocManager::install('TreeSummarizedExperiment')"
    )
  }

  ct <- round(as.matrix(abund_table))
  ct[ct < 0] <- 0
  sd <- data.frame(group = group, row.names = colnames(ct))

  tryCatch(
    {
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = ct),
        colData = sd
      )
      ar <- ANCOMBC::ancombc2(
        data = tse,
        fix_formula = "group",
        p_adj_method = p_adjust_method,
        group = "group",
        global = FALSE,
        ...
      )
      rd <- ar$res
      lfc_cols <- grep("^lfc_group", colnames(rd), value = TRUE)
      p_cols <- grep("^p_group", colnames(rd), value = TRUE)
      q_cols <- grep("^q_group", colnames(rd), value = TRUE)
      se_cols <- grep("^se_group", colnames(rd), value = TRUE)

      results <- lapply(seq_along(lfc_cols), function(j) {
        data.frame(
          Feature = rd$taxon,
          Comparison = gsub("^lfc_", "", lfc_cols[j]),
          Method = "ANCOM-BC2",
          EffectSize = round(rd[[lfc_cols[j]]], 4),
          EffectSize_type = "lfc (bias-corrected)",
          SE = if (length(se_cols) >= j) {
            round(rd[[se_cols[j]]], 4)
          } else {
            NA_real_
          },
          P.unadj = rd[[p_cols[j]]],
          P.adj = if (length(q_cols) >= j) {
            rd[[q_cols[j]]]
          } else {
            p.adjust(rd[[p_cols[j]]], method = p_adjust_method)
          },
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, results)
    },
    error = function(e) stop("ANCOM-BC2 运行失败: ", e$message)
  )
}


# ============================================================================
# 14. LinDA + log2FoldChange
# ============================================================================
.diff_linda <- function(abund_table, group, p_adjust_method, ...) {
  if (!requireNamespace("MicrobiomeStat", quietly = TRUE)) {
    stop("需要安装 MicrobiomeStat: install.packages('MicrobiomeStat')")
  }

  ct <- as.data.frame(round(as.matrix(abund_table)))
  ct[ct < 0] <- 0
  meta <- data.frame(group = group, row.names = colnames(ct))

  tryCatch(
    {
      lr <- MicrobiomeStat::linda(
        feature.dat = ct,
        meta.dat = meta,
        formula = "~ group",
        feature.dat.type = "count",
        p.adj.method = p_adjust_method,
        ...
      )
      results <- lapply(names(lr$output), function(nm) {
        rr <- lr$output[[nm]]
        data.frame(
          Feature = rownames(rr),
          Comparison = nm,
          Method = "LinDA",
          EffectSize = round(rr$log2FoldChange, 4),
          EffectSize_type = "log2FoldChange",
          stat = round(rr$stat, 4),
          P.unadj = rr$pvalue,
          P.adj = rr$padj,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, results)
    },
    error = function(e) stop("LinDA 运行失败: ", e$message)
  )
}


# ============================================================================
# 辅助: 丰度汇总
# ============================================================================
#' 计算各组丰度统计量
#' @export
abund_summary <- function(abund_table, group) {
  group <- as.factor(group)
  do.call(
    rbind,
    lapply(levels(group), function(g) {
      sub <- abund_table[, group == g, drop = FALSE]
      data.frame(
        Feature = rownames(abund_table),
        Group = g,
        Mean = round(rowMeans(sub, na.rm = TRUE), 6),
        SD = round(apply(sub, 1, sd, na.rm = TRUE), 6),
        Median = round(apply(sub, 1, median, na.rm = TRUE), 6),
        SE = round(
          apply(sub, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
          6
        ),
        N = ncol(sub),
        stringsAsFactors = FALSE
      )
    })
  )
}


# ============================================================================
# 可视化: 火山图 (适用于含 log2FC 类效应量的方法)
# ============================================================================
#' @export
plot_volcano <- function(
  res,
  log2fc_cutoff = 1,
  pvalue_cutoff = 0.05,
  label_top_n = 10,
  colors = c("#e74c3c", "#3498db", "gray80")
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要 ggplot2")
  }
  if (!"EffectSize" %in% colnames(res)) {
    stop("结果中无 EffectSize 列")
  }

  res$neglog10p <- -log10(res$P.adj + 1e-300)
  res$Type <- ifelse(
    is.na(res$P.adj) | res$P.adj >= pvalue_cutoff,
    "ns",
    ifelse(
      res$EffectSize > log2fc_cutoff,
      "Up",
      ifelse(res$EffectSize < -log2fc_cutoff, "Down", "ns")
    )
  )
  res$Type <- factor(res$Type, levels = c("Up", "Down", "ns"))

  es_type <- res$EffectSize_type[1]

  p <- ggplot2::ggplot(
    res,
    ggplot2::aes(x = EffectSize, y = neglog10p, color = Type)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::scale_color_manual(
      values = stats::setNames(colors, c("Up", "Down", "ns"))
    ) +
    ggplot2::geom_vline(
      xintercept = c(-log2fc_cutoff, log2fc_cutoff),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(pvalue_cutoff),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::labs(x = es_type, y = "-log10(Adjusted P)", color = "Type") +
    ggplot2::theme_bw()

  if (label_top_n > 0) {
    sig <- res[res$Type != "ns" & !is.na(res$Type), ]
    top <- utils::head(sig[order(sig$P.adj), ], label_top_n)
    if (nrow(top) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = top,
          ggplot2::aes(label = Feature),
          size = 3,
          max.overlaps = 20,
          show.legend = FALSE
        )
    }
  }
  p
}


# ============================================================================
# 可视化: 效应量柱状图 (通用，适配所有方法)
# ============================================================================
#' @export
plot_effect_bar <- function(res, top_n = 20, colors = NULL, sig_only = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要 ggplot2")
  }
  if (!"EffectSize" %in% colnames(res)) {
    stop("结果中无 EffectSize 列")
  }

  df <- res[!is.na(res$EffectSize), ]
  if (sig_only && "Significance" %in% colnames(df)) {
    df <- df[df$Significance != "ns", ]
  }
  df <- df[order(abs(df$EffectSize), decreasing = TRUE), ]
  df <- utils::head(df, top_n)
  df$Feature <- factor(df$Feature, levels = rev(df$Feature))

  es_label <- if (!is.null(df$EffectSize_type[1])) {
    df$EffectSize_type[1]
  } else {
    "Effect Size"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Feature, y = EffectSize)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = es_label, x = NULL)

  if ("Group" %in% colnames(df) && !all(is.na(df$Group))) {
    p <- p + ggplot2::geom_bar(ggplot2::aes(fill = Group), stat = "identity")
    if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
  } else {
    p <- p + ggplot2::geom_bar(stat = "identity", fill = "steelblue")
  }

  if ("Significance" %in% colnames(df)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          label = Significance,
          hjust = ifelse(EffectSize >= 0, -0.2, 1.2)
        ),
        size = 3.5
      )
  }
  p
}
