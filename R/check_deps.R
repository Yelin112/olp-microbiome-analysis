# ============================================================================
# check_deps.R — 检查当前机器上函数库所需的 R 包
# ============================================================================
#
# 用法（每台机器换环境或首次 clone 后跑一次）:
#   Rscript R/check_deps.R
#   # 或在 R 中:
#   source("R/check_deps.R"); check_deps()
#   check_deps(only = "utils")        # 只检查某一组
#
# 只报告缺失情况，不会自动安装。
# 新增函数依赖新包时，请同步更新下面的 DEPS 表。

DEPS <- list(
  # 通用绘图基础设施（R/utils/）
  utils = list(
    required = c("ggplot2"),
    optional = c(
      RColorBrewer = "get_colors() 使用 Brewer 色板",
      paletteer = "get_colors() 使用 paletteer 色板",
      officer = "PPTX 导出",
      rvg = "PPTX 导出"
    )
  ),
  compare_plot = list(
    required = c("ggplot2", "ggpubr", "dplyr")
  ),
  pcoa_plot = list(
    required = c("ggplot2"),
    optional = c(
      ggExtra = "marginal = 'density'",
      aplot = "marginal = 'boxplot'",
      RColorBrewer = "Brewer 色板"
    )
  ),
  beta_calc = list(
    required = c("vegan"),
    optional = c(phyloseq = "phyloseq 输入 / UniFrac 距离")
  ),
  abundance_processing = list(
    required = c("vegan"),
    optional = c(
      SRS = "norm: SRS",
      metagenomeSeq = "norm: CSS",
      edgeR = "norm: TMM/RLE",
      DESeq2 = "norm: DESeq2",
      Wrench = "norm: Wrench",
      MASS = "diff: 部分检验",
      FSA = "diff: KW + Dunn",
      betareg = "diff: beta 回归",
      ANCOMBC = "diff: ANCOM-BC",
      ALDEx2 = "diff: ALDEx2",
      MicrobiomeStat = "diff: LinDA 等",
      TreeSummarizedExperiment = "diff: ANCOM-BC 输入"
    )
  ),
  network_analysis = list(
    required = c("igraph", "R6"),
    optional = c(
      ggraph = "网络绘图",
      ggrepel = "标签避让",
      visNetwork = "交互式网络",
      rgexf = "导出 Gephi",
      tidyr = "plot.NetworkComparison(type = 'roles')"
    )
  )
)

check_deps <- function(only = NULL) {
  groups <- if (is.null(only)) names(DEPS) else only
  bad <- setdiff(groups, names(DEPS))
  if (length(bad) > 0) {
    stop("未知分组: ", paste(bad, collapse = ", "), "；可选: ", paste(names(DEPS), collapse = ", "))
  }

  has <- function(p) requireNamespace(p, quietly = TRUE)
  n_missing_required <- 0

  cat(sprintf("R %s | %s\n\n", getRversion(), R.version$platform))
  for (g in groups) {
    d <- DEPS[[g]]
    req <- d$required %||% character(0)
    opt <- d$optional %||% character(0)

    req_miss <- req[!vapply(req, has, logical(1))]
    opt_miss <- opt[!vapply(names(opt), has, logical(1))]
    n_missing_required <- n_missing_required + length(req_miss)

    status <- if (length(req_miss) > 0) "缺必需包" else if (length(opt_miss) > 0) "可用（部分功能受限）" else "OK"
    cat(sprintf("[%s] %s\n", g, status))
    if (length(req_miss) > 0) {
      cat("  必需缺失: ", paste(req_miss, collapse = ", "), "\n", sep = "")
    }
    if (length(opt_miss) > 0) {
      cat(sprintf("  可选缺失: %s → %s\n", names(opt_miss), opt_miss), sep = "")
    }
  }

  cat("\n")
  if (n_missing_required == 0) {
    cat("所有必需包均已安装。\n")
  } else {
    cat(n_missing_required, " 个必需包缺失。\n", sep = "")
  }
  invisible(n_missing_required == 0)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0 && !interactive()) {
  ok <- check_deps()
  quit(status = if (isTRUE(ok)) 0 else 1)
}
