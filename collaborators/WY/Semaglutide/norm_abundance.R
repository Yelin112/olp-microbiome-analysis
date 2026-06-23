#' ============================================================================
#' norm_abundance: 通用物种丰度标准化/转换函数
#' ============================================================================
#'
#' 独立于 microeco 包，直接对丰度矩阵进行标准化或转换。
#' 输入: 行为物种(features)，列为样本(samples) 的矩阵或 data.frame。
#'
#' @param abund_table matrix 或 data.frame，行=物种/特征, 列=样本。
#' @param method 标准化/转换方法，见 Details。
#' @param sample.size 抽平深度 (rarefy/SRS)，默认 NULL 使用最小文库。
#' @param rngseed 随机种子，默认 123。
#' @param replace 抽平是否放回抽样，默认 TRUE。
#' @param pseudocount CLR 中零值伪计数，默认 1。
#' @param intersect.no GMPR 最少共有特征数，默认 10。
#' @param ct.min GMPR 计算比值的最小计数，默认 1。
#' @param condition Wrench 方法的分组条件向量。
#' @param MARGIN vegan 方法的操作维度: 1=样本(行), 2=特征(列)。NULL 使用各方法默认值。
#' @param logbase 对数底数 (log 方法)，默认 2。
#' @param CSS_p CSS 的分位数，默认 NULL 自动计算。
#' @param ... 额外参数传递给底层函数。
#'
#' @details
#' 支持的标准化方法 (Normalization):
#' \describe{
#'   \item{rarefy}{经典抽平法, 基于 sample() 函数}
#'   \item{SRS}{Scaled Ranked Subsampling (Beule & Karlovsky 2020), 需要 SRS 包}
#'   \item{TSS}{Total Sum Scaling, 每个样本除以总和 (即相对丰度)}
#'   \item{CSS}{Cumulative Sum Scaling (metagenomeSeq), 需要 metagenomeSeq 包}
#'   \item{TMM}{Trimmed Mean of M-values (edgeR), 需要 edgeR 包}
#'   \item{RLE}{Relative Log Expression (edgeR), 需要 edgeR 包}
#'   \item{DESeq2}{DESeq2 Median-of-ratios 标准化, 需要 DESeq2 包}
#'   \item{GMPR}{Geometric Mean of Pairwise Ratios}
#'   \item{clr}{Centered Log-Ratio 变换}
#'   \item{rclr}{Robust CLR (忽略零值)}
#'   \item{Wrench}{Group-wise compositional bias (Wrench 包)}
#'   \item{eBay}{Empirical Bayes 标准化}
#' }
#'
#' 支持的转换方法 (Transformation, 基于 vegan::decostand 或自实现):
#' \describe{
#'   \item{total}{除以行/列总和 (默认 MARGIN=1, 即样本)}
#'   \item{max}{除以行/列最大值 (默认 MARGIN=2, 即特征)}
#'   \item{frequency}{除以总和再乘以非零项数 (默认 MARGIN=2)}
#'   \item{normalize}{使平方和为1 (默认 MARGIN=1)}
#'   \item{range}{标准化到 0-1 (默认 MARGIN=2)}
#'   \item{standardize}{零均值单位方差 (默认 MARGIN=2)}
#'   \item{pa}{转为存在/缺失 (0/1)}
#'   \item{chi.square}{卡方变换}
#'   \item{hellinger}{Hellinger 变换 (total 的平方根)}
#'   \item{log}{对数变换}
#'   \item{AST}{Arc sine square root 变换}
#' }
#'
#' @return 与输入同维度的 data.frame，行=物种, 列=样本。
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(200, 50), nrow = 20, ncol = 10,
#'               dimnames = list(paste0("Sp", 1:20), paste0("S", 1:10)))
#'
#' # 相对丰度
#' rel <- norm_abundance(mat, method = "TSS")
#'
#' # CLR 变换
#' clr <- norm_abundance(mat, method = "clr")
#'
#' # 抽平
#' rare <- norm_abundance(mat, method = "rarefy", sample.size = 200)
#'
#' @export

norm_abundance <- function(
    abund_table,
    method = c("rarefy", "SRS", "TSS", "CSS", "TMM", "RLE", "DESeq2",
               "GMPR", "clr", "rclr", "Wrench", "eBay",
               "total", "max", "frequency", "normalize", "range",
               "standardize", "pa", "chi.square", "hellinger", "log", "AST"),
    sample.size = NULL,
    rngseed = 123,
    replace = TRUE,
    pseudocount = 1,
    intersect.no = 10,
    ct.min = 1,
    condition = NULL,
    MARGIN = NULL,
    logbase = 2,
    CSS_p = NULL,
    ...
) {

  method <- match.arg(method)
  abund_table <- as.data.frame(abund_table)

  if (is.null(rownames(abund_table)))
    rownames(abund_table) <- paste0("Feature_", seq_len(nrow(abund_table)))
  if (is.null(colnames(abund_table)))
    colnames(abund_table) <- paste0("Sample_", seq_len(ncol(abund_table)))

  # 内部用转置表 (行=样本, 列=特征) 仅供 vegan 系方法
  t_table <- t(abund_table)

  # %||% 运算符: 若 x 为 NULL 则返回 y
  `%||%` <- function(x, y) if (is.null(x)) y else x

  res <- switch(method,
    # --- 标准化方法 ---
    "rarefy"      = .norm_rarefy(abund_table, sample.size, rngseed, replace),
    "SRS"         = .norm_srs(abund_table, sample.size, rngseed),
    "TSS"         = .norm_tss(abund_table),
    "CSS"         = .norm_css(abund_table, CSS_p, ...),
    "TMM"         = .norm_tmm_rle(abund_table, "TMM", ...),
    "RLE"         = .norm_tmm_rle(abund_table, "RLE", ...),
    "DESeq2"      = .norm_deseq2(abund_table, ...),
    "GMPR"        = .norm_gmpr(abund_table, intersect.no, ct.min),
    "clr"         = .norm_clr(abund_table, pseudocount),
    "rclr"        = .norm_rclr(abund_table),
    "Wrench"      = .norm_wrench(abund_table, condition, ...),
    "eBay"        = .norm_ebay(abund_table),
    # --- vegan decostand 方法 ---
    "total"       = .norm_decostand(t_table, "total",       MARGIN %||% 1, ...),
    "max"         = .norm_decostand(t_table, "max",         MARGIN %||% 2, ...),
    "frequency"   = .norm_decostand(t_table, "frequency",   MARGIN %||% 2, ...),
    "normalize"   = .norm_decostand(t_table, "normalize",   MARGIN %||% 1, ...),
    "range"       = .norm_decostand(t_table, "range",       MARGIN %||% 2, ...),
    "standardize" = .norm_decostand(t_table, "standardize", MARGIN %||% 2, ...),
    "pa"          = .norm_decostand(t_table, "pa",          MARGIN %||% 1, ...),
    "chi.square"  = .norm_decostand(t_table, "chi.square",  MARGIN %||% 1, ...),
    "hellinger"   = .norm_decostand(t_table, "hellinger",   MARGIN %||% 1, ...),
    "log"         = .norm_decostand(t_table, "log",         MARGIN %||% 1, logbase = logbase, ...),
    # --- 其他转换 ---
    "AST"         = .norm_ast(abund_table)
  )

  as.data.frame(res)
}


# ============================================================================
# 1. Rarefy (经典抽平)
# ============================================================================
.norm_rarefy <- function(abund_table, sample.size, rngseed, replace) {
  set.seed(rngseed)
  lib_sizes <- colSums(abund_table)

  if (is.null(sample.size)) {
    sample.size <- min(lib_sizes)
    message("Rarefy: 使用最小文库深度 = ", sample.size)
  }

  keep <- lib_sizes >= sample.size
  if (sum(keep) == 0) stop("所有样本文库深度均小于 sample.size = ", sample.size)
  if (sum(!keep) > 0)
    message("Rarefy: 移除 ", sum(!keep), " 个深度不足的样本: ",
            paste(colnames(abund_table)[!keep], collapse = ", "))

  abund_table <- abund_table[, keep, drop = FALSE]

  rarefied <- sapply(seq_len(ncol(abund_table)), function(j) {
    x <- as.integer(abund_table[, j])
    pool <- rep(seq_along(x), times = x)
    sampled <- sample(pool, size = sample.size, replace = replace)
    tabulate(sampled, nbins = length(x))
  })

  rownames(rarefied) <- rownames(abund_table)
  colnames(rarefied) <- colnames(abund_table)
  as.data.frame(rarefied)
}


# ============================================================================
# 2. SRS (Scaled Ranked Subsampling)
# ============================================================================
.norm_srs <- function(abund_table, sample.size, rngseed) {
  if (!requireNamespace("SRS", quietly = TRUE))
    stop("需要安装 SRS 包: install.packages('SRS')")

  if (is.null(sample.size)) {
    sample.size <- min(colSums(abund_table))
    message("SRS: 使用最小文库深度 = ", sample.size)
  }

  set.seed(rngseed)
  res <- SRS::SRS(abund_table, Cmin = sample.size)
  rownames(res) <- rownames(abund_table)
  res
}


# ============================================================================
# 3. TSS (Total Sum Scaling / 相对丰度)
# ============================================================================
.norm_tss <- function(abund_table) {
  lib_sizes <- colSums(abund_table)
  lib_sizes[lib_sizes == 0] <- 1
  sweep(abund_table, 2, lib_sizes, "/")
}


# ============================================================================
# 4. CSS (Cumulative Sum Scaling)
# ============================================================================
.norm_css <- function(abund_table, CSS_p, ...) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE))
    stop("需要安装 metagenomeSeq: BiocManager::install('metagenomeSeq')")
  if (!requireNamespace("Biobase", quietly = TRUE))
    stop("需要安装 Biobase: BiocManager::install('Biobase')")

  ct <- as.matrix(abund_table)
  pheno <- Biobase::AnnotatedDataFrame(
    data.frame(sample = colnames(ct), row.names = colnames(ct)))
  mgs <- metagenomeSeq::newMRexperiment(ct, phenoData = pheno)

  if (is.null(CSS_p)) {
    p_val <- metagenomeSeq::cumNormStatFast(mgs)
  } else {
    p_val <- CSS_p
  }
  mgs <- metagenomeSeq::cumNorm(mgs, p = p_val)
  res <- metagenomeSeq::MRcounts(mgs, norm = TRUE, log = FALSE)
  as.data.frame(res)
}


# ============================================================================
# 5. TMM / RLE (edgeR)
# ============================================================================
.norm_tmm_rle <- function(abund_table, norm_method, ...) {
  if (!requireNamespace("edgeR", quietly = TRUE))
    stop("需要安装 edgeR: BiocManager::install('edgeR')")

  ct <- as.matrix(round(abund_table))
  ct[ct < 0] <- 0

  dge <- edgeR::DGEList(counts = ct)
  dge <- edgeR::normLibSizes(dge, method = norm_method, ...)

  norm_factors <- dge$samples$norm.factors
  lib_sizes <- dge$samples$lib.size
  effective_lib <- lib_sizes * norm_factors
  scale_factors <- effective_lib / mean(effective_lib)

  res <- sweep(ct, 2, scale_factors, "/")
  as.data.frame(res)
}


# ============================================================================
# 6. DESeq2 Median-of-ratios
# ============================================================================
.norm_deseq2 <- function(abund_table, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("需要安装 DESeq2: BiocManager::install('DESeq2')")

  ct <- as.matrix(round(abund_table))
  ct[ct < 0] <- 0

  col_data <- data.frame(condition = rep("A", ncol(ct)), row.names = colnames(ct))
  dds <- DESeq2::DESeqDataSetFromMatrix(ct, col_data, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  res <- DESeq2::counts(dds, normalized = TRUE)
  as.data.frame(res)
}


# ============================================================================
# 7. GMPR (Geometric Mean of Pairwise Ratios)
# ============================================================================
.norm_gmpr <- function(abund_table, intersect.no, ct.min) {
  ct <- as.matrix(abund_table)
  n_samples <- ncol(ct)
  size_factors <- rep(NA_real_, n_samples)

  for (i in seq_len(n_samples)) {
    ratios_list <- numeric(0)
    for (j in seq_len(n_samples)) {
      if (i == j) next
      both_nonzero <- (ct[, i] >= ct.min) & (ct[, j] >= ct.min)
      if (sum(both_nonzero) < intersect.no) next
      r <- ct[both_nonzero, i] / ct[both_nonzero, j]
      ratios_list <- c(ratios_list, median(r))
    }
    if (length(ratios_list) > 0) {
      size_factors[i] <- exp(mean(log(ratios_list)))
    }
  }

  if (any(is.na(size_factors))) {
    warning("GMPR: ", sum(is.na(size_factors)), " 个样本无法计算 size factor, 使用 1 替代")
    size_factors[is.na(size_factors)] <- 1
  }

  res <- sweep(ct, 2, size_factors, "/")
  as.data.frame(res)
}


# ============================================================================
# 8. CLR (Centered Log-Ratio)
#    clr_ki = log(x_ki / g(x_i)), g = geometric mean
# ============================================================================
.norm_clr <- function(abund_table, pseudocount) {
  ct <- as.matrix(abund_table) + pseudocount

  res <- apply(ct, 2, function(x) {
    gm <- exp(mean(log(x)))
    log(x / gm)
  })

  rownames(res) <- rownames(abund_table)
  as.data.frame(res)
}


# ============================================================================
# 9. Robust CLR (零值保留为零, 只用非零值算几何平均)
#    rclr_ki = log(x_ki / g(x_i > 0)), 零值保持为 0
# ============================================================================
.norm_rclr <- function(abund_table) {
  ct <- as.matrix(abund_table)

  res <- apply(ct, 2, function(x) {
    nonzero <- x[x > 0]
    if (length(nonzero) == 0) return(x)
    gm <- exp(mean(log(nonzero)))
    out <- rep(0, length(x))
    out[x > 0] <- log(x[x > 0] / gm)
    out
  })

  rownames(res) <- rownames(abund_table)
  as.data.frame(res)
}


# ============================================================================
# 10. Wrench
# ============================================================================
.norm_wrench <- function(abund_table, condition, ...) {
  if (!requireNamespace("Wrench", quietly = TRUE))
    stop("需要安装 Wrench: BiocManager::install('Wrench')")
  if (is.null(condition))
    stop("Wrench 方法需要提供 condition 参数 (分组向量)")

  ct <- as.matrix(abund_table)
  W <- Wrench::wrench(ct, condition = condition, ...)
  norm_factors <- W$nf

  res <- sweep(ct, 2, norm_factors, "/")
  as.data.frame(res)
}


# ============================================================================
# 11. eBay (Empirical Bayes)
#     shrink sample proportions towards global mean
# ============================================================================
.norm_ebay <- function(abund_table) {
  ct <- as.matrix(abund_table)
  n_features <- nrow(ct)

  # TSS -> proportions
  col_totals <- colSums(ct)
  col_totals[col_totals == 0] <- 1
  props <- sweep(ct, 2, col_totals, "/")

  # global mean proportion per feature
  global_mean <- rowMeans(props)

  # EB shrinkage per sample
  res <- sapply(seq_len(ncol(ct)), function(j) {
    total_j <- col_totals[j]
    lambda <- n_features / (n_features + total_j)
    shrunk <- lambda * global_mean + (1 - lambda) * props[, j]
    shrunk / sum(shrunk)
  })

  rownames(res) <- rownames(abund_table)
  colnames(res) <- colnames(abund_table)
  as.data.frame(res)
}


# ============================================================================
# 12. vegan::decostand 封装
#     输入 t_table: 行=样本, 列=特征
#     输出转回: 行=特征, 列=样本
# ============================================================================
.norm_decostand <- function(t_table, method, MARGIN, ...) {
  if (!requireNamespace("vegan", quietly = TRUE))
    stop("需要安装 vegan 包: install.packages('vegan')")

  res <- vegan::decostand(t_table, method = method, MARGIN = MARGIN, ...)
  as.data.frame(t(res))
}


# ============================================================================
# 13. AST (Arc Sine Square Root)
#     要求输入为 [0, 1] 比例数据, 否则自动先做 TSS
# ============================================================================
.norm_ast <- function(abund_table) {
  ct <- as.matrix(abund_table)

  if (max(ct, na.rm = TRUE) > 1) {
    message("AST: 检测到非比例数据 (max > 1), 先转为相对丰度再做 AST 变换")
    lib_sizes <- colSums(ct)
    lib_sizes[lib_sizes == 0] <- 1
    ct <- sweep(ct, 2, lib_sizes, "/")
  }

  res <- asin(sqrt(ct))
  as.data.frame(res)
}


# ============================================================================
# 便捷函数: 批量标准化
# ============================================================================
#' 对同一数据批量运行多种标准化方法
#'
#' @param abund_table 丰度矩阵 (行=特征, 列=样本)
#' @param methods 方法向量
#' @return named list
#' @export
norm_batch <- function(abund_table,
                       methods = c("TSS", "clr", "rclr", "CSS", "hellinger", "AST")) {
  results <- list()
  for (m in methods) {
    cat("运行标准化:", m, "... ")
    results[[m]] <- tryCatch({
      res <- norm_abundance(abund_table, method = m)
      cat("完成\n")
      res
    }, error = function(e) {
      cat("失败:", e$message, "\n")
      NULL
    })
  }
  results
}


# ============================================================================
# 可视化: 标准化前后文库大小对比
# ============================================================================
#' @export
plot_norm_comparison <- function(raw, normalized, method_name = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要 ggplot2")

  lib_raw  <- data.frame(Sample = colnames(raw), LibSize = colSums(raw), Type = "Raw")
  lib_norm <- data.frame(Sample = colnames(normalized), LibSize = colSums(normalized), Type = method_name)
  lib_df <- rbind(lib_raw, lib_norm)
  lib_df$Type <- factor(lib_df$Type, levels = c("Raw", method_name))

  ggplot2::ggplot(lib_df, ggplot2::aes(x = reorder(Sample, LibSize), y = LibSize, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Sample", y = "Library Size / Sum",
                  title = paste0("Library Size: Raw vs ", method_name)) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = c("Raw" = "#e74c3c", "#3498db"))
}