# ============================================================================
# beta_calc() — 通用多元距离计算与排序
# 开发日期：2026-04-06
# 依赖：vegan（必需）；phyloseq（可选，UniFrac 距离）
#
# 支持三种输入方式：
#   1. data + meta     通用"样本×特征"矩阵（转录组、代谢组等）
#   2. otu + meta      微生物组 OTU 表（特征×样本，函数内自动转置）
#   3. ps              phyloseq 对象（包含 OTU 表、metadata、进化树）
#
# 支持三种排序方法：
#   PCoA（默认）、NMDS、PCA
#
# 输出可直接对接 pcoa_plot()
# ============================================================================

library(vegan)

# ── 内部辅助函数 ──────────────────────────────────────────────────────────────

#' 解析三种输入方式，返回标准化 (mat, meta_df)
#' 统一规范：mat 为"样本×特征"，meta_df 行名为样本名，已对齐
#' @noRd
.resolve_input <- function(data, otu, tree, ps, meta, group) {

  # ── 路径1：phyloseq 对象 ────────────────────────────────────────────────────
  if (!is.null(ps)) {
    if (!requireNamespace("phyloseq", quietly = TRUE))
      stop("ps 输入需要 phyloseq 包：BiocManager::install('phyloseq')")

    otu_mat <- as.matrix(phyloseq::otu_table(ps))
    if (phyloseq::taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    mat     <- as.data.frame(otu_mat)   # 样本 × 特征
    meta_df <- as.data.frame(phyloseq::sample_data(ps))
    phylo_tree <- tryCatch(phyloseq::phy_tree(ps), error = function(e) NULL)

    if (!group %in% names(meta_df))
      stop(sprintf("列 '%s' 不在 phyloseq metadata 中", group))

    return(list(mat = mat, meta = meta_df, tree = phylo_tree, ps = ps))
  }

  # ── 路径2：OTU 表（特征×样本）─────────────────────────────────────────────
  if (!is.null(otu)) {
    if (is.null(meta)) stop("使用 otu 输入时需提供 meta 参数")
    mat     <- as.data.frame(t(as.matrix(otu)))   # 转置 → 样本×特征
    meta_df <- as.data.frame(meta)

    common <- intersect(rownames(mat), rownames(meta_df))
    if (length(common) == 0)
      stop("otu 列名与 meta 行名无交集，请检查样本名是否一致")
    mat     <- mat[common, , drop = FALSE]
    meta_df <- meta_df[common, , drop = FALSE]

    if (!group %in% names(meta_df))
      stop(sprintf("列 '%s' 不在 meta 中", group))

    return(list(mat = mat, meta = meta_df, tree = tree, ps = NULL))
  }

  # ── 路径3：通用矩阵（样本×特征）──────────────────────────────────────────
  if (!is.null(data)) {
    if (is.null(meta)) stop("使用 data 输入时需提供 meta 参数")
    mat     <- as.data.frame(data)
    meta_df <- as.data.frame(meta)

    if (!all(sapply(mat, is.numeric)))
      stop("data 中含非数值列，请仅传入数值型特征矩阵")

    common <- intersect(rownames(mat), rownames(meta_df))
    if (length(common) > 0) {
      mat     <- mat[common, , drop = FALSE]
      meta_df <- meta_df[common, , drop = FALSE]
    } else if (nrow(mat) == nrow(meta_df)) {
      message("行名不匹配，假设 data 与 meta 行顺序一致")
    } else {
      stop("data 行数（", nrow(mat), "）与 meta 行数（", nrow(meta_df),
           "）不一致，且行名无法对齐")
    }

    if (!group %in% names(meta_df))
      stop(sprintf("列 '%s' 不在 meta 中", group))

    return(list(mat = mat, meta = meta_df, tree = NULL, ps = NULL))
  }

  stop("必须提供 data、otu、ps 三者之一")
}

# ─────────────────────────────────────────────────────────────────────────────

#' 计算距离矩阵
#' @noRd
.calc_dist <- function(mat, dist_method, tree = NULL, ps = NULL) {

  unifrac_dists <- c("unifrac", "wunifrac", "dpcoa")

  if (dist_method %in% unifrac_dists) {
    if (is.null(ps))
      stop(sprintf(
        "'%s' 距离需要 phyloseq 对象（含进化树），请使用 ps 参数输入",
        dist_method
      ))
    if (!requireNamespace("phyloseq", quietly = TRUE))
      stop("需要 phyloseq 包：BiocManager::install('phyloseq')")
    ps_rela <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    return(phyloseq::distance(ps_rela, method = dist_method))
  }

  tryCatch(
    vegan::vegdist(mat, method = dist_method),
    error = function(e)
      stop(sprintf("vegdist 计算 '%s' 距离失败：%s", dist_method, e$message))
  )
}

# ─────────────────────────────────────────────────────────────────────────────

#' 运行排序，返回坐标 + 轴标签 + stress
#' @noRd
.do_ordination <- function(mat, dist_mat, method, ndim, scale, seed, verbose) {

  if (method == "PCoA") {
    res    <- cmdscale(dist_mat, k = ndim, eig = TRUE)
    coords <- as.data.frame(res$points)
    rownames(coords) <- rownames(mat)
    names(coords)    <- paste0("PC", seq_len(ndim))

    eig     <- res$eig
    pos_eig <- eig[eig > 0]
    var_pct <- pmax(round(eig[seq_len(ndim)] / sum(pos_eig) * 100, 2), 0)
    axis_labs <- sprintf("PCoA%d (%.1f%%)", seq_len(ndim), var_pct)
    return(list(coords = coords, axis_labs = axis_labs, stress = NULL))
  }

  if (method == "NMDS") {
    set.seed(seed)
    res    <- vegan::metaMDS(dist_mat, k = ndim, trace = if (verbose) 1L else 0L)
    coords <- as.data.frame(res$points)
    rownames(coords) <- rownames(mat)
    names(coords)    <- paste0("NMDS", seq_len(ndim))

    stress    <- round(res$stress, 4)
    axis_labs <- paste0("NMDS", seq_len(ndim))
    if (verbose) message("NMDS stress = ", stress)
    return(list(coords = coords, axis_labs = axis_labs, stress = stress))
  }

  if (method == "PCA") {
    res    <- prcomp(mat, center = TRUE, scale. = scale)
    coords <- as.data.frame(res$x[, seq_len(ndim), drop = FALSE])
    rownames(coords) <- rownames(mat)
    names(coords)    <- paste0("PC", seq_len(ndim))

    var_pct   <- round(summary(res)$importance[2L, seq_len(ndim)] * 100, 2)
    axis_labs <- sprintf("PC%d (%.1f%%)", seq_len(ndim), var_pct)
    return(list(coords = coords, axis_labs = axis_labs, stress = NULL))
  }
}

# ─────────────────────────────────────────────────────────────────────────────

#' 运行单次统计检验，返回标准化 data.frame
#' @noRd
.run_one_stat <- function(dist_mat, group_vec, stat_method, seed) {
  set.seed(seed)

  if (stat_method == "adonis2") {
    df  <- data.frame(group = group_vec)
    res <- vegan::adonis2(dist_mat ~ group, data = df, permutations = 999)
    return(data.frame(
      method         = "adonis2",
      statistic_name = "R2",
      statistic      = round(res$R2[1L], 4),
      p_value        = res$`Pr(>F)`[1L],
      stringsAsFactors = FALSE
    ))
  }

  if (stat_method == "anosim") {
    res <- vegan::anosim(dist_mat, group_vec, permutations = 999)
    return(data.frame(
      method         = "anosim",
      statistic_name = "R",
      statistic      = round(res$statistic, 4),
      p_value        = res$signif,
      stringsAsFactors = FALSE
    ))
  }

  if (stat_method == "mrpp") {
    res <- vegan::mrpp(dist_mat, group_vec, permutations = 999)
    return(data.frame(
      method         = "mrpp",
      statistic_name = "delta",
      statistic      = round(res$delta, 4),
      p_value        = res$Pvalue,
      stringsAsFactors = FALSE
    ))
  }

  NULL
}

# ─────────────────────────────────────────────────────────────────────────────

#' 将 stat_table 格式化为可直接传给 pcoa_plot(stat_text) 的字符向量
#' @noRd
.format_stat_text <- function(stat_table, stress = NULL) {
  lines <- character(0)

  if (!is.null(stress))
    lines <- c(lines, sprintf("Stress = %.4f", stress))

  if (!is.null(stat_table) && nrow(stat_table) > 0) {
    for (i in seq_len(nrow(stat_table))) {
      row   <- stat_table[i, ]
      p_str <- if (row$p_value < 0.001) "P < 0.001"
                else sprintf("P = %.3f", row$p_value)
      line  <- switch(row$method,
        "adonis2" = sprintf("Adonis: R\u00b2 = %.4f, %s", row$statistic, p_str),
        "anosim"  = sprintf("ANOSIM: R = %.4f, %s",       row$statistic, p_str),
        "mrpp"    = sprintf("MRPP: \u03b4 = %.4f, %s",    row$statistic, p_str),
        sprintf("%s = %.4f, %s", row$statistic_name, row$statistic, p_str)
      )
      lines <- c(lines, line)
    }
  }

  if (length(lines) == 0L) NULL else lines
}

# ── 主函数 ───────────────────────────────────────────────────────────────────

#' 通用多元距离计算与排序（Beta 多样性）
#'
#' 对"样本×特征"矩阵计算样本间距离并进行降维，输出可直接对接 \code{pcoa_plot()}。
#' 不仅限于微生物组：转录组、代谢组、蛋白质组、生态群落等任意"样本×特征"数据
#' 均可使用。
#'
#' @section 输入方式（三选一，优先级：ps > otu > data）：
#' \describe{
#'   \item{\code{data + meta}}{通用矩阵，样本为行，特征为列}
#'   \item{\code{otu + meta}}{微生物组 OTU 表，特征为行、样本为列（自动转置）}
#'   \item{\code{ps}}{phyloseq 对象，包含 OTU 表、metadata、进化树}
#' }
#'
#' @param data     样本×特征矩阵或 data.frame（通用输入，样本为行）
#' @param otu      特征×样本 OTU 表（微生物组输入，函数内自动转置）
#' @param tree     进化树（phylo 对象，计算 unifrac/wunifrac 时需要）
#' @param ps       phyloseq 对象（优先级最高）
#' @param meta     样本元数据 data.frame，行名为样本名
#' @param group    meta 中的分组列名（字符串）
#' @param method   降维方法：\code{"PCoA"}（默认）、\code{"NMDS"}、\code{"PCA"}
#' @param dist     距离类型，传给 \code{vegan::vegdist()}。
#'   常用：\code{"bray"}（默认）、\code{"euclidean"}、\code{"jaccard"}、
#'   \code{"manhattan"}、\code{"horn"}、\code{"canberra"}。
#'   微生物组专用（需 phyloseq 及进化树）：\code{"unifrac"}、\code{"wunifrac"}。
#'   \strong{method="PCA" 时此参数被忽略}（统计检验自动使用欧氏距离）。
#' @param ndim     保留的排序维度数，默认 2
#' @param stat     统计检验方法，可传入向量使用多种方法：
#'   \code{"adonis2"}（默认）、\code{"anosim"}、\code{"mrpp"}、\code{"none"}
#' @param pairwise 是否进行两两组间比较，默认 \code{FALSE}
#' @param scale    PCA 时是否对特征标准化（\code{scale.=TRUE}）；
#'   其他方法时无效，默认 \code{FALSE}
#' @param seed     NMDS 迭代和置换检验的随机种子，默认 42
#' @param verbose  是否打印过程信息，默认 \code{TRUE}
#'
#' @return 命名列表：
#' \describe{
#'   \item{\code{coords}}{data.frame，降维坐标 + 分组列；行名为样本名；
#'     可直接传给 \code{pcoa_plot()} 的 \code{data} 参数}
#'   \item{\code{dist_mat}}{dist 对象，距离矩阵
#'     （PCA 时为欧氏距离矩阵，仅用于统计）}
#'   \item{\code{stat_table}}{data.frame，整体检验结果；
#'     列：method / statistic_name / statistic / p_value}
#'   \item{\code{pair_table}}{data.frame，两两比较结果
#'     （\code{pairwise=FALSE} 时为 \code{NULL}）；
#'     列：group1 / group2 / method / statistic_name / statistic / p_value}
#'   \item{\code{stat_text}}{字符向量，格式化注释文字，
#'     可直接传给 \code{pcoa_plot()} 的 \code{stat_text} 参数}
#'   \item{\code{axis_labs}}{字符向量，格式化轴标签，
#'     如 \code{c("PCoA1 (39.9\%)", "PCoA2 (12.7\%)")}}
#'   \item{\code{method}}{字符串，降维方法}
#'   \item{\code{dist_method}}{字符串，距离类型}
#'   \item{\code{stress}}{数值，NMDS stress（非 NMDS 时为 \code{NULL}）}
#' }
#'
#' @examples
#' library(vegan)
#' data(dune)
#' data(dune.env)
#'
#' # 基本用法：PCoA + Bray-Curtis
#' res <- beta_calc(data = dune, meta = dune.env, group = "Management",
#'                  method = "PCoA", dist = "bray")
#' res$stat_text    # 直接传给 pcoa_plot()
#' res$axis_labs    # 格式化轴标签
#'
#' # 与 pcoa_plot() 对接
#' # pcoa_plot(res$coords, x = "PC1", y = "PC2", group.by = "Management",
#' #           xlab = res$axis_labs[1], ylab = res$axis_labs[2],
#' #           stat_text = res$stat_text)
#'
#' # PCA（转录组场景，先标准化）
#' expr <- matrix(rnorm(400), nrow = 40,
#'                dimnames = list(paste0("S", 1:40), paste0("Gene", 1:10)))
#' meta <- data.frame(group = rep(c("A", "B"), each = 20),
#'                    row.names = rownames(expr))
#' res_pca <- beta_calc(data = expr, meta = meta, group = "group",
#'                      method = "PCA", scale = TRUE)
#'
#' # 多种统计方法 + 两两比较
#' res2 <- beta_calc(data = dune, meta = dune.env, group = "Management",
#'                   stat = c("adonis2", "anosim"), pairwise = TRUE)
#' res2$pair_table
#'
#' @export
beta_calc <- function(
  data     = NULL,
  otu      = NULL,
  tree     = NULL,
  ps       = NULL,
  meta     = NULL,
  group    = NULL,
  method   = c("PCoA", "NMDS", "PCA"),
  dist     = "bray",
  ndim     = 2L,
  stat     = "adonis2",
  pairwise = FALSE,
  scale    = FALSE,
  seed     = 42L,
  verbose  = TRUE
) {

  # ── 参数验证 ────────────────────────────────────────────────────────────────
  method <- match.arg(method)
  if (is.null(group)) stop("必须指定 group 参数（分组列名）")

  valid_stats <- c("adonis2", "anosim", "mrpp", "none")
  bad_stat <- setdiff(stat, valid_stats)
  if (length(bad_stat) > 0)
    stop(sprintf("无效的 stat 选项：%s\n  可选：%s",
                 paste(bad_stat, collapse = ", "),
                 paste(valid_stats, collapse = ", ")))

  # ── 输入解析 ────────────────────────────────────────────────────────────────
  if (verbose) message("── 解析输入...")
  inp     <- .resolve_input(data, otu, tree, ps, meta, group)
  mat     <- inp$mat
  meta_df <- inp$meta

  # 缺失值处理
  na_rows <- apply(mat, 1L, function(r) any(is.na(r)))
  if (any(na_rows)) {
    message("移除 ", sum(na_rows), " 个含 NA 的样本")
    mat     <- mat[!na_rows, , drop = FALSE]
    meta_df <- meta_df[!na_rows, , drop = FALSE]
  }
  if (nrow(mat) < 2L) stop("有效样本数不足（< 2），无法计算")

  # 分组因子化
  meta_df[[group]] <- factor(meta_df[[group]])
  group_vec <- meta_df[[group]]
  n_grp     <- nlevels(group_vec)

  if (verbose)
    message(sprintf("   样本：%d  |  分组：%d  |  特征：%d",
                    nrow(mat), n_grp, ncol(mat)))

  # ── 距离计算 ────────────────────────────────────────────────────────────────
  if (method == "PCA") {
    if (verbose)
      message("── PCA 模式：排序在原始矩阵上进行；统计检验自动使用欧氏距离")
    # 用与 PCA 一致的缩放做欧氏距离
    scaled_mat <- scale(mat, center = TRUE, scale = scale)
    dist_mat   <- dist(scaled_mat)
    dist_used  <- "euclidean (PCA)"
  } else {
    if (verbose) message("── 计算距离矩阵（", dist, "）...")
    dist_mat  <- .calc_dist(mat, dist, inp$tree, inp$ps)
    dist_used <- dist
  }

  # ── 排序 ────────────────────────────────────────────────────────────────────
  if (verbose) message("── 运行 ", method, " 排序...")
  ord    <- .do_ordination(mat, dist_mat, method, ndim, scale, seed, verbose)
  coords <- ord$coords
  coords[[group]] <- group_vec   # 添加分组列（行序与 mat 对齐）

  # ── 统计检验 ────────────────────────────────────────────────────────────────
  stat_table <- NULL
  pair_table <- NULL

  run_stat <- !("none" %in% stat)

  if (run_stat) {
    if (verbose) message("── 统计检验（", paste(stat, collapse = ", "), "）...")

    stat_rows <- lapply(stat, function(s) {
      tryCatch(
        .run_one_stat(dist_mat, group_vec, s, seed),
        error = function(e) {
          warning(sprintf("'%s' 检验失败：%s", s, e$message))
          NULL
        }
      )
    })
    stat_table <- do.call(rbind, Filter(Negate(is.null), stat_rows))

    # 两两比较
    if (pairwise && n_grp >= 2L) {
      if (verbose) message("── 两两比较...")
      grp_levels  <- levels(group_vec)
      pair_combns <- combn(grp_levels, 2L, simplify = FALSE)

      pair_rows <- lapply(pair_combns, function(pair) {
        idx      <- group_vec %in% pair
        sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
        sub_grp  <- droplevels(group_vec[idx])

        rows <- lapply(stat, function(s) {
          tryCatch({
            r          <- .run_one_stat(sub_dist, sub_grp, s, seed)
            r$group1   <- pair[1L]
            r$group2   <- pair[2L]
            r
          }, error = function(e) NULL)
        })
        do.call(rbind, Filter(Negate(is.null), rows))
      })

      pair_table <- do.call(rbind, Filter(Negate(is.null), pair_rows))
      if (!is.null(pair_table) && nrow(pair_table) > 0L) {
        pair_table <- pair_table[, c("group1", "group2", "method",
                                     "statistic_name", "statistic", "p_value")]
        rownames(pair_table) <- NULL
      }
    }
  }

  # ── 格式化输出 ──────────────────────────────────────────────────────────────
  stat_text <- .format_stat_text(stat_table, ord$stress)

  if (verbose && !is.null(stat_text)) {
    message("── 统计结果：")
    for (txt in stat_text) message("   ", txt)
  }

  # ── 返回 ────────────────────────────────────────────────────────────────────
  list(
    coords      = coords,
    dist_mat    = dist_mat,
    stat_table  = stat_table,
    pair_table  = pair_table,
    stat_text   = stat_text,
    axis_labs   = ord$axis_labs,
    method      = method,
    dist_method = dist_used,
    stress      = ord$stress
  )
}
