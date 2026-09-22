#' ============================================================================
#' aggregate_taxa: 微生物丰度表分类学水平聚合函数
#' ============================================================================
#'
#' 将 OTU/ASV/Gene 丰度表按指定分类学水平 (Phylum/Genus/Species 等) 聚合。
#' 独立于 microeco 包，与同目录下 filter_abundance() / norm_abundance() /
#' diff_abundance() 衔接: 聚合 -> 过滤 -> 标准化 -> 差异分析。
#'
#' 输入: 行=特征(features)，列=样本(samples) 的丰度矩阵或 data.frame，
#' 以及行=特征、列=分类等级的分类注释表。两表通过行名 (Feature ID) 匹配。
#'
#' @param abund_table matrix 或 data.frame，行=Feature ID (OTU/ASV/Gene)，列=样本。
#'   行名必须是 Feature ID。
#' @param tax_table matrix 或 data.frame，行=Feature ID，列=分类等级
#'   (如 Kingdom, Phylum, ..., Species)。行名必须与 abund_table 的行名同源。
#' @param target_rank 字符串，目标聚合水平 (如 "Phylum" 或 "Genus")，
#'   必须是 tax_table 的列名之一。
#' @param na_action 字符串，如何处理分类缺失的特征:
#'   \describe{
#'     \item{"unclassified"}{(默认) 保留并重命名为 Unclassified 或利用上一级填充}
#'     \item{"remove"}{直接丢弃这些特征}
#'   }
#' @param fill_na_with_rank (可选) 字符串。当 target_rank 缺失时，
#'   尝试用上一级分类名填充，如 target="Species", fill="Genus" 时，
#'   缺失值将命名为 "Unclassified_属名"；若上一级也缺失则记为 "Unassigned"。
#' @param keep_all_abundance 布尔值 (默认 TRUE)。
#'   \describe{
#'     \item{TRUE}{保留丰度表中的所有特征；在 tax_table 中找不到的记为 Unclassified}
#'     \item{FALSE}{仅保留丰度表与 tax_table 中 Feature ID 均存在的行 (交集)}
#'   }
#' @param verbose 是否打印匹配与聚合过程信息，默认 TRUE。
#'
#' @return data.frame，行=聚合后的分类名 (target_rank 各水平)，列=样本，
#'   可直接传给 filter_abundance() / norm_abundance() / diff_abundance()。
#'
#' @examples
#' abund <- matrix(1:12, nrow = 4, dimnames = list(paste0("OTU", 1:4), paste0("S", 1:3)))
#' tax <- data.frame(
#'   Phylum = c("Firmicutes", "Firmicutes", NA, "Bacteroidetes"),
#'   row.names = paste0("OTU", 1:4)
#' )
#' aggregate_taxa(abund, tax, target_rank = "Phylum")
#'
#' @export

aggregate_taxa <- function(
  abund_table,
  tax_table,
  target_rank,
  na_action = c("unclassified", "remove"),
  fill_na_with_rank = NULL,
  keep_all_abundance = TRUE,
  verbose = TRUE
) {
  na_action <- match.arg(na_action)
  abund_table <- as.data.frame(abund_table)
  tax_table <- as.data.frame(tax_table)

  if (is.null(rownames(abund_table)) || is.null(rownames(tax_table))) {
    stop("abund_table 与 tax_table 都必须有行名 (Feature ID) 用于匹配")
  }
  if (!target_rank %in% colnames(tax_table)) {
    stop("错误: 在 tax_table 中找不到目标列 '", target_rank, "'")
  }
  if (
    !is.null(fill_na_with_rank) && !fill_na_with_rank %in% colnames(tax_table)
  ) {
    warning(
      "警告: 填充列 '",
      fill_na_with_rank,
      "' 不存在, 已忽略智能填充功能。"
    )
    fill_na_with_rank <- NULL
  }

  # ======================== ID 匹配健康检查 ========================
  feature_ids <- rownames(abund_table)
  tax_ids <- rownames(tax_table)
  common_ids <- intersect(feature_ids, tax_ids)
  only_in_abund <- setdiff(feature_ids, tax_ids)

  if (length(common_ids) == 0) {
    stop("abund_table 与 tax_table 之间没有任何匹配的 Feature ID")
  }

  if (verbose) {
    cat(">>> [aggregate_taxa] 数据匹配报告:\n")
    cat(sprintf("    丰度表特征数: %d\n", length(feature_ids)))
    cat(sprintf("    分类表特征数: %d\n", length(tax_ids)))
    cat(sprintf("    成功匹配特征数: %d\n", length(common_ids)))
    if (length(only_in_abund) == 0) {
      cat("    完美匹配: 所有丰度表特征均有分类注释。\n")
    }
  }

  if (length(only_in_abund) > 0) {
    if (keep_all_abundance) {
      warning(sprintf(
        "注意: 有 %d 个特征仅存在于丰度表中。策略: 保留并标记为 Unclassified。",
        length(only_in_abund)
      ))
    } else {
      warning(sprintf(
        "注意: 有 %d 个特征仅存在于丰度表中。策略: 已丢弃 (keep_all_abundance=FALSE)。",
        length(only_in_abund)
      ))
    }
  }

  # ======================== 交并集控制 ========================
  if (!keep_all_abundance) {
    keep_idx <- feature_ids %in% tax_ids
    abund_table <- abund_table[keep_idx, , drop = FALSE]
    feature_ids <- rownames(abund_table)
  }

  # ======================== 取分类等级并处理缺失 ========================
  m <- match(feature_ids, tax_ids)
  rank_vals <- as.character(tax_table[[target_rank]][m])
  is_missing <- is.na(rank_vals) | rank_vals == "" | rank_vals == "NA"

  if (na_action == "remove") {
    keep2 <- !is_missing
    if (verbose && sum(!keep2) > 0) {
      cat(sprintf("    已移除 %d 行分类缺失的数据。\n", sum(!keep2)))
    }
    abund_table <- abund_table[keep2, , drop = FALSE]
    rank_vals <- rank_vals[keep2]
  } else if (any(is_missing)) {
    if (!is.null(fill_na_with_rank)) {
      # 场景 A: 智能填充 (如 Unclassified_Bacteroides)
      higher_vals <- as.character(tax_table[[fill_na_with_rank]][m])
      is_higher_missing <- is.na(higher_vals) |
        higher_vals == "" |
        higher_vals == "NA"
      fill_vals <- paste0("Unclassified_", higher_vals)
      fill_vals[is_higher_missing] <- "Unassigned"
      rank_vals[is_missing] <- fill_vals[is_missing]
    } else {
      # 场景 B: 简单重命名
      rank_vals[is_missing] <- "Unclassified"
    }
  }

  # ======================== 聚合计算 ========================
  # rowsum() 按 rank_vals 分组对丰度矩阵逐样本 (列) 求和, 避免额外的 dplyr/tibble 依赖
  agg_mat <- rowsum(as.matrix(abund_table), group = rank_vals, reorder = FALSE)
  agg_df <- as.data.frame(agg_mat)

  if (verbose) {
    cat(sprintf(
      "    聚合完成: %d 个特征 -> %d 个 %s 水平分类。\n",
      nrow(abund_table),
      nrow(agg_df),
      target_rank
    ))
  }

  agg_df
}
