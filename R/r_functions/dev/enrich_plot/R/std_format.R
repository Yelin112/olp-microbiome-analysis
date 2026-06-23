# ==============================================================================
# EnrichStd: 富集分析标准化数据格式
# 创建: 2026-04-24
#
# 职责:
#   将各种来源的富集分析结果统一转换为内部标准格式 (EnrichStd)，
#   供所有可视化函数使用。
#   数据转换与可视化解耦，新增数据来源只需扩展 as_enrich_std()。
# ==============================================================================

# ------------------------------------------------------------------------------
# EnrichStd 标准格式说明
# ------------------------------------------------------------------------------
#
# EnrichStd 是一个带有 "EnrichStd" 类标记的 data.frame，包含以下列：
#
# 必需列：
#   term_id     <chr>  通路/term 的唯一标识符
#                      (如 "hsa04064", "GO:0006954", "K00001", 不限ID类型)
#   term_name   <chr>  通路/term 的显示名称
#   gene_ids    <chr>  富集到的基因/特征，"/" 分隔字符串
#                      (基因symbol/ENTREZID/KO号/OTU等均可，与上游一致)
#   p_value     <dbl>  原始 p 值
#   p_adjust    <dbl>  校正后 p 值
#   gene_count  <int>  富集到的基因/特征数量
#
# 可选列（由转换函数尽力填充）：
#   bg_count    <int>  该 term 的背景基因总数
#   gene_ratio  <dbl>  gene_count / 输入基因总数（气泡图 x 轴默认值）
#   rich_factor <dbl>  gene_count / bg_count（富集因子）
#
# 额外列：
#   转换过程中原始数据的其他列会原样保留，不会丢失。

# ------------------------------------------------------------------------------
# 标准列名：自动识别的候选列名映射
# ------------------------------------------------------------------------------

.ENRICH_COL_CANDIDATES <- list(
  term_id = c(
    "ID",
    "id",
    "term_id",
    "pathway_id",
    "KEGGID",
    "go_id"
  ),
  term_name = c(
    "Description",
    "description",
    "term_name",
    "pathway",
    "Pathway",
    "pathway_name",
    "name",
    "Name",
    "Term",
    "term"
  ),
  gene_ids = c(
    "geneID",
    "gene_ids",
    "genes",
    "gene",
    "GeneID",
    "gene_list",
    "core_enrichment"
  ),
  p_value = c(
    "pvalue",
    "p_value",
    "pval",
    "PValue",
    "p.value",
    "Pvalue"
  ),
  p_adjust = c(
    "p.adjust",
    "p_adjust",
    "padj",
    "FDR",
    "fdr",
    "BH",
    "qvalue",
    "qval",
    "p.adj"
  ),
  gene_count = c(
    "Count",
    "count",
    "gene_count",
    "overlap",
    "Overlap",
    "k"
  ),
  bg_count = c(
    "bg_count",
    "BgCount",
    "background",
    "M",
    "setSize",
    "set_size"
  )
)


# ------------------------------------------------------------------------------
# as_enrich_std(): 通用入口（S3 泛型）
# ------------------------------------------------------------------------------

#' 将富集分析结果转换为 EnrichStd 标准格式
#'
#' @param x 输入数据，支持以下类型：
#'   - enrichResult  (clusterProfiler::enrichKEGG/enrichGO 的结果)
#'   - gseaResult    (clusterProfiler::gseKEGG/gseGO 的结果)
#'   - data.frame    (需要满足列要求，或通过 col_map 指定列名)
#' @param col_map 列名映射，named list，仅在自动识别失败时使用。
#'   格式：list(标准列名 = "实际列名")
#'   可用的标准列名：term_id, term_name, gene_ids, p_value, p_adjust,
#'                   gene_count, bg_count
#'   示例：col_map = list(term_name = "my_pathway_col", gene_ids = "gene_list")
#' @param n_input_genes 输入基因总数，用于计算 gene_ratio。
#'   对于 enrichResult 对象会自动提取，data.frame 输入时若提供则更准确。
#' @param ... 传递给具体 S3 方法的其他参数
#'
#' @return 带有 "EnrichStd" 类标记的 data.frame
#'
#' @examples
#' \dontrun{
#' # clusterProfiler 结果直接转换（无需任何映射）
#' std <- as_enrich_std(kegg_result)
#'
#' # 自定义 data.frame（列名不标准时）
#' std <- as_enrich_std(my_df, col_map = list(
#'   term_name = "pathway_name",
#'   gene_ids  = "gene_list"
#' ))
#' }
#'
#' @export
as_enrich_std <- function(x, col_map = NULL, n_input_genes = NULL, ...) {
  UseMethod("as_enrich_std")
}


# ------------------------------------------------------------------------------
# as_enrich_std.enrichResult(): clusterProfiler enrichKEGG/enrichGO 结果
# ------------------------------------------------------------------------------

#' @export
#' @rdname as_enrich_std
as_enrich_std.enrichResult <- function(
  x,
  col_map = NULL,
  n_input_genes = NULL,
  ...
) {
  df <- as.data.frame(x)

  # 提取输入基因总数（用于 gene_ratio 分母）
  if (is.null(n_input_genes)) {
    n_input_genes <- length(x@gene)
  }

  # BgRatio 格式: "M/N"，提取 bg_count (M) 和背景总数 (N)
  if ("BgRatio" %in% colnames(df)) {
    bg_parts <- strsplit(df$BgRatio, "/")
    df$bg_count <- as.integer(sapply(bg_parts, `[`, 1))
  }

  # 直接映射标准列（enrichResult 的列名是固定的）
  std_df <- data.frame(
    term_id = .get_col(df, "term_id", col_map, "ID"),
    term_name = .get_col(df, "term_name", col_map, "Description"),
    gene_ids = .get_col(df, "gene_ids", col_map, "geneID"),
    p_value = .get_col(df, "p_value", col_map, "pvalue"),
    p_adjust = .get_col(df, "p_adjust", col_map, "p.adjust"),
    gene_count = as.integer(.get_col(df, "gene_count", col_map, "Count")),
    stringsAsFactors = FALSE
  )

  # 填充可选列
  if ("bg_count" %in% colnames(df)) {
    std_df$bg_count <- df$bg_count
  }

  # 计算 gene_ratio 和 rich_factor
  std_df <- .compute_derived_cols(std_df, n_input_genes)

  # 保留原始数据中其他列（附加在末尾，不覆盖标准列）
  # 注意：已映射为标准列的原始列名、以及函数内部新增的列都要排除
  extra_cols <- setdiff(
    colnames(df),
    c(
      # enrichResult 原始列
      "ID", "Description", "geneID", "pvalue", "p.adjust",
      "Count", "BgRatio", "GeneRatio", "qvalue", "RichFactor",
      # 函数内部新增列（避免重复）
      "bg_count"
    )
  )
  if (length(extra_cols) > 0) {
    std_df <- cbind(std_df, df[extra_cols])
  }

  .new_enrich_std(std_df)
}


# ------------------------------------------------------------------------------
# as_enrich_std.gseaResult(): clusterProfiler GSEA 结果
# ------------------------------------------------------------------------------

#' @export
#' @rdname as_enrich_std
as_enrich_std.gseaResult <- function(
  x,
  col_map = NULL,
  n_input_genes = NULL,
  ...
) {
  df <- as.data.frame(x)

  # GSEA 结果使用 core_enrichment 而非 geneID
  std_df <- data.frame(
    term_id = .get_col(df, "term_id", col_map, "ID"),
    term_name = .get_col(df, "term_name", col_map, "Description"),
    gene_ids = .get_col(df, "gene_ids", col_map, "core_enrichment"),
    p_value = .get_col(df, "p_value", col_map, "pvalue"),
    p_adjust = .get_col(df, "p_adjust", col_map, "p.adjust"),
    gene_count = as.integer(
      sapply(
        strsplit(.get_col(df, "gene_ids", col_map, "core_enrichment"), "/"),
        length
      )
    ),
    stringsAsFactors = FALSE
  )

  # GSEA 特有列
  if ("setSize" %in% colnames(df)) {
    std_df$bg_count <- df$setSize
  }
  if ("NES" %in% colnames(df)) {
    std_df$NES <- df$NES
  }

  if (is.null(n_input_genes)) {
    n_input_genes <- tryCatch(length(x@geneList), error = function(e) NULL)
  }

  std_df <- .compute_derived_cols(std_df, n_input_genes)

  extra_cols <- setdiff(
    colnames(df),
    c(
      "ID",
      "Description",
      "core_enrichment",
      "pvalue",
      "p.adjust",
      "setSize",
      "NES",
      "GeneRatio",
      "enrichmentScore"
    )
  )
  if (length(extra_cols) > 0) {
    std_df <- cbind(std_df, df[extra_cols])
  }

  .new_enrich_std(std_df)
}


# ------------------------------------------------------------------------------
# as_enrich_std.data.frame(): 自定义 data.frame 输入
# ------------------------------------------------------------------------------

#' @export
#' @rdname as_enrich_std
as_enrich_std.data.frame <- function(
  x,
  col_map = NULL,
  n_input_genes = NULL,
  ...
) {
  df <- x

  # 对每个标准列：先查 col_map，再自动识别，再报错
  resolved <- list()
  missing_required <- character(0)

  required_cols <- c(
    "term_name",
    "gene_ids",
    "p_value",
    "p_adjust",
    "gene_count"
  )
  optional_cols <- c("term_id", "bg_count")

  all_cols <- c(required_cols, optional_cols)

  for (col in all_cols) {
    result <- .resolve_col(df, col, col_map)

    if (is.null(result)) {
      if (col %in% required_cols) {
        missing_required <- c(missing_required, col)
      }
      # 可选列找不到时跳过
    } else {
      resolved[[col]] <- result
    }
  }

  # 报告缺失的必需列
  if (length(missing_required) > 0) {
    candidates_hint <- sapply(missing_required, function(col) {
      cands <- .ENRICH_COL_CANDIDATES[[col]]
      paste0("  - ", col, ": 可识别的列名包括 ", paste(cands, collapse = ", "))
    })
    stop(
      "缺少必需列，请通过 col_map 指定:\n",
      paste(candidates_hint, collapse = "\n"),
      "\n\n",
      "示例: as_enrich_std(df, col_map = list(term_name = '你的列名', ...))\n",
      "当前数据列名: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  # 构建标准 data.frame
  std_df <- data.frame(
    term_id = if (!is.null(resolved$term_id)) {
      df[[resolved$term_id]]
    } else {
      paste0("term_", seq_len(nrow(df)))
    },
    term_name = df[[resolved$term_name]],
    gene_ids = df[[resolved$gene_ids]],
    p_value = as.numeric(df[[resolved$p_value]]),
    p_adjust = as.numeric(df[[resolved$p_adjust]]),
    gene_count = as.integer(df[[resolved$gene_count]]),
    stringsAsFactors = FALSE
  )

  if (!is.null(resolved$bg_count)) {
    std_df$bg_count <- as.integer(df[[resolved$bg_count]])
  }

  # 尝试从 GeneRatio 字符列（"k/n" 格式）解析数值型 gene_ratio
  # 常见于 clusterProfiler 导出为 data.frame 后再传入的场景
  gene_ratio_src <- NULL
  if ("GeneRatio" %in% colnames(df) && is.character(df[["GeneRatio"]])) {
    parts <- strsplit(df[["GeneRatio"]], "/")
    ratios <- sapply(parts, function(p) {
      if (length(p) < 2) return(NA_real_)
      k     <- suppressWarnings(as.numeric(p[1]))
      n_val <- suppressWarnings(as.numeric(p[2]))
      if (!is.na(k) && !is.na(n_val) && n_val > 0) k / n_val else NA_real_
    })
    if (any(!is.na(ratios))) {
      std_df$gene_ratio <- ratios
      gene_ratio_src <- "GeneRatio"
    }
  }

  std_df <- .compute_derived_cols(std_df, n_input_genes)

  # 保留其他未参与映射的列（已解析的 GeneRatio 字符列不再附加）
  used_src_cols <- c(unlist(resolved), gene_ratio_src)
  extra_cols <- setdiff(colnames(df), used_src_cols)
  if (length(extra_cols) > 0) {
    std_df <- cbind(std_df, df[extra_cols])
  }

  .new_enrich_std(std_df)
}


# ------------------------------------------------------------------------------
# as_enrich_std.EnrichStd(): 已经是标准格式，直接返回
# ------------------------------------------------------------------------------

#' @export
#' @rdname as_enrich_std
as_enrich_std.EnrichStd <- function(
  x,
  col_map = NULL,
  n_input_genes = NULL,
  ...
) {
  x
}


# ------------------------------------------------------------------------------
# 内部辅助函数
# ------------------------------------------------------------------------------

# 构造 EnrichStd 对象，附加类标记并验证
.new_enrich_std <- function(df) {
  # 确保必需列存在且类型正确
  stopifnot(
    "term_name" %in% colnames(df),
    "gene_ids" %in% colnames(df),
    "p_value" %in% colnames(df),
    "p_adjust" %in% colnames(df),
    "gene_count" %in% colnames(df)
  )

  class(df) <- c("EnrichStd", "data.frame")
  df
}

# 计算 gene_ratio 和 rich_factor
.compute_derived_cols <- function(std_df, n_input_genes) {
  # gene_ratio: gene_count / n_input_genes
  if (!is.null(n_input_genes) && n_input_genes > 0) {
    std_df$gene_ratio <- std_df$gene_count / n_input_genes
  } else if (!"gene_ratio" %in% colnames(std_df)) {
    std_df$gene_ratio <- NA_real_
  }

  # rich_factor: gene_count / bg_count
  if ("bg_count" %in% colnames(std_df)) {
    std_df$rich_factor <- ifelse(
      std_df$bg_count > 0,
      std_df$gene_count / std_df$bg_count,
      NA_real_
    )
  } else if (!"rich_factor" %in% colnames(std_df)) {
    std_df$rich_factor <- NA_real_
  }

  std_df
}

# 从 df 中获取某列：先查 col_map，再用默认列名
.get_col <- function(df, std_name, col_map, default_name) {
  actual_name <- if (!is.null(col_map) && std_name %in% names(col_map)) {
    col_map[[std_name]]
  } else {
    default_name
  }

  if (!actual_name %in% colnames(df)) {
    stop(
      "找不到列 '",
      actual_name,
      "'（对应标准字段: ",
      std_name,
      "）\n",
      "当前数据列名: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  df[[actual_name]]
}

# 自动识别列名（用于 data.frame 方法）
# 返回实际列名，找不到则返回 NULL
.resolve_col <- function(df, std_name, col_map) {
  # 1. 先查 col_map（用户明确指定）
  if (!is.null(col_map) && std_name %in% names(col_map)) {
    user_col <- col_map[[std_name]]
    if (!user_col %in% colnames(df)) {
      stop(
        "col_map 指定的列 '",
        user_col,
        "' 不存在于数据中\n",
        "当前数据列名: ",
        paste(colnames(df), collapse = ", ")
      )
    }
    return(user_col)
  }

  # 2. 自动识别候选列名
  candidates <- .ENRICH_COL_CANDIDATES[[std_name]]
  match_col <- intersect(candidates, colnames(df))

  if (length(match_col) > 0) {
    if (length(match_col) > 1) {
      message(
        "列 '",
        std_name,
        "' 匹配到多个候选列: ",
        paste(match_col, collapse = ", "),
        "，使用第一个: ",
        match_col[1]
      )
    }
    return(match_col[1])
  }

  NULL # 未找到
}


# 检测是否为双组 list 输入（bar / lollipop / scatter 共用）
.is_dual_input <- function(x) {
  is.list(x) && !inherits(x, c("data.frame", "EnrichStd", "enrichResult"))
}


# ------------------------------------------------------------------------------
# 内部：计算 x 轴数值列（供 bar/lollipop/combined 共用）
# x_var: "p_adjust" → -log10(p_adjust)；"gene_ratio" → gene_ratio；"gene_count" → gene_count
# ------------------------------------------------------------------------------

.compute_x_col <- function(df, x_var = "p_adjust") {
  if (x_var == "gene_ratio") {
    col <- if ("gene_ratio" %in% colnames(df)) df$gene_ratio else NULL
    if (is.null(col) || all(is.na(col))) {
      warning("gene_ratio 不可用，退回 -log10(p_adjust)")
      return(-log10(df$p_adjust))
    }
    return(col)
  }
  if (x_var == "gene_count") {
    return(as.numeric(df$gene_count))
  }
  -log10(df$p_adjust)
}

# ------------------------------------------------------------------------------
# 工具函数：合并两组富集结果为双向发散格式（供 bar/lollipop 使用）
# ------------------------------------------------------------------------------

#' 将两组富集结果合并为双向发散图所需的数据格式
#'
#' 用于双向 bar/lollipop 图，将两个富集结果（如 Up/Down 组）合并为单一
#' data.frame，正负方向分别代表两组，x 轴为 -log10(p.adjust) * direction。
#'
#' @param data_list 长度为 2 的命名 list，每个元素为 enrichResult / EnrichStd / data.frame。
#'   名称作为分组标签，如 list(Up = kegg_up, Down = kegg_down)。
#' @param top_n    每组各取前 N 个 term，默认 10（双向共最多 2*top_n 个 term）
#' @param filter_padj p_adjust 过滤阈值，NULL 表示不过滤
#' @param col_map  列名映射，传递给 as_enrich_std()
#'
#' @return data.frame，包含标准列 + 以下新增列：
#'   - group:     分组标签（factor，levels 按 data_list 名称顺序）
#'   - direction: +1 或 -1（第一组为 +1，第二组为 -1）
#'   - x:         -log10(p_adjust) * direction（x 轴值）
#'
#' @export
prep_enrich_dual <- function(
  data_list,
  top_n       = 10,
  filter_padj = NULL,
  col_map     = NULL,
  x_var       = "p_adjust"
) {
  if (!is.list(data_list) || length(data_list) != 2) {
    stop("data_list 必须是长度为 2 的 list")
  }
  if (is.null(names(data_list)) || any(nchar(names(data_list)) == 0)) {
    stop("data_list 的两个元素必须有名称，如 list(Up = ..., Down = ...)")
  }

  group_names <- names(data_list)

  parts <- lapply(seq_along(data_list), function(i) {
    std <- as_enrich_std(data_list[[i]], col_map = col_map)
    std <- prep_enrich_data(std, top_n = top_n, filter_padj = filter_padj)
    std$group     <- group_names[i]
    std$direction <- if (i == 1) 1L else -1L
    std
  })

  combined <- do.call(rbind, parts)
  combined$group <- factor(combined$group, levels = group_names)
  combined$x     <- .compute_x_col(combined, x_var) * combined$direction

  # y 轴排序：按 x 值从小到大排列（发散图左右对称）
  combined <- combined[order(combined$x), ]
  combined$term_name <- factor(combined$term_name, levels = unique(combined$term_name))

  rownames(combined) <- NULL
  combined
}


# ------------------------------------------------------------------------------
# 工具函数：合并多分类富集结果（供 combined 样式使用）
# ------------------------------------------------------------------------------

#' 将多组富集结果合并为带分类标签的统一格式
#'
#' 用于 plot_enrich_combined()，支持任意分类（GO本体/KEGG/自定义数据库）。
#'
#' @param data 输入数据，两种形式：
#'   - **named list**：每个元素为 enrichResult / EnrichStd / data.frame，
#'     名称作为分类标签，如 `list(BP = go_bp, KEGG = kegg_res)`
#'   - **data.frame**：已合并的数据，需同时指定 `ontology_col`
#' @param ontology_col 分类列名，仅 data.frame 输入时有效
#' @param col_map 列名映射，传给 as_enrich_std()
#' @param top_n 每个分类取前 N 个 term（按 p_adjust 升序），默认 5
#' @param category_order 分类排列顺序（决定 y 轴从下到上的顺序），
#'   NULL 时使用输入顺序的逆序（与模板一致）
#'
#' @return data.frame，包含标准列 + 新增列：
#'   - category: 分类标签（factor）
#'   - index:    行序号（用于 y 轴数值定位）
#'
#' @export
prep_enrich_combined <- function(
  data,
  ontology_col   = NULL,
  col_map        = NULL,
  top_n          = 5,
  category_order = NULL
) {
  # 接受 named list 输入
  if (is.list(data) && !is.data.frame(data) &&
      !inherits(data, "enrichResult") && !inherits(data, "EnrichStd")) {
    if (is.null(names(data)) || any(nchar(names(data)) == 0)) {
      stop("list 输入的每个元素必须有名称（作为分类标签）")
    }
    parts <- lapply(names(data), function(cat_name) {
      std <- as_enrich_std(data[[cat_name]], col_map = col_map)
      std$category <- cat_name
      std
    })
    # 各分类可能有不同的额外列（如 GO 有 ONTOLOGY 而 KEGG 没有），
    # 对齐所有列后再合并，缺失列补 NA
    all_cols <- unique(unlist(lapply(parts, colnames)))
    parts <- lapply(parts, function(df) {
      missing <- setdiff(all_cols, colnames(df))
      if (length(missing) > 0) {
        df[missing] <- NA
      }
      df[, all_cols, drop = FALSE]
    })
    combined <- do.call(rbind, parts)

  # 接受 data.frame 输入（含分类列）
  } else {
    combined <- as_enrich_std(data, col_map = col_map)
    if (is.null(ontology_col)) {
      stop("data.frame 输入时需指定 ontology_col（分类列名）")
    }
    if (!ontology_col %in% colnames(combined)) {
      stop("列 '", ontology_col, "' 不存在，当前列名: ",
           paste(colnames(combined), collapse = ", "))
    }
    combined$category <- combined[[ontology_col]]
  }

  # 每个分类取 top_n（按 p_adjust 升序）
  combined <- do.call(rbind, lapply(
    split(combined, combined$category),
    function(grp) {
      grp <- grp[order(grp$p_adjust), , drop = FALSE]
      head(grp, top_n)
    }
  ))

  if (nrow(combined) == 0) {
    stop("筛选后所有分类均无数据，请检查输入数据或放宽 top_n / filter_padj")
  }

  # 分类排列顺序
  all_cats <- unique(as.character(combined$category))
  if (is.null(category_order)) {
    category_order <- rev(all_cats)  # 逆序：最后一个分类在图的最下方
  } else {
    missing_cats <- setdiff(all_cats, category_order)
    if (length(missing_cats) > 0) {
      warning("category_order 未包含以下分类，将附加到末尾: ",
              paste(missing_cats, collapse = ", "))
      category_order <- c(category_order, missing_cats)
    }
    category_order <- intersect(category_order, all_cats)
  }

  combined$category <- factor(combined$category, levels = category_order)

  # 排序：按 category（factor 顺序）然后 p_adjust
  combined <- combined[order(combined$category, combined$p_adjust), , drop = FALSE]

  # 添加行索引（y 轴数值定位使用）
  combined$index <- seq_len(nrow(combined))

  # 固定 term_name 显示顺序
  combined$term_name <- factor(combined$term_name, levels = combined$term_name)

  rownames(combined) <- NULL
  combined
}


# ------------------------------------------------------------------------------
# 工具函数：展开 gene_ids 为长格式（供热图/网络图使用）
# ------------------------------------------------------------------------------

#' 将 EnrichStd 的 gene_ids 列展开为长格式
#'
#' @param std_df EnrichStd 对象
#' @param gene_col gene_ids 列名，默认 "gene_ids"
#' @param sep 分隔符，默认 "/"
#' @return data.frame，每行一个 gene-term 对
#'
#' @export
expand_gene_ids <- function(std_df, gene_col = "gene_ids", sep = "/") {
  if (!gene_col %in% colnames(std_df)) {
    stop("列 '", gene_col, "' 不存在")
  }

  rows <- lapply(seq_len(nrow(std_df)), function(i) {
    genes <- unlist(strsplit(std_df[[gene_col]][i], sep, fixed = TRUE))
    genes <- trimws(genes)
    genes <- genes[nchar(genes) > 0]
    if (length(genes) == 0) {
      return(NULL)
    }
    row_data <- std_df[rep(i, length(genes)), , drop = FALSE]
    row_data[[gene_col]] <- genes
    row_data
  })

  result <- do.call(rbind, rows[!sapply(rows, is.null)])
  rownames(result) <- NULL
  # 展开后是长格式，不再是合法的 EnrichStd，去掉类标记避免 print 误触发
  class(result) <- "data.frame"
  result
}


# ------------------------------------------------------------------------------
# 工具函数：筛选和预处理（可视化函数的公共前处理）
# ------------------------------------------------------------------------------

#' 对 EnrichStd 数据进行筛选和排序，供绘图函数调用
#'
#' @param std_df EnrichStd 对象
#' @param top_n 选取前 N 个 term（按 p_adjust 升序），NULL 表示不限制
#' @param filter_padj p_adjust 阈值过滤，NULL 表示不过滤
#' @param terms 手动指定 term_name 向量，会覆盖 top_n
#' @param order_by 排序列，默认 "p_adjust"；可选 "gene_count", "rich_factor", "gene_ratio"
#' @param decreasing 排序方向，默认 FALSE（升序）
#'
#' @return 筛选排序后的 EnrichStd 对象
#'
#' @export
prep_enrich_data <- function(
  std_df,
  top_n = 10,
  filter_padj = NULL,
  terms = NULL,
  order_by = "p_adjust",
  decreasing = FALSE
) {
  std_df <- as_enrich_std(std_df) # 确保输入合法

  # p_adjust 过滤
  if (!is.null(filter_padj)) {
    std_df <- std_df[std_df$p_adjust <= filter_padj, , drop = FALSE]
    if (nrow(std_df) == 0) {
      stop(
        "过滤后没有剩余的 term，请放宽 filter_padj 阈值（当前: ",
        filter_padj,
        "）"
      )
    }
  }

  # 手动指定 terms
  if (!is.null(terms)) {
    std_df <- std_df[std_df$term_name %in% terms, , drop = FALSE]
    not_found <- setdiff(terms, std_df$term_name)
    if (length(not_found) > 0) {
      warning("以下 term 未找到: ", paste(not_found, collapse = ", "))
    }
    if (nrow(std_df) == 0) {
      stop("指定的 terms 均未找到，请检查 term_name 拼写")
    }
    return(.new_enrich_std(std_df))
  }

  # 排序
  if (order_by %in% colnames(std_df)) {
    std_df <- std_df[
      order(std_df[[order_by]], decreasing = decreasing),
      ,
      drop = FALSE
    ]
  } else {
    warning("排序列 '", order_by, "' 不存在，跳过排序")
  }

  # 取 top_n
  if (!is.null(top_n) && top_n < nrow(std_df)) {
    std_df <- std_df[seq_len(top_n), , drop = FALSE]
  }

  .new_enrich_std(std_df)
}


# ------------------------------------------------------------------------------
# print / summary 方法
# ------------------------------------------------------------------------------

#' @export
print.EnrichStd <- function(x, ...) {
  cat("== EnrichStd 富集分析标准格式 ==\n")
  cat("Terms:", nrow(x), "\n")
  cat("列:", paste(colnames(x), collapse = ", "), "\n")
  if (nrow(x) > 0) {
    if ("p_adjust" %in% colnames(x) && any(!is.na(x$p_adjust))) {
      cat(
        "\np_adjust 范围:",
        round(min(x$p_adjust, na.rm = TRUE), 4),
        "~",
        round(max(x$p_adjust, na.rm = TRUE), 4),
        "\n"
      )
    }
    if ("gene_count" %in% colnames(x) && any(!is.na(x$gene_count))) {
      cat(
        "gene_count 范围:",
        min(x$gene_count, na.rm = TRUE),
        "~",
        max(x$gene_count, na.rm = TRUE),
        "\n"
      )
    }
    # 预览列：只选实际存在的列
    preview_cols <- intersect(
      c("term_name", "p_adjust", "gene_count"),
      colnames(x)
    )
    if (length(preview_cols) > 0) {
      cat("\n前", min(3, nrow(x)), "行预览:\n")
      preview <- x[seq_len(min(3, nrow(x))), preview_cols, drop = FALSE]
      class(preview) <- "data.frame"
      print(preview, row.names = FALSE)
    }
  }
  invisible(x)
}

#' @export
summary.EnrichStd <- function(object, ...) {
  cat("== EnrichStd 摘要 ==\n")
  cat("Terms 总数:", nrow(object), "\n")
  cat("p_adjust < 0.05:", sum(object$p_adjust < 0.05, na.rm = TRUE), "\n")
  cat("p_adjust < 0.01:", sum(object$p_adjust < 0.01, na.rm = TRUE), "\n")
  cat("gene_count 均值:", round(mean(object$gene_count, na.rm = TRUE), 1), "\n")
  if (!all(is.na(object$rich_factor))) {
    cat(
      "rich_factor 均值:",
      round(mean(object$rich_factor, na.rm = TRUE), 3),
      "\n"
    )
  }
  invisible(object)
}
