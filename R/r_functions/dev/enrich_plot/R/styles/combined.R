# ==============================================================================
# plot_enrich_combined(): 多分类富集分析综合可视化
# 创建: 2026-04-24
#
# 依赖包:
#   必需: ggplot2, dplyr
#   必需: gground   (geom_round_col, geom_round_rect)
#         install: remotes::install_github("gground 的仓库地址")
#   可选: ggprism   (theme_prism)
#         install: install.packages("ggprism")
# ==============================================================================

# 依赖
# - core/std_format.R: as_enrich_std(), prep_enrich_combined()


# ------------------------------------------------------------------------------
# 内部工具：截断基因列表字符串
# ------------------------------------------------------------------------------

library(ggplot2)

.truncate_gene_ids <- function(gene_str, sep = "/", max_genes = 5) {
  genes <- strsplit(gene_str, sep, fixed = TRUE)[[1]]
  if (length(genes) <= max_genes) {
    return(paste(genes, collapse = ", "))
  }
  paste0(paste(genes[seq_len(max_genes)], collapse = ", "),
         " ... (+", length(genes) - max_genes, ")")
}


# ------------------------------------------------------------------------------
# plot_enrich_combined()
# ------------------------------------------------------------------------------

#' 多分类富集分析综合可视化
#'
#' 将多个富集分析结果（如 GO 三个本体 + KEGG）整合为单一图形：
#' - 左侧圆角色块标注分类，气泡显示基因数量
#' - 右侧圆角条形图展示 -log10(p.adjust)
#' - 条内/条下显示通路名称和基因列表
#'
#' @param data 输入数据，支持两种形式：
#'   - **named list**：元素为 enrichResult / EnrichStd / data.frame，
#'     名称作为分类标签。
#'     示例：`list(BP = go_bp, CC = go_cc, MF = go_mf, KEGG = kegg_res)`
#'   - **data.frame**：已合并数据，需配合 `ontology_col` 指定分类列
#' @param ontology_col 分类列名，仅 data.frame 输入时使用
#' @param col_map 列名映射，见 as_enrich_std()
#' @param top_n 每个分类取前 N 个 term，默认 5
#' @param category_order 分类排列顺序（y 轴从下到上），
#'   NULL 时使用输入顺序的逆序
#' @param category_colors 分类配色，named vector。
#'   NULL 时自动生成
#' @param show_gene_ids 是否在条形下方显示基因列表，默认 TRUE
#' @param max_genes_shown 每个条形显示的最多基因数，默认 5
#' @param label_size 通路名称字号，默认 4.5
#' @param gene_label_size 基因列表字号，默认 3
#' @param bar_width 条形宽度，默认 0.6
#' @param bar_alpha 条形透明度，默认 0.8
#' @param bubble_width 左侧气泡的 x 位置（相对单位），默认 0.8
#' @param rect_width 左侧分类色块的宽度（相对单位），默认 0.5
#' @param title 图标题，默认 NULL
#' @param use_round_geom 是否使用圆角 geom（需要 gground 包），默认 TRUE，
#'   若 gground 未安装会自动降级为普通 geom
#' @param use_prism_theme 是否使用 theme_prism()（需要 ggprism），默认 TRUE，
#'   若未安装会降级为 theme_bw()
#'
#' @return ggplot 对象
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler); library(org.Hs.eg.db)
#' source("core/std_format.R"); source("styles/combined.R")
#'
#' # 方式1: named list（最常用）
#' go  <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "ALL")
#' keg <- enrichKEGG(gene = gene_ids, organism = "hsa")
#'
#' plot_enrich_combined(
#'   list(BP = subset(go, ONTOLOGY=="BP"),
#'        CC = subset(go, ONTOLOGY=="CC"),
#'        MF = subset(go, ONTOLOGY=="MF"),
#'        KEGG = keg),
#'   top_n = 5
#' )
#'
#' # 方式2: 已合并 data.frame，指定分类列
#' combined_df$my_category <- ...
#' plot_enrich_combined(combined_df, ontology_col = "my_category", top_n = 4)
#' }
#'
#' @export
plot_enrich_combined <- function(
  data,
  ontology_col      = NULL,
  col_map           = NULL,
  top_n             = 5,
  category_order    = NULL,
  category_colors   = NULL,
  show_gene_ids     = TRUE,
  max_genes_shown   = 5,
  label_size        = 4.5,
  gene_label_size   = 3,
  bar_width         = 0.6,
  bar_alpha         = 0.8,
  bubble_width      = 0.20,
  rect_width        = 0.13,
  title             = NULL,
  x_var             = "p_adjust",
  base_size         = 11,
  use_round_geom    = TRUE,
  use_prism_theme   = TRUE
) {
  # 1. 检查可选依赖 ----
  has_gground <- requireNamespace("gground", quietly = TRUE)
  has_ggprism <- requireNamespace("ggprism", quietly = TRUE)

  if (use_round_geom && !has_gground) {
    message("gground 包未安装，已降级使用普通 geom_col / geom_rect。",
            "\n安装方法请咨询包的维护者或使用 remotes::install_github()")
    use_round_geom <- FALSE
  }
  if (use_prism_theme && !has_ggprism) {
    message("ggprism 包未安装，已降级使用 theme_bw()。",
            "\n安装: install.packages('ggprism')")
    use_prism_theme <- FALSE
  }

  # 2. 数据准备 ----
  plot_data <- prep_enrich_combined(
    data           = data,
    ontology_col   = ontology_col,
    col_map        = col_map,
    top_n          = top_n,
    category_order = category_order
  )

  cats   <- levels(plot_data$category)
  n_cats <- length(cats)

  # 3. 分类配色 ----
  if (is.null(category_colors)) {
    default_pal <- c(
      '#4393C3', '#D6604D', '#74C476', '#FD8D3C',
      '#9E9AC8', '#41AB5D', '#FC8D59', '#80B1D3'
    )
    category_colors <- setNames(
      default_pal[seq_len(n_cats)],
      cats
    )
  } else {
    # 补全缺失的分类颜色
    missing_cats <- setdiff(cats, names(category_colors))
    if (length(missing_cats) > 0) {
      extra_cols <- grDevices::rainbow(length(missing_cats))
      category_colors[missing_cats] <- extra_cols
    }
  }

  # 3.5 计算 x 轴值 ----
  plot_data$.x_val <- .compute_x_col(plot_data, x_var)
  x_label <- switch(x_var,
    gene_ratio = "GeneRatio",
    gene_count = "Gene Count",
    "-log10(p.adjust)"
  )

  # 4. 布局参数 ----
  # 所有左侧位置按 xaxis_max 等比缩放，确保 x 轴量级无论是 p_adjust 还是 gene_ratio 都协调
  xaxis_max   <- max(plot_data$.x_val, na.rm = TRUE) * 1.1
  left_bubble <- -bubble_width  * xaxis_max          # 基因数量气泡 x 位置
  left_rect_r <- -(bubble_width + rect_width * 0.4) * xaxis_max   # 分类色块右边界
  left_rect_l <- -(bubble_width + rect_width * 1.4) * xaxis_max   # 分类色块左边界
  label_x0    <- xaxis_max * 0.02    # 条内标签起始 x（替代硬编码 0.05）
  gene_x0     <- xaxis_max * 0.03    # 基因列表起始 x（替代硬编码 0.1）

  # 5. 处理基因列表显示 ----
  if (show_gene_ids) {
    plot_data$gene_label <- sapply(
      plot_data$gene_ids,
      .truncate_gene_ids,
      max_genes = max_genes_shown
    )
  }

  # 6. 构建左侧分类色块数据 ----
  rect_data <- do.call(rbind, lapply(cats, function(cat) {
    rows <- which(plot_data$category == cat)
    data.frame(
      category = cat,
      ymin     = min(rows) - 0.4,
      ymax     = max(rows) + 0.4,
      xmin     = left_rect_l,
      xmax     = left_rect_r,
      stringsAsFactors = FALSE
    )
  }))
  rect_data$category <- factor(rect_data$category, levels = cats)

  # 7. 选择 geom 函数 ----
  if (use_round_geom) {
    geom_bar_fn  <- gground::geom_round_col
    geom_rect_fn <- gground::geom_round_rect
  } else {
    geom_bar_fn  <- geom_col
    geom_rect_fn <- function(..., radius = NULL) geom_rect(...)
  }

  # 8. 绘图 ----
  p <- ggplot(plot_data,
              aes(x = .x_val, y = index, fill = category)) +

    # 主体条形
    geom_bar_fn(
      aes(y = term_name),
      width = bar_width,
      alpha = bar_alpha
    ) +

    # 通路名称标签（条内左对齐）
    geom_text(
      aes(x = label_x0, label = as.character(term_name)),
      hjust = 0,
      size  = label_size,
      color = "black"
    ) +

    # 基因列表标签（条底边下方，用显式 y 偏移避免被相邻条遮挡）
    {
      if (show_gene_ids) {
        geom_text(
          aes(x     = gene_x0,
              y     = index - bar_width / 2 - 0.05,
              label = gene_label,
              colour = category),
          hjust       = 0,
          vjust       = 1,        # 文字顶端贴在计算出的 y 位置
          size        = gene_label_size,
          fontface    = "italic",
          show.legend = FALSE
        )
      }
    } +

    # 左侧基因数量气泡
    geom_point(
      aes(x = left_bubble, size = gene_count),
      shape = 21
    ) +
    geom_text(
      aes(x = left_bubble, label = gene_count),
      size = gene_label_size + 0.5
    ) +
    scale_size_continuous(name = "Gene Count", range = c(5, 16)) +

    # 左侧分类色块
    geom_rect_fn(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
      data         = rect_data,
      inherit.aes  = FALSE
    ) +
    geom_text(
      aes(
        x     = (xmin + xmax) / 2,
        y     = (ymin + ymax) / 2,
        label = category
      ),
      data        = rect_data,
      inherit.aes = FALSE,
      size        = label_size - 0.5,
      color       = "white",
      fontface    = "bold"
    ) +

    # 手动 x 轴线（用 annotate 避免数据行数不匹配的警告）
    annotate(
      "segment",
      x = 0, y = 0.5, xend = xaxis_max, yend = 0.5,
      linewidth = 1.5, color = "black"
    ) +

    # 坐标轴和图例
    scale_fill_manual(name = "Category", values = category_colors) +
    scale_colour_manual(values = category_colors, guide = "none") +
    scale_x_continuous(
      breaks = pretty(c(0, xaxis_max)),
      expand = expansion(c(0, 0))
    ) +
    labs(x = x_label, y = NULL, title = title)

  # 9. 主题 ----
  base_theme <- if (use_prism_theme) {
    ggprism::theme_prism(base_size = base_size)
  } else {
    theme_bw(base_size = base_size)
  }

  p <- p + base_theme +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank(),
      legend.title = element_text(),
      plot.title   = element_text(size = base_size + 3, face = "bold", hjust = 0.5)
    )

  p
}


# 注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot(
    type        = "combined",
    fn          = plot_enrich_combined,
    description = "多分类综合图（左侧分类标签+气泡，右侧圆角条形+基因列表）",
    overwrite   = TRUE
  )
}
