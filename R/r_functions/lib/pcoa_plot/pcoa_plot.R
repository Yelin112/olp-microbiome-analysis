# ============================================================================
# pcoa_plot() — PCoA / 二维散点图（含边际图）
# 开发日期：2026-04-06
# 依赖：ggplot2（必需）、utils/palette_system.R、ggExtra（density 边际）、aplot（boxplot 边际）
# ============================================================================

library(ggplot2)

# ---- 加载工具模块（配色、%||%）----
if (!exists("%||%") || !exists("get_colors")) {
  # 推荐：脚本先 source R/init.R 再 load_utils()。此处为兜底，依赖 OLP_ROOT 环境变量
  .init <- file.path(Sys.getenv("OLP_ROOT"), "R", "init.R")
  if (!file.exists(.init)) {
    stop("未找到 utils 函数（%||% / get_colors）。请先 source R/init.R 并调用 load_utils()，",
         "或设置 OLP_ROOT 环境变量。")
  }
  source(.init)
  load_utils()
  rm(.init)
}

# ── 内部工具 ─────────────────────────────────────────────────────────────────

#' 获取分组颜色：命名向量原样返回（由主函数按组名校验），其余交给 palette_system
#' @noRd
.pcoa_get_colors <- function(palette, n) {
  if (length(palette) > 1 && !is.null(names(palette))) {
    return(palette)
  }
  get_colors(palette, n = n)
}

# ── 主函数 ───────────────────────────────────────────────────────────────────

#' PCoA / 二维散点图（含边际图）
#'
#' 绘制二维坐标散点图（PCoA、PCA、UMAP 等），支持：
#' - 分组着色 + 置信椭圆
#' - 统计结果注释文字
#' - 密度曲线边际图（类似 ggExtra::ggMarginal）
#' - 箱线图边际图（顶部 + 右侧，类似发表图）
#'
#' @param data         data.frame，包含坐标列和分组列
#' @param x            字符串，x 轴坐标列名（如 "PCoA1"）
#' @param y            字符串，y 轴坐标列名（如 "PCoA2"）
#' @param group.by     字符串，分组列名
#' @param marginal     边际图类型：
#'   \itemize{
#'     \item \code{"density"}：顶部+右侧密度曲线（需要 ggExtra 包）
#'     \item \code{"boxplot"}：顶部+右侧箱线图（需要 aplot 包）
#'     \item \code{"none"}：不添加边际图，返回纯 ggplot 对象
#'   }
#' @param add_ellipse  是否绘制置信椭圆，默认 TRUE
#' @param ellipse_level 椭圆置信水平，默认 0.95
#' @param ellipse_type  椭圆统计类型："norm"（默认）、"t"、"euclid"
#' @param ellipse_alpha 椭圆填充透明度，默认 0.15
#' @param stat_text    统计注释文字（字符向量，每个元素一行）。例如：
#'   \code{c("Adonis R² = 0.33", "P = 0.001")}。NULL 则不显示。
#' @param stat_x       注释 x 坐标。Inf = 右边缘（默认），-Inf = 左边缘，
#'   或传入数据坐标值
#' @param stat_y       注释 y 坐标。Inf = 上边缘（默认），-Inf = 下边缘，
#'   或传入数据坐标值
#' @param stat_hjust   注释水平对齐：1 = 右对齐（默认），0 = 左对齐，0.5 = 居中
#' @param stat_vjust   注释垂直对齐：1.2（默认），值越大文字越向下偏移
#' @param stat_size    注释字体大小，默认 3.5
#' @param stat_italic  注释是否斜体，默认 TRUE
#' @param sig_letters  显著性字母（仅 boxplot 边际图有效），命名字符向量。
#'   格式：\code{c(groupA = "a", groupB = "b", groupC = "b")}，
#'   names 需与分组水平一致
#' @param point_size   散点大小，默认 2.5
#' @param point_alpha  散点透明度，默认 0.85
#' @param point_shape  散点形状，默认 16（实心圆）；21-25 为有描边的填充形状
#' @param palette      配色方案，交给 palette_system.R 的 \code{get_colors()} 解析：
#'   预设名（"NPG" 默认、"AAAS"、"JAMA" 等）、RColorBrewer 色板名，
#'   或颜色向量（如 \code{c("#E64B35", "#4DBBD5")}）；
#'   命名向量（如 \code{c(Control = "#E64B35", Treatment = "#4DBBD5")}）按组名映射
#' @param marginal_size 边际图宽/高占主图的比例，默认 0.25（即 25\%）
#' @param xlab         x 轴标签，NULL 则使用列名
#' @param ylab         y 轴标签，NULL 则使用列名
#' @param title        图标题，默认 NULL
#' @param theme_use    ggplot2 主题函数，默认 theme_classic
#' @param ...          传递给 \code{geom_point()} 的其他参数
#'   （如 \code{stroke = 0.5}，仅对 shape 21-25 有效）
#'
#' @return
#' \itemize{
#'   \item \code{marginal = "none"}：ggplot 对象（可直接用 \code{+} 叠加图层）
#'   \item \code{marginal = "density"}：ggExtraPlot 对象（支持 ggsave）
#'   \item \code{marginal = "boxplot"}：aplot 对象（支持 ggsave）
#' }
#'
#' @examples
#' # 基本用法：密度边际图
#' set.seed(42)
#' df <- data.frame(
#'   PC1 = c(rnorm(30, -0.2), rnorm(30, 0.2)),
#'   PC2 = c(rnorm(30, 0.1), rnorm(30, -0.1)),
#'   group = rep(c("Control", "Treatment"), each = 30)
#' )
#' pcoa_plot(df, x = "PC1", y = "PC2", group.by = "group")
#'
#' # 箱线图边际图 + 统计注释 + 显著性字母
#' pcoa_plot(df, x = "PC1", y = "PC2", group.by = "group",
#'           marginal = "boxplot",
#'           stat_text = c("Adonis R² = 0.33; P = 0.001"),
#'           sig_letters = c(Control = "b", Treatment = "a"))
#'
#' @export
pcoa_plot <- function(
  data,
  x,
  y,
  group.by,
  marginal = c("density", "boxplot", "none"),
  add_ellipse = TRUE,
  ellipse_level = 0.95,
  ellipse_type = "norm",
  ellipse_alpha = 0.15,
  stat_text = NULL,
  stat_x = Inf,
  stat_y = Inf,
  stat_hjust = 1,
  stat_vjust = 1.5,
  stat_size = 3.5,
  stat_italic = TRUE,
  sig_letters = NULL,
  point_size = 2.5,
  point_alpha = 0.85,
  point_shape = 16,
  palette = "NPG",
  marginal_size = 0.25,
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  theme_use = theme_classic,
  ...
) {
  # ── 输入验证 ────────────────────────────────────────────────────────────────
  if (!is.data.frame(data)) {
    stop("`data` 必须是 data.frame")
  }
  for (col in c(x, y, group.by)) {
    if (!col %in% names(data)) {
      stop(sprintf("列 '%s' 不存在于 data 中", col))
    }
  }
  if (!is.numeric(data[[x]]) || !is.numeric(data[[y]])) {
    stop(sprintf("列 '%s' 和 '%s' 必须是数值型", x, y))
  }

  marginal <- match.arg(marginal)

  # ── 缺失值处理 ──────────────────────────────────────────────────────────────
  na_mask <- is.na(data[[x]]) | is.na(data[[y]]) | is.na(data[[group.by]])
  if (any(na_mask)) {
    message("移除 ", sum(na_mask), " 行含 NA 的数据")
    data <- data[!na_mask, , drop = FALSE]
  }
  if (nrow(data) == 0) {
    stop("移除 NA 后数据为空")
  }

  # ── 分组因子化 ──────────────────────────────────────────────────────────────
  data[[group.by]] <- factor(data[[group.by]])
  grp_levels <- levels(data[[group.by]])
  n_grp <- length(grp_levels)

  # ── 配色 ────────────────────────────────────────────────────────────────────
  colors <- .pcoa_get_colors(palette, n_grp)
  # 用户传入命名向量时保留其 names（按组名映射颜色）；否则按分组水平顺序命名
  if (is.null(names(colors))) {
    names(colors) <- grp_levels
  } else {
    missing_grp <- setdiff(grp_levels, names(colors))
    if (length(missing_grp) > 0) {
      warning(
        "palette 缺少以下分组的颜色：",
        paste(missing_grp, collapse = ", "),
        "，将使用 NA（灰色）填充"
      )
    }
  }

  # ── 构建主图 ────────────────────────────────────────────────────────────────
  p <- ggplot(
    data,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      color = .data[[group.by]],
      fill = .data[[group.by]]
    )
  ) +
    geom_point(
      size = point_size,
      alpha = point_alpha,
      shape = point_shape,
      ...
    ) +
    scale_color_manual(values = colors, name = group.by) +
    scale_fill_manual(
      values = colors,
      name = group.by,
      na.value = "transparent"
    ) +
    labs(x = xlab %||% x, y = ylab %||% y, title = title) +
    theme_use()

  # ── 置信椭圆 ────────────────────────────────────────────────────────────────
  if (add_ellipse) {
    p <- p +
      stat_ellipse(
        aes(fill = .data[[group.by]]),
        geom = "polygon",
        level = ellipse_level,
        type = ellipse_type,
        alpha = ellipse_alpha,
        color = NA,
        show.legend = FALSE
      ) +
      stat_ellipse(
        aes(color = .data[[group.by]]),
        level = ellipse_level,
        type = ellipse_type,
        linetype = "dashed",
        linewidth = 0.5,
        show.legend = FALSE
      )
  }

  # ── 统计注释 ────────────────────────────────────────────────────────────────
  if (!is.null(stat_text)) {
    label <- paste(stat_text, collapse = "\n")
    fontface <- if (stat_italic) "italic" else "plain"
    p <- p +
      annotate(
        "text",
        x = stat_x,
        y = stat_y,
        label = label,
        hjust = stat_hjust,
        vjust = stat_vjust,
        size = stat_size,
        fontface = fontface
      )
  }

  # ── 边际图 ──────────────────────────────────────────────────────────────────
  if (marginal == "none") {
    return(p)
  }

  # ── 密度边际图（ggExtra） ────────────────────────────────────────────────────
  if (marginal == "density") {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
      stop("density 边际图需要 ggExtra 包，请运行：install.packages('ggExtra')")
    }
    # ggMarginal 的 size 参数表示主图是边际图的倍数
    # marginal_size = 0.25 → 主图是边际图的 4 倍 → size = 4
    marg_size_int <- max(2L, round(1 / marginal_size))
    result <- ggExtra::ggMarginal(
      p,
      type = "density",
      groupFill = TRUE,
      groupColour = TRUE,
      alpha = 0.5,
      size = marg_size_int
    )
    return(result)
  }

  # ── 箱线图边际图（aplot） ────────────────────────────────────────────────────
  if (marginal == "boxplot") {
    if (!requireNamespace("aplot", quietly = TRUE)) {
      stop("boxplot 边际图需要 aplot 包，请运行：install.packages('aplot')")
    }

    # 共用箱线图样式
    box_theme <- list(
      theme_use(),
      theme(
        legend.position = "none",
        axis.title = element_blank()
      )
    )

    # 右侧面板：y 变量分组箱线图（y 轴与主图对齐）
    p_right <- ggplot(
      data,
      aes(
        x = .data[[group.by]],
        y = .data[[y]],
        fill = .data[[group.by]]
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.6,
        alpha = 0.7,
        color = "black",
        linewidth = 0.5
      ) +
      scale_fill_manual(values = colors) +
      box_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank()
      )

    # 显著性字母叠加到右侧面板
    if (!is.null(sig_letters)) {
      # 验证分组名
      missing_grp <- setdiff(names(sig_letters), grp_levels)
      if (length(missing_grp) > 0) {
        warning(
          "sig_letters 中的分组名 [",
          paste(missing_grp, collapse = ", "),
          "] 不存在于数据中，已忽略"
        )
        sig_letters <- sig_letters[names(sig_letters) %in% grp_levels]
      }
      if (length(sig_letters) > 0) {
        y_range <- range(data[[y]], na.rm = TRUE)
        y_offset <- diff(y_range) * 0.06
        # 以 75th 百分位 + offset 作为字母位置
        y_q75 <- tapply(
          data[[y]],
          data[[group.by]],
          quantile,
          probs = 0.75,
          na.rm = TRUE
        )
        letter_df <- data.frame(
          grp = names(sig_letters),
          letter = as.character(sig_letters),
          ypos = as.numeric(y_q75[names(sig_letters)]) + y_offset,
          stringsAsFactors = FALSE
        )
        names(letter_df)[1] <- group.by
        p_right <- p_right +
          geom_text(
            data = letter_df,
            aes(
              x = .data[[group.by]],
              y = .data[["ypos"]],
              label = .data[["letter"]]
            ),
            inherit.aes = FALSE,
            size = 3.5,
            fontface = "bold"
          )
      }
    }

    # 顶部面板：x 变量分组箱线图（coord_flip 后 x 轴与主图 x 轴对齐）
    p_top <- ggplot(
      data,
      aes(
        x = .data[[group.by]],
        y = .data[[x]],
        fill = .data[[group.by]]
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.6,
        alpha = 0.7,
        color = "black",
        linewidth = 0.5
      ) +
      scale_fill_manual(values = colors) +
      coord_flip() +
      box_theme +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )

    # 用 aplot 拼合：先右后上
    result <- p |>
      aplot::insert_right(p_right, width = marginal_size) |>
      aplot::insert_top(p_top, height = marginal_size)

    return(result)
  }
}
