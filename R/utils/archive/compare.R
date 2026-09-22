library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

#' 组间比较可视化函数 (compare_plot)
#'
#' @description
#' 一个高度集成的 ggplot2 包装器，用于快速绘制高质量的组间比较图。
#' 支持箱线图、小提琴图、柱状图等，内置统计检验、斑马纹背景、趋势线及智能配色。
#'
#' @param data 数据框 (data.frame)。
#' @param value.var Y轴数值列名 (字符串)。
#' @param group.by X轴分组列名 (字符串)。
#' @param fill.by 填充颜色列名 (字符串)。默认与 group.by 相同。若指定不同列，则开启“组内簇状比较”模式。
#' @param split.by 分面列名 (字符串)。用于 facet_wrap 分面。
#' @param add_bg 是否添加斑马纹背景 (逻辑值)。默认 FALSE。
#' @param bg_color 背景颜色 (字符串)。默认为 "#0000000D" (5%透明度的黑色)，形成极淡的阴影。
#' @param plot_type 主图类型: "violin", "box", "bar", "dot"。
#' @param add_box 是否叠加箱线图 (逻辑值)。
#' @param add_point 是否叠加散点 (逻辑值)。
#' @param add_trend 是否叠加趋势线 (逻辑值)。连线为实线。
#' @param add_stat 统计方法: "none", "mean", "median", "t.test", "wilcox.test", "anova", "kruskal.test"。
#' @param comparisons 指定比较组 (List)。例如 list(c("A","B"), c("A","C"))。
#' @param palette 配色方案。可以是预设名 ("NPG", "AAAS", "Paired"...) 或自定义颜色向量。
#' @param theme_use 使用的主题函数。例如 theme_classic, theme_bw, theme_minimal。
#' @param title,xlab,ylab 标题与轴标签。
#' @param ... 传递给 theme_use 的其他参数 (如 base_size)。
#' @export
compare_plot <- function(
  data,
  value.var,
  group.by,
  fill.by = NULL,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "#0000000D",
  plot_type = c("violin", "box", "bar", "dot"),
  add_box = FALSE,
  add_point = FALSE,
  add_trend = FALSE,
  add_stat = "none",
  comparisons = NULL,
  palette = "NPG",
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  theme_use = theme_classic,
  ...
) {
  # --- 0. 参数校验 ---
  plot_type <- match.arg(plot_type)
  if (is.null(xlab)) {
    xlab <- group.by
  }
  if (is.null(ylab)) {
    ylab <- value.var
  }

  fill_var <- if (is.null(fill.by)) group.by else fill.by

  cols_needed <- c(value.var, group.by, fill_var)
  if (!all(cols_needed %in% colnames(data))) {
    stop("指定的列名在数据框中不存在，请检查拼写。")
  }

  for (col in unique(c(group.by, fill_var))) {
    if (!is.factor(data[[col]])) {
      data[[col]] <- factor(data[[col]], levels = unique(data[[col]]))
    }
  }

  # --- 1. 颜色解析 ---
  n_groups <- nlevels(data[[fill_var]])
  sci_palettes <- list(
    "NPG" = c(
      "#E64B35",
      "#4DBBD5",
      "#00A087",
      "#3C5488",
      "#F39B7F",
      "#8491B4",
      "#91D1C2",
      "#DC0000",
      "#7E6148"
    ),
    "AAAS" = c(
      "#3B4992",
      "#EE0000",
      "#008B45",
      "#631879",
      "#008280",
      "#BB0021",
      "#5F559B",
      "#A20056",
      "#808180"
    ),
    "NEJM" = c(
      "#BC3C29",
      "#0072B5",
      "#E18727",
      "#20854E",
      "#7876B1",
      "#6F99AD",
      "#FFDC91",
      "#EE4C97"
    ),
    "Lancet" = c(
      "#00468B",
      "#ED0000",
      "#42B540",
      "#0099B4",
      "#925E9F",
      "#FDAF91",
      "#AD002A",
      "#ADB6B6"
    ),
    "JCO" = c(
      "#0073C2",
      "#EFC000",
      "#868686",
      "#CD534C",
      "#7AA6DC",
      "#003C67",
      "#8F7700",
      "#3B3B3B"
    )
  )

  if (length(palette) == 1 && is.character(palette)) {
    if (palette %in% names(sci_palettes)) {
      base_cols <- sci_palettes[[palette]]
    } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      base_cols <- RColorBrewer::brewer.pal(min(n_groups, max_n), palette)
    } else {
      warning("Palette not found. Using NPG.")
      base_cols <- sci_palettes[["NPG"]]
    }
  } else {
    base_cols <- palette
  }

  if (n_groups > length(base_cols)) {
    final_colors <- colorRampPalette(base_cols)(n_groups)
  } else {
    final_colors <- base_cols[1:n_groups]
  }

  # --- 2. 绘图初始化 ---
  p <- ggplot(
    data,
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[fill_var]])
  )

  is_grouped <- (fill_var != group.by)
  dodge_width <- 0.8

  # --- 3. 背景层 ---
  if (isTRUE(add_bg)) {
    n_x_groups <- nlevels(data[[group.by]])
    bg_df <- data.frame(x_pos = seq(1, n_x_groups, 2)) %>%
      mutate(xmin = x_pos - 0.5, xmax = x_pos + 0.5)

    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = bg_color,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  # --- 4. 主数据图层 ---
  pos <- if (is_grouped) {
    position_dodge(width = dodge_width)
  } else {
    position_dodge(width = 0)
  }

  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3,
        position = pos
      )
  } else if (plot_type == "box") {
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black",
        position = pos
      ) +
      # [修复] 显式指定 group = fill_var，确保白点能跟着箱子一起错位
      stat_summary(
        aes(group = .data[[fill_var]]),
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2.5,
        fill = "white",
        color = "black",
        position = pos,
        show.legend = FALSE
      )
  } else if (plot_type == "bar") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black",
        position = pos
      ) +
      stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.2,
        position = pos
      )
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black",
        position = pos
      ) +
      stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.2,
        position = pos
      )
  }

  # --- 5. 叠加层 ---
  # (A) 趋势线
  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(
          group = if (is_grouped) .data[[fill_var]] else 1,
          color = if (is_grouped) .data[[fill_var]] else NULL
        ),
        linewidth = 1,
        linetype = "solid",
        position = pos,
        show.legend = FALSE
      )
    if (!is_grouped) p <- p + aes(color = NULL)
  }

  # (B) 叠加 Box
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8,
        position = pos
      ) +
      # [修复] 叠加 Box 的白点也需要指定 group
      stat_summary(
        aes(group = .data[[fill_var]]),
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2,
        fill = "white",
        color = "black",
        position = pos,
        show.legend = FALSE
      )
  }

  # (C) 散点
  if (isTRUE(add_point)) {
    if (is_grouped) {
      p <- p +
        geom_point(
          position = position_jitterdodge(
            jitter.width = 0.2,
            dodge.width = dodge_width
          ),
          size = 1.5,
          alpha = 0.6,
          color = "grey30",
          shape = 16
        )
    } else {
      p <- p +
        geom_jitter(
          width = 0.2,
          size = 1.5,
          alpha = 0.6,
          color = "grey30",
          shape = 16
        )
    }
  }

  # --- 6. 统计层 ---
  if (add_stat %in% c("mean", "median")) {
    # [修复] 用户手动开启的统计点也需要指定 group
    p <- p +
      stat_summary(
        aes(group = .data[[fill_var]]),
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black",
        position = pos,
        show.legend = FALSE
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    if (is_grouped) {
      p <- p +
        stat_compare_means(
          aes(group = .data[[fill_var]]),
          method = add_stat,
          label = "p.signif"
        )
    } else {
      label_format <- if (
        is.null(comparisons) && nlevels(data[[group.by]]) > 2
      ) {
        "p.format"
      } else {
        "p.signif"
      }
      p <- p +
        stat_compare_means(
          comparisons = comparisons,
          method = add_stat,
          label = label_format
        )
    }
  }

  # --- 7. 收尾 ---
  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  if (is.function(theme_use)) {
    p <- p + theme_use(...)
  } else {
    p <- p + theme_classic(...)
  }

  if (add_bg) {
    p <- p + theme(panel.grid = element_blank())
  }

  p <- p + scale_fill_manual(values = final_colors)
  p <- p + labs(title = title, x = xlab, y = ylab, fill = fill_var)

  if (is_grouped && add_trend) {
    p <- p + scale_color_manual(values = rep("black", n_groups), guide = "none")
  }

  return(p)
}
