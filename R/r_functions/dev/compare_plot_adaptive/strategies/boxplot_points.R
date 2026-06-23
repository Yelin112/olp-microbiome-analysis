#' 箱线图/小提琴图 + 散点策略（适用于 10 ≤ n < 20）
#'
#' 箱线图（或小提琴图）叠加全部数据散点，不重复显示离群点
#'
#' @param data 数据框
#' @param group_col 分组列名（字符串）
#' @param value_col 数值列名（字符串）
#' @param colors 颜色向量，长度须 >= 分组数
#' @param use_violin 是否使用小提琴图替代箱线图，默认 FALSE
#' @param box_width 箱宽，默认 0.6
#' @param box_alpha 箱体透明度，默认 0.7
#' @param point_size 散点大小，默认 1.5
#' @param point_alpha 散点透明度，默认 0.6
#' @param point_jitter 散点抖动幅度，默认 0.15
#'
#' @return ggplot 对象
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_jitter stat_summary
#' @importFrom ggplot2 scale_fill_manual theme_classic
#' @importFrom rlang sym
#'
#' @examples
#' data <- data.frame(group = rep(c("A","B","C","D"), each=15),
#'                    value = rnorm(60, rep(c(35,38,25,30), each=15), 5))
#' colors <- c("#4DBBD5", "#E64B35", "#00A087", "#F39B7F")
#' plot_boxplot_points(data, "group", "value", colors)
#'
#' @export
plot_boxplot_points <- function(data,
                                 group_col,
                                 value_col,
                                 colors,
                                 use_violin   = FALSE,
                                 box_width    = 0.6,
                                 box_alpha    = 0.7,
                                 point_size   = 1.5,
                                 point_alpha  = 0.6,
                                 point_jitter = 0.15) {

  p <- ggplot2::ggplot(data,
    ggplot2::aes(
      x    = !!rlang::sym(group_col),
      y    = !!rlang::sym(value_col),
      fill = !!rlang::sym(group_col)
    )
  )

  if (use_violin) {
    # 小提琴图模式
    p <- p +
      ggplot2::geom_violin(
        scale     = "width",
        trim      = FALSE,
        width     = box_width,
        alpha     = box_alpha,
        color     = "black",
        linewidth = 0.5
      ) +
      # 内置箱线图（细）
      ggplot2::geom_boxplot(
        width          = box_width * 0.25,
        fill           = "white",
        outlier.shape  = NA,
        color          = "black",
        linewidth      = 0.4,
        alpha          = 0.9
      )
  } else {
    # 箱线图模式：隐藏离群点（由散点层统一展示）
    p <- p +
      ggplot2::geom_boxplot(
        width          = box_width,
        alpha          = box_alpha,
        outlier.shape  = NA,
        color          = "black",
        linewidth      = 0.5
      )
  }

  p <- p +
    # 散点层：展示所有数据点
    ggplot2::geom_jitter(
      shape  = 21,
      size   = point_size,
      stroke = 0.3,
      fill   = "white",
      color  = "black",
      alpha  = point_alpha,
      width  = point_jitter,
      height = 0
    ) +
    ggplot2::scale_fill_manual(values = colors, name = group_col) +
    ggplot2::theme_classic()

  return(p)
}
