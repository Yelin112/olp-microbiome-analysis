#' 纯散点图策略（适用于 n < 5）
#'
#' 每个数据点全部显示，叠加均值 ± 标准误差线（crossbar 风格）
#'
#' @param data 数据框
#' @param group_col 分组列名（字符串）
#' @param value_col 数值列名（字符串）
#' @param colors 颜色向量，长度须 >= 分组数
#' @param point_size 散点大小，默认 3
#' @param point_stroke 散点边框粗细，默认 0.5
#' @param errorbar_size 误差线粗细，默认 1.2
#' @param errorbar_width 误差线帽宽，默认 0.3
#'
#' @return ggplot 对象
#'
#' @importFrom ggplot2 ggplot aes geom_jitter stat_summary
#' @importFrom ggplot2 scale_fill_manual scale_color_manual theme_classic
#' @importFrom rlang sym
#'
#' @examples
#' data <- data.frame(group = rep(c("A","B","C"), each=3), value = rnorm(9, c(3,5,4), 0.5))
#' colors <- c("#E64B35", "#4DBBD5", "#00A087")
#' plot_pure_scatter(data, "group", "value", colors)
#'
#' @export
plot_pure_scatter <- function(data,
                               group_col,
                               value_col,
                               colors,
                               point_size    = 3,
                               point_stroke  = 0.5,
                               errorbar_size = 1.2,
                               errorbar_width = 0.3) {

  ggplot2::ggplot(data,
    ggplot2::aes(
      x    = !!rlang::sym(group_col),
      y    = !!rlang::sym(value_col),
      fill = !!rlang::sym(group_col)
    )
  ) +
    # 散点层：所有数据点，轻微抖动避免重叠
    ggplot2::geom_jitter(
      shape  = 21,
      size   = point_size,
      stroke = point_stroke,
      color  = "black",
      width  = 0.15,
      height = 0
    ) +
    # 均值线（crossbar）
    ggplot2::stat_summary(
      fun      = mean,
      geom     = "crossbar",
      width    = errorbar_width,
      linewidth = errorbar_size,
      color    = "black",
      fatten   = 1
    ) +
    # 标准误差线
    ggplot2::stat_summary(
      fun.data = mean_se,
      geom     = "errorbar",
      width    = errorbar_width * 0.6,
      linewidth = errorbar_size * 0.7,
      color    = "black"
    ) +
    ggplot2::scale_fill_manual(values = colors, name = group_col) +
    ggplot2::theme_classic()
}
