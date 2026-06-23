#' 柱状图 + 散点策略（适用于 5 ≤ n < 10）
#'
#' 均值柱状图叠加标准误差线，全部数据点以黑色散点叠加显示
#'
#' @param data 数据框
#' @param group_col 分组列名（字符串）
#' @param value_col 数值列名（字符串）
#' @param colors 颜色向量，长度须 >= 分组数
#' @param bar_width 柱宽，默认 0.7
#' @param bar_alpha 柱子透明度，默认 0.85
#' @param point_size 散点大小，默认 2
#' @param point_alpha 散点透明度，默认 0.8
#' @param point_jitter 散点抖动幅度，默认 0.1
#'
#' @return ggplot 对象
#'
#' @importFrom ggplot2 ggplot aes stat_summary geom_jitter
#' @importFrom ggplot2 scale_fill_manual theme_classic
#' @importFrom rlang sym
#'
#' @examples
#' data <- data.frame(group = rep(c("A","B","C","D"), each=7),
#'                    value = rnorm(28, rep(c(35,38,25,30), each=7), 4))
#' colors <- c("#4DBBD5", "#E64B35", "#00A087", "#F39B7F")
#' plot_bar_points(data, "group", "value", colors)
#'
#' @export
plot_bar_points <- function(data,
                             group_col,
                             value_col,
                             colors,
                             bar_width    = 0.7,
                             bar_alpha    = 0.85,
                             point_size   = 2,
                             point_alpha  = 0.8,
                             point_jitter = 0.1) {

  ggplot2::ggplot(data,
    ggplot2::aes(
      x    = !!rlang::sym(group_col),
      y    = !!rlang::sym(value_col),
      fill = !!rlang::sym(group_col)
    )
  ) +
    # 均值柱状图
    ggplot2::stat_summary(
      fun      = mean,
      geom     = "bar",
      width    = bar_width,
      alpha    = bar_alpha,
      color    = "black",
      linewidth = 0.5
    ) +
    # 标准误差线
    ggplot2::stat_summary(
      fun.data  = mean_se,
      geom      = "errorbar",
      width     = bar_width * 0.35,
      linewidth = 0.7,
      color     = "black"
    ) +
    # 散点层：黑色填充，覆盖在柱状图上方
    ggplot2::geom_jitter(
      shape  = 21,
      size   = point_size,
      stroke = 0.4,
      fill   = "black",
      color  = "white",
      alpha  = point_alpha,
      width  = point_jitter,
      height = 0
    ) +
    ggplot2::scale_fill_manual(values = colors, name = group_col) +
    ggplot2::theme_classic()
}
