#' 纯箱线图策略（适用于 n ≥ 20）
#'
#' 大样本量时使用箱线图或小提琴图，可选叠加少量代表性散点
#'
#' @param data 数据框
#' @param group_col 分组列名（字符串）
#' @param value_col 数值列名（字符串）
#' @param colors 颜色向量，长度须 >= 分组数
#' @param use_violin 是否使用小提琴图，默认 FALSE
#' @param show_points 是否叠加散点（大数据量时建议关闭），默认 FALSE
#' @param box_width 箱宽，默认 0.7
#' @param box_alpha 箱体透明度，默认 0.8
#' @param point_size 散点大小（show_points=TRUE时），默认 1
#' @param point_alpha 散点透明度（show_points=TRUE时），默认 0.4
#'
#' @return ggplot 对象
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_jitter stat_summary
#' @importFrom ggplot2 scale_fill_manual theme_classic
#' @importFrom rlang sym
#'
#' @examples
#' data <- data.frame(group = rep(c("A","B","C","D"), each=30),
#'                    value = rnorm(120, rep(c(35,38,25,30), each=30), 5))
#' colors <- c("#4DBBD5", "#E64B35", "#00A087", "#F39B7F")
#' plot_pure_boxplot(data, "group", "value", colors)
#'
#' @export
plot_pure_boxplot <- function(data,
                               group_col,
                               value_col,
                               colors,
                               use_violin   = FALSE,
                               show_points  = FALSE,
                               box_width    = 0.7,
                               box_alpha    = 0.8,
                               point_size   = 1,
                               point_alpha  = 0.4) {

  p <- ggplot2::ggplot(data,
    ggplot2::aes(
      x    = !!rlang::sym(group_col),
      y    = !!rlang::sym(value_col),
      fill = !!rlang::sym(group_col)
    )
  )

  if (use_violin) {
    p <- p +
      ggplot2::geom_violin(
        scale     = "width",
        trim      = FALSE,
        width     = box_width,
        alpha     = box_alpha,
        color     = "black",
        linewidth = 0.5
      ) +
      ggplot2::geom_boxplot(
        width          = box_width * 0.2,
        fill           = "white",
        outlier.shape  = NA,
        color          = "black",
        linewidth      = 0.4,
        alpha          = 0.9
      )
  } else {
    p <- p +
      ggplot2::geom_boxplot(
        width          = box_width,
        alpha          = box_alpha,
        outlier.shape  = if (show_points) NA else 19,
        outlier.size   = 1.5,
        outlier.alpha  = 0.5,
        color          = "black",
        linewidth      = 0.5
      )
  }

  # 叠加均值点
  p <- p +
    ggplot2::stat_summary(
      fun       = mean,
      geom      = "point",
      shape     = 21,
      size      = 2.5,
      fill      = "white",
      color     = "black",
      stroke    = 0.8,
      show.legend = FALSE
    )

  # 可选：叠加散点（大数据量慎用）
  if (show_points) {
    p <- p +
      ggplot2::geom_jitter(
        shape  = 16,
        size   = point_size,
        alpha  = point_alpha,
        color  = "black",
        width  = 0.2,
        height = 0
      )
  }

  p <- p +
    ggplot2::scale_fill_manual(values = colors, name = group_col) +
    ggplot2::theme_classic()

  return(p)
}
