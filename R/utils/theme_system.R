# ============================================================================
# theme_system.R — 统一主题系统（主题工厂 + 期刊预设）
# 版本: 1.0
# ============================================================================
#
# 用法:
#   source("R/utils/helpers.R")
#   source("R/utils/theme_system.R")
#
# 依赖: ggplot2, helpers.R
#
# 设计原则:
#   - 三层架构: 基石主题 → 场景变体 → 期刊预设
#   - 本文件不定义调色板，颜色由 palette_system.R 负责
#   - 所有主题函数返回 theme 对象，用 + theme_xxx() 方式使用

library(ggplot2)

# ============================================================================
# Layer 1: 基石主题
# ============================================================================

#' 基础主题
#'
#' @description
#' 所有主题的底层模板。基于 theme_minimal 重建，
#' 提供了统一的白底、黑字、可控边框和网格。
#'
#' @param base_size 基础字号，默认 12。
#' @param border 是否显示全黑边框 (TRUE) 或仅显示坐标轴线 (FALSE)。
#' @param grid 是否显示灰色虚线网格。
#' @param ... 传递给 theme() 的其他参数，用于微调。
#'
#' @return ggplot2 theme 对象
#' @export
#'
#' @examples
#' ggplot(iris, aes(Species, Sepal.Length)) +
#'   geom_boxplot() +
#'   theme_pub_base()
theme_pub_base <- function(base_size = 12, border = TRUE, grid = TRUE, ...) {
  scale <- base_size / 12

  t <- theme_minimal(base_size = base_size) +
    theme(
      # ── 文本 ──
      text = element_text(color = "black"),
      plot.title = element_text(size = 14 * scale, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11 * scale, color = "grey40"),
      plot.caption = element_text(
        size = 9 * scale,
        color = "grey60",
        hjust = 1
      ),
      axis.title = element_text(size = 13 * scale),
      axis.text = element_text(size = 12 * scale, color = "black"),
      legend.title = element_text(size = 12 * scale, face = "bold"),
      legend.text = element_text(size = 11 * scale),
      strip.text = element_text(size = 12 * scale, face = "bold"),

      # ── 背景 ──
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.key = element_rect(fill = "transparent", color = NA),

      # ── 刻度线 ──
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),

      # ── 间距 ──
      plot.margin = margin(10, 10, 10, 10)
    )

  # 边框 vs 轴线
  if (border) {
    t <- t +
      theme(
        panel.border = element_rect(
          fill = "transparent",
          color = "black",
          linewidth = 1.5
        ),
        axis.line = element_blank()
      )
  } else {
    t <- t +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1.2)
      )
  }

  # 网格控制
  if (grid) {
    t <- t +
      theme(
        panel.grid.major = element_line(
          color = "grey90",
          linetype = "dashed",
          linewidth = 0.4
        ),
        panel.grid.minor = element_blank()
      )
  } else {
    t <- t + theme(panel.grid = element_blank())
  }

  t + theme(...)
}

# ============================================================================
# Layer 2: 场景变体
# ============================================================================

#' 统计图主题
#'
#' @description
#' 适合柱状图、箱线图、散点图等带有明确坐标轴参考线的图形。
#' 虚线网格 + 右侧图例。
#'
#' @inheritParams theme_pub_base
#' @export
theme_pub_stat <- function(base_size = 12, ...) {
  theme_pub_base(base_size, border = TRUE, grid = TRUE) +
    theme(
      panel.grid.major = element_line(
        color = "grey85",
        linetype = "dashed",
        linewidth = 0.5
      ),
      legend.position = "right"
    ) +
    theme(...)
}

#' 演示/汇报主题
#'
#' @description
#' 适合 PPT/海报/论文展示。大字体、无网格、底部图例。
#' 可通过调整 base_size 快速适配不同观看距离。
#'
#' @inheritParams theme_pub_base
#' @export
theme_pub_present <- function(base_size = 18, ...) {
  theme_pub_base(base_size, border = FALSE, grid = FALSE) +
    theme(
      axis.line = element_line(color = "grey30", linewidth = 1.2),
      legend.position = "bottom",
      plot.margin = margin(15, 15, 15, 15)
    ) +
    theme(...)
}

#' 海报主题
#'
#' @description 演示主题的超大字体变体，适合学术海报。
#' @inheritParams theme_pub_base
#' @export
theme_pub_poster <- function(base_size = 24, ...) {
  theme_pub_present(base_size, ...)
}

# ============================================================================
# Layer 3: 期刊风格预设（主题 + 可选的默认调色板建议）
# ============================================================================

#' Nature 系列期刊风格
#'
#' @description
#' 黑边框 + 右侧图例，适合 Nature Communications / Scientific Reports 等。
#' 建议搭配 palette = "NPG" 使用。
#'
#' @inheritParams theme_pub_base
#' @export
theme_nature <- function(base_size = 12, ...) {
  theme_pub_base(base_size, border = TRUE, grid = FALSE) +
    theme(
      legend.position = "right",
      panel.border = element_rect(
        fill = "transparent",
        color = "black",
        linewidth = 0.8
      ),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      strip.background = element_rect(
        fill = "grey95",
        color = "black",
        linewidth = 0.5
      ),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    theme(...)
}

#' Cell 系列期刊风格
#'
#' @description
#' 极简风格: 无边框、浅灰网格线、底部图例。
#' 适合 Cell Reports / Molecular Cell 等。
#'
#' @inheritParams theme_pub_base
#' @export
theme_cell <- function(base_size = 12, ...) {
  theme_pub_base(base_size, border = FALSE, grid = FALSE) +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey40", linewidth = 0.8),
      axis.ticks = element_line(color = "grey40"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 11),
      plot.title = element_text(size = 14, face = "plain"),
      legend.position = "bottom"
    ) +
    theme(...)
}

#' JAMA 系列期刊风格
#'
#' @description
#' 无边框、粗轴线、底部图例。适合 JAMA Network Open 等。
#' 建议搭配 palette = "JAMA" 使用。
#'
#' @inheritParams theme_pub_base
#' @export
theme_jama <- function(base_size = 12, ...) {
  theme_pub_base(base_size, border = FALSE, grid = FALSE) +
    theme(
      axis.line = element_line(color = "grey30", linewidth = 1),
      axis.ticks = element_line(color = "grey30"),
      legend.position = "bottom"
    ) +
    theme(...)
}

#' Science 系列期刊风格
#'
#' @description
#' 紧凑型，适合 Science / Science Advances 等高信息密度期刊。
#' 薄边框 + 小边距 + 上方图例。
#'
#' @inheritParams theme_pub_base
#' @export
theme_science <- function(base_size = 10, ...) {
  theme_pub_base(base_size, border = TRUE, grid = FALSE) +
    theme(
      panel.border = element_rect(
        fill = "transparent",
        color = "black",
        linewidth = 0.6
      ),
      legend.position = "top",
      plot.margin = margin(5, 5, 5, 5),
      axis.ticks.length = unit(0.15, "cm")
    ) +
    theme(...)
}
