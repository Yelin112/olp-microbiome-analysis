# 绘图尺寸固定函数 ---------------------------------------------------------------

library(ggplot2)
library(grid)
library(gtable)

# ==============================================================================
# 1. 定义用户接口函数 (就像 theme() 一样使用)
# ==============================================================================
#' 通过 + 号固定绘图区大小
#' @param width 宽度 (默认单位 cm)
#' @param height 高度 (默认单位 cm)
#' @param margin 边距 (可选)
fix_panel <- function(width = 4, height = 4, margin = NULL) {
  # 创建一个特殊的列表对象，存储参数
  structure(
    list(width = width, height = height, margin = margin),
    class = "fix_panel_params"
  )
}

# ==============================================================================
# 2. 定义 + 号的加法逻辑 (S3 方法)
# ==============================================================================
#' @method ggplot_add fix_panel_params
#' @export
ggplot_add.fix_panel_params <- function(object, plot, object_name) {
  # 1. 将参数存入 plot 对象的一个新槽位中 (不影响原有结构)
  plot$fix_panel_data <- object

  # 2. 修改 plot 的类名 (Class)
  #    这是关键！我们给它打上 "fixed_panel_plot" 的标签，
  #    这样 R 在打印它时，就会调用我们在下面定义的特殊 print 方法。
  class(plot) <- c("fixed_panel_plot", class(plot))

  return(plot)
}

# ==============================================================================
# 3. 劫持 Print 方法 (渲染时生效)
# ==============================================================================
#' @method print fixed_panel_plot
#' @export
print.fixed_panel_plot <- function(x, ...) {
  # 提取参数
  params <- x$fix_panel_data
  w <- unit(params$width, "cm")
  h <- unit(params$height, "cm")

  # 1. 转换为 gtable (图形表格)
  g <- ggplotGrob(x)

  # 2. 查找 Panel 位置 (逻辑同之前)
  panels <- grep("panel", g$layout$name)
  if (length(panels) > 0) {
    panel_index_w <- unique(g$layout$l[panels])
    panel_index_h <- unique(g$layout$t[panels])

    # 3. 强制修改尺寸
    # 宽度
    if (length(panel_index_w) == 1) {
      g$widths[panel_index_w] <- w
    } else {
      g$widths[panel_index_w] <- rep(w, length(panel_index_w))
    }
    # 高度
    if (length(panel_index_h) == 1) {
      g$heights[panel_index_h] <- h
    } else {
      g$heights[panel_index_h] <- rep(h, length(panel_index_h))
    }
  }

  # 4. 绘制
  grid.newpage()
  grid.draw(g)

  # 返回 gtable 以便兼容 ggsave (ggsave 默认会打印并捕捉输出)
  return(invisible(x))
}


# 主题设置系统 -----------------------------------------------------------------

library(ggplot2)
library(scales) # 用于颜色插值

# ==============================================================================
# 1. Level 1: 基石主题 (Base Theme)
# ==============================================================================
#' 全局基础主题
#' @param base_size 基础字号，默认12
#' @param border 是否显示全黑框 (TRUE) 或仅显示坐标轴线 (FALSE)
#' @param grid 是否显示网格 (通常在基石层关闭，交给场景层决定)
#' @param ... 传递给 theme() 的其他参数
theme_my_base <- function(base_size = 12, border = TRUE, grid = FALSE, ...) {
  # 1. 字体缩放比例 (核心)
  scale <- base_size / 12

  # 2. 基础定义
  t <- ggplot2::theme_minimal(base_size = base_size) +
    theme(
      # 文本
      text = element_text(color = "black"),
      plot.title = element_text(size = 14 * scale, face = "bold", hjust = 0),
      axis.title = element_text(size = 13 * scale),
      axis.text = element_text(size = 12 * scale, color = "black"),
      legend.title = element_text(size = 12 * scale, face = "bold"),
      legend.text = element_text(size = 11 * scale),

      # 背景
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.key = element_rect(fill = "transparent", color = NA),

      # 刻度线（theme_minimal 默认清除，需显式恢复）
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm")
    )

  # 3. 边框逻辑
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

  # 4. 网格逻辑
  if (!grid) {
    t <- t + theme(panel.grid = element_blank())
  }

  # 5. 允许透传
  return(t + theme(...))
}

# ==============================================================================
# 2. Level 2: 场景变体 (Scenario Variants)
# ==============================================================================
#' 场景：SCP 统计图风格 (StatPlot)
#' 特点：带灰色虚线网格，适合柱状图、折线图、散点图
theme_my_stat <- function(base_size = 12, grid_type = "major", ...) {
  # 继承基石 (开启边框)
  t <- theme_my_base(base_size = base_size, border = TRUE, grid = FALSE) +
    theme(
      # 特有的虚线网格
      panel.grid.major = element_line(
        colour = "grey85",
        linetype = "dashed",
        linewidth = 0.5
      ),
      # 图例默认在右侧
      legend.position = "right"
    )

  # 网格控制
  if (grid_type == "horizontal") {
    t <- t + theme(panel.grid.major.x = element_blank())
  } else if (grid_type == "none") {
    t <- t + theme(panel.grid.major = element_blank())
  }

  # 允许用户覆盖 (比如把图例改到底部)
  return(t + theme(...))
}

# ==============================================================================
# 3. Level 3: 智能配色 (Color Scales)
# ==============================================================================
# 内置色板库
my_palettes <- list(
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
  )
)

#' 智能填充色 (自动插值)
scale_fill_scp <- function(palette = "NPG", reverse = FALSE, ...) {
  pal <- my_palettes[[palette]] %||% my_palettes[["NPG"]]
  if (reverse) {
    pal <- rev(pal)
  }

  # 插值逻辑：颜色不够时自动渐变
  pal_fun <- function(n) {
    if (n <= length(pal)) pal[1:n] else colorRampPalette(pal)(n)
  }
  discrete_scale("fill", paste0("scp_", palette), palette = pal_fun, ...)
}

#' 智能线条色 (自动插值)
scale_color_scp <- function(palette = "NPG", reverse = FALSE, ...) {
  pal <- my_palettes[[palette]] %||% my_palettes[["NPG"]]
  if (reverse) {
    pal <- rev(pal)
  }
  pal_fun <- function(n) {
    if (n <= length(pal)) pal[1:n] else colorRampPalette(pal)(n)
  }
  discrete_scale("colour", paste0("scp_", palette), palette = pal_fun, ...)
}
