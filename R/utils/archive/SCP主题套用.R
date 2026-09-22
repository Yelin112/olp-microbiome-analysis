library(ggplot2)
library(grid)

# ==============================================================================
# 1. SCP Theme (主题)
# ==============================================================================
theme_scp <- function(aspect.ratio = NULL, base_size = 12, ...) {
  text_size_scale <- base_size / 12

  # 定义基础样式
  t <- theme(
    # 基础文本设置
    text = element_text(size = 12 * text_size_scale, color = "black"),

    # 标题与副标题
    plot.title = element_text(
      size = 14 * text_size_scale,
      colour = "black",
      vjust = 1,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 13 * text_size_scale,
      hjust = 0,
      margin = margin(b = 5)
    ),

    # 背景与边框 (核心风格：白底 + 黑框)
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 1
    ), # 使用全框
    axis.line = element_blank(), # 隐藏默认坐标轴线，使用panel.border代替

    # 坐标轴
    axis.title = element_text(size = 13 * text_size_scale, colour = "black"),
    axis.text = element_text(size = 12 * text_size_scale, colour = "black"),

    # 分面 (Facets)
    strip.text = element_text(
      size = 12.5 * text_size_scale,
      colour = "black",
      margin = margin(4, 4, 4, 4)
    ),
    strip.background = element_rect(fill = "transparent", linetype = 0), # 透明分面背景
    strip.placement = "outside",

    # 图例
    legend.title = element_text(
      size = 12 * text_size_scale,
      colour = "black",
      hjust = 0
    ),
    legend.text = element_text(size = 11 * text_size_scale, colour = "black"),
    legend.key = element_rect(fill = "transparent", color = "transparent"), # 图例key透明
    legend.background = element_blank(),
    legend.key.size = unit(12, "pt"),

    aspect.ratio = aspect.ratio
  )

  # 允许传入额外的 theme 参数覆盖默认设置
  return(t + theme(...))
}

# ==============================================================================
# 2. SCP Palette (配色)
# ==============================================================================

# 定义一个本地的常用色板列表 (替代原包中的 SCP::palette_list)
my_palettes <- list(
  # 离散色板 (Discrete)
  "Paired" = c(
    "#A6CEE3",
    "#1F78B4",
    "#B2DF8A",
    "#33A02C",
    "#FB9A99",
    "#E31A1C",
    "#FDBF6F",
    "#FF7F00",
    "#CAB2D6",
    "#6A3D9A",
    "#FFFF99",
    "#B15928"
  ),
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
    "#808180",
    "#1B1919"
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

  # 连续色板 (Continuous) - 适合热图或梯度
  "Spectral" = c(
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FFFFBF",
    "#FEE08B",
    "#FDAE61",
    "#F46D43",
    "#D53E4F",
    "#9E0142"
  ),
  "RdBu" = c(
    "#2166AC",
    "#4393C3",
    "#92C5DE",
    "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7",
    "#F4A582",
    "#D6604D",
    "#B2182B"
  ),
  "Viridis" = c(
    "#440154",
    "#482878",
    "#3E4A89",
    "#31688E",
    "#26828E",
    "#1F9E89",
    "#35B779",
    "#6DCD59",
    "#B4DE2C",
    "#FDE725"
  )
)

# 核心配色函数
palette_scp <- function(
  x = NULL,
  palette = "NPG",
  n = NULL,
  type = "auto",
  reverse = FALSE
) {
  # 获取色板颜色
  if (!palette %in% names(my_palettes)) {
    warning(paste("Palette", palette, "not found. Using 'NPG'."))
    pal_colors <- my_palettes[["NPG"]]
  } else {
    pal_colors <- my_palettes[[palette]]
  }

  if (reverse) {
    pal_colors <- rev(pal_colors)
  }

  # 确定数据类型
  if (type == "auto") {
    if (is.numeric(x) && length(unique(x)) > 20) {
      # 简单判断连续变量
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }

  # 逻辑处理
  if (type == "discrete") {
    if (is.null(x)) {
      n_needed <- if (is.null(n)) length(pal_colors) else n
    } else {
      x <- as.factor(x)
      n_needed <- nlevels(x)
    }

    if (n_needed <= length(pal_colors)) {
      out_cols <- pal_colors[1:n_needed]
    } else {
      # 如果需要的颜色多于色板已有颜色，则进行插值
      out_cols <- colorRampPalette(pal_colors)(n_needed)
    }

    if (!is.null(x)) {
      names(out_cols) <- levels(x)
    }
    return(out_cols)
  } else if (type == "continuous") {
    # 连续型返回一个插值函数或指定数量的颜色
    if (is.null(n)) {
      n <- 100
    }
    return(colorRampPalette(pal_colors)(n))
  }
}


# --- 准备数据 ---
library(ggplot2)
data("mpg")

# --- 1. 散点图示例 (离散配色) ---
p1 <- ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "常规 ggplot2 默认样式", subtitle = "Scatter Plot")

p2 <- ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 3, alpha = 0.8) +
  # 使用提取的主题
  theme_scp(base_size = 14, legend.position = "right") +
  # 使用提取的配色 (自动匹配 class 的水平)
  scale_color_manual(values = palette_scp(mpg$class, palette = "NPG")) +
  labs(title = "SCP 样式美化", subtitle = "Scatter Plot with NPG Palette")

# --- 2. 箱线图示例 (填充配色) ---
p3 <- ggplot(mpg, aes(x = manufacturer, y = hwy, fill = manufacturer)) +
  geom_boxplot() +
  theme_scp() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 额外调整x轴标签
  # 这里的 manufacturer 种类很多(15个)，NPG色板只有9个色，函数会自动插值生成足够颜色
  scale_fill_manual(values = palette_scp(mpg$manufacturer, palette = "NPG")) +
  labs(title = "自动插值配色", subtitle = "Levels > Palette Colors")

# --- 3. 热图/连续变量示例 ---
# 构造一个虚构数据
df_heat <- expand.grid(X = 1:10, Y = 1:10)
df_heat$Val <- runif(100)

p4 <- ggplot(df_heat, aes(x = X, y = Y, fill = Val)) +
  geom_tile() +
  theme_scp() +
  # 使用连续型配色 (Spectral)
  scale_fill_gradientn(
    colors = palette_scp(
      type = "continuous",
      palette = "Spectral",
      reverse = TRUE
    )
  ) +
  labs(title = "连续型配色", subtitle = "Heatmap Style")

# --- 拼图展示 ---
# 如果您安装了 patchwork，可以组合查看效果
if (require(patchwork)) {
  (p1 | p2) / (p3 | p4)
} else {
  print(p2)
  print(p3)
}


library(ggplot2)

#' 我的核心主题函数
#' @param base_size 基础字号，默认12。改为 20 即可用于 PPT。
#' @param border 是否显示全框（学术风），默认为 TRUE。FALSE 则只显示坐标轴线。
#' @param grid 是否显示网格线，默认为 FALSE (SCP风格偏向无网格)。
#' @param aspect.ratio 强制长宽比，默认为 NULL。
#' @param ... 其他传递给 theme() 的参数
theme_my_base <- function(
  base_size = 12,
  border = TRUE,
  grid = FALSE,
  aspect.ratio = NULL,
  ...
) {
  # 1. 核心逻辑：计算缩放比例 (SCP 的精髓)
  # 所有的元素大小都基于这个 scale，保证整体协调
  scale <- base_size / 12

  # 2. 定义基础 element (白底)
  # 我们基于 theme_minimal 进行修改，你也可以改用 theme_classic
  t <- theme_minimal(base_size = base_size) +
    theme(
      # --- 文本设置 (自动缩放) ---
      text = element_text(color = "black", size = base_size),
      plot.title = element_text(size = 14 * scale, face = "bold", hjust = 0),
      plot.subtitle = element_text(
        size = 12 * scale,
        color = "grey30",
        margin = margin(b = 10)
      ),
      axis.title = element_text(size = 13 * scale),
      axis.text = element_text(size = 12 * scale, color = "black"),
      legend.title = element_text(size = 12 * scale, face = "bold"),
      legend.text = element_text(size = 11 * scale),

      # --- 背景与边框 ---
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_blank(), # 图例背景透明
      legend.key = element_rect(fill = "transparent", color = NA),

      # --- 长宽比 ---
      aspect.ratio = aspect.ratio
    )

  # 3. 条件控制：边框风格 (仿 SCP 的 panel_fix)
  if (border) {
    t <- t +
      theme(
        panel.border = element_rect(
          fill = "transparent",
          color = "black",
          linewidth = 1
        ),
        axis.line = element_blank() # 有框就不需要轴线了
      )
  } else {
    t <- t +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8)
      )
  }

  # 4. 条件控制：网格线
  if (!grid) {
    t <- t + theme(panel.grid = element_blank())
  }

  # 5. 允许用户通过 ... 传入额外的 theme 设置进行覆盖
  return(t + theme(...))
}


# 场景 1: 幻灯片深色模式 (Dark Mode for PPT)
theme_my_ppt <- function(base_size = 16, ...) {
  theme_my_base(base_size = base_size, border = FALSE, ...) +
    theme(
      plot.background = element_rect(fill = "#2b2b2b", color = NA),
      panel.background = element_rect(fill = "#2b2b2b", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "grey90"),
      panel.border = element_rect(color = "white", fill = NA),
      legend.text = element_text(color = "white")
    )
}

# 场景 2: 降维图专用 (Clean Mode)
# 类似于 Seurat 的 DimPlot 风格，无轴、无字、无框
theme_my_clean <- function(...) {
  theme_void() +
    theme(
      aspect.ratio = 1, # 强制正方形
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "right"
    )
}

# 场景 3: 只有横向网格线 (适合做 BarPlot 比较)
theme_my_bar <- function(...) {
  theme_my_base(grid = FALSE, border = FALSE, ...) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      axis.ticks.x = element_blank()
    )
}


# 1. 旋转 X 轴标签 (常用！)
# 用法: + rotate_x(45)
rotate_x <- function(angle = 45, hjust = 1, vjust = 1) {
  theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust))
}

# 2. 移除图例 (常用！)
# 用法: + no_legend()
no_legend <- function() {
  theme(legend.position = "none")
}

# 3. 调整图例位置到图内部 (节省空间)
# 用法: + legend_inside(0.8, 0.8)
legend_inside <- function(x = 0.8, y = 0.8) {
  theme(
    legend.position = c(x, y),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA) # 半透明白底防遮挡
  )
}

theme_my_stat <- function(base_size = 12, grid_type = "major", ...) {
  # 1. 继承基石主题 (强制开启边框 border=TRUE)
  t <- theme_my_base(base_size = base_size, border = TRUE, grid = FALSE) +
    theme(
      # SCP 风格的核心：灰色虚线网格
      panel.grid.major = element_line(
        colour = "grey85",
        linetype = "dashed",
        linewidth = 0.5
      ),
      panel.grid.minor = element_blank(),

      # 坐标轴标签微调 (StatPlot 默认稍微紧凑一点)
      axis.text = element_text(color = "black", size = base_size * 0.9),

      # 图例优化
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  # 2. 网格线控制逻辑 (让它更灵活)
  if (grid_type == "horizontal") {
    t <- t + theme(panel.grid.major.x = element_blank())
  } else if (grid_type == "none") {
    t <- t + theme(panel.grid.major = element_blank())
  } else if (grid_type == "both") {
    t <- t +
      theme(
        panel.grid.minor = element_line(colour = "grey90", linetype = "dotted")
      )
  }

  # 3. 允许用户覆盖
  return(t + theme(...))
}


# 定义 SCP 常用色板库
scp_palettes <- list(
  # 经典离散色板
  "Paired" = c(
    "#A6CEE3",
    "#1F78B4",
    "#B2DF8A",
    "#33A02C",
    "#FB9A99",
    "#E31A1C",
    "#FDBF6F",
    "#FF7F00",
    "#CAB2D6",
    "#6A3D9A",
    "#FFFF99",
    "#B15928"
  ),
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
    "#808180",
    "#1B1919"
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
  "JCO" = c(
    "#0073C2",
    "#EFC000",
    "#868686",
    "#CD534C",
    "#7AA6DC",
    "#003C67",
    "#8F7700",
    "#3B3B3B",
    "#A73030",
    "#4A6990"
  )
)

#' SCP 风格填充色 Scale
#' @param palette 色板名称，默认为 "NPG" (Nature 风格)。可选 "Paired", "AAAS", "NEJM", "JCO"。
#' @param reverse 是否反转颜色顺序
#' @param ... 传递给 discrete_scale 的其他参数
scale_fill_scp <- function(palette = "NPG", reverse = FALSE, ...) {
  # 获取色值
  pal <- scp_palettes[[palette]]
  if (is.null(pal)) {
    warning("Palette not found. Using NPG.")
    pal <- scp_palettes[["NPG"]]
  }
  if (reverse) {
    pal <- rev(pal)
  }

  # 定义插值函数 (这是 SCP 不会报错的核心)
  pal_fun <- function(n) {
    if (n <= length(pal)) {
      pal[1:n]
    } else {
      colorRampPalette(pal)(n)
    }
  }

  discrete_scale("fill", paste0("scp_", palette), palette = pal_fun, ...)
}

#' SCP 风格线条色 Scale (同上，用于 color=...)
scale_color_scp <- function(palette = "NPG", reverse = FALSE, ...) {
  pal <- scp_palettes[[palette]]
  if (is.null(pal)) {
    pal <- scp_palettes[["NPG"]]
  }
  if (reverse) {
    pal <- rev(pal)
  }

  pal_fun <- function(n) {
    if (n <= length(pal)) pal[1:n] else colorRampPalette(pal)(n)
  }

  discrete_scale("colour", paste0("scp_", palette), palette = pal_fun, ...)
}


library(ggplot2)

# 准备数据：模拟一个多组柱状图
df_bar <- data.frame(
  Group = rep(c("Control", "TreatA", "TreatB", "TreatC"), each = 3),
  Sample = rep(c("S1", "S2", "S3"), 4),
  Value = c(10, 12, 11, 15, 18, 16, 8, 9, 7, 20, 22, 21)
)

# --- 绘图 ---
p <- ggplot(df_bar, aes(x = Group, y = Value, fill = Group)) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    color = "black",
    width = 0.7
  ) +

  # 1. 应用 SCP StatPlot 风格主题
  #    base_size = 14: 稍微调大一点字体
  #    grid_type = "horizontal": 柱状图通常只需要横向网格线
  theme_my_stat(base_size = 14, grid_type = "horizontal") +

  # 2. 应用 SCP NPG 配色
  scale_fill_scp(palette = "NPG") +

  # 3. 额外微调 (可选)
  labs(
    title = "Bar Plot in SCP Style",
    subtitle = "Recreated using modular themes"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) # 像 StatPlot 一样让柱子贴底

print(p)


p3 + theme_my_stat(base_size = 14, grid_type = "horizontal")


# 定义您的“经济学人”风格主题
theme_economist_style <- function(base_size = 12, ...) {
  # 1. 继承：先拿一个带网格的基石
  #    注意：我们传入 grid=FALSE，因为我们要自定义特殊的白色网格
  theme_my_base(base_size = base_size, border = FALSE, grid = FALSE) +

    # 2. 覆盖：加上经济学人的特色
    theme(
      # 背景色
      plot.background = element_rect(fill = "#d5e4eb", color = NA),
      panel.background = element_rect(fill = "#d5e4eb", color = NA),

      # 特殊的白色网格线 (Grid lines)
      panel.grid.major.y = element_line(color = "white", linewidth = 1.2),
      panel.grid.major.x = element_blank(), # 只要横线

      # 坐标轴文字
      axis.text = element_text(color = "grey30"),
      axis.title.y = element_text(margin = margin(r = 10)),

      # 图例
      legend.position = "top",
      legend.justification = "left", # 图例靠左对齐

      ... # 接收外部参数
    )
}


p3 + theme_economist_style(base_size = 14)

p3 + theme_my_base()


#' 旋转 X 轴文字
#' @param angle 旋转角度，推荐 45 或 90
#' @param hjust 水平对齐 (自动调整，通常不需要动)
#' @param vjust 垂直对齐 (自动调整，通常不需要动)
rotate_x_text <- function(angle = 45, hjust = 1, vjust = 1) {
  if (angle == 90) {
    vjust <- 0.5
  } # 90度时垂直居中更好看
  theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust))
}

#' 移除 X 轴或 Y 轴 (用于热图或特定展示)
remove_x_axis <- function() {
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
}
remove_y_axis <- function() {
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
}

#' 仅移除轴标题 (保留刻度)
remove_titles <- function() {
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
}

p1 + theme_my_stat()


library(ggplot2)
library(grid)
library(gtable)

#' 固定 ggplot2 绘图面板 (Panel) 的物理尺寸
#'
#' @description
#' 这是一个管道友好的函数，用于强制设定 ggplot 绘图核心区域的高度和宽度。
#' 无论轴标签、标题、图例占据多少空间，Panel 的尺寸都将严格固定。
#'
#' @param p ggplot2 对象
#' @param width Panel 的宽度 (默认单位 cm，也可以传入 unit 对象)
#' @param height Panel 的高度 (默认单位 cm，也可以传入 unit 对象)
#' @param margin 是否添加额外的边距 (默认为 NULL)
#'
#' @return 返回一个 gtable 对象，可以直接 print 显示或用于 ggsave
#' @export
set_panel_size <- function(
  p = NULL,
  width = unit(4, "cm"),
  height = unit(4, "cm"),
  margin = NULL
) {
  # 1. 确保输入是 ggplot 对象并转换为 gtable
  if (inherits(p, "ggplot")) {
    g <- ggplotGrob(p)
  } else if (inherits(p, "gtable")) {
    g <- p
  } else {
    stop("Input must be a ggplot object or a gtable.")
  }

  # 2. 处理单位 (如果用户输入数字，默认当做 cm)
  if (!is.unit(width)) {
    width <- unit(width, "cm")
  }
  if (!is.unit(height)) {
    height <- unit(height, "cm")
  }

  # 3. 查找 Panel 在 gtable 中的位置
  #    g$layout 包含了布局信息，我们需要找到名为 "panel" 的行和列
  panels <- grep("panel", g$layout$name)
  if (length(panels) == 0) {
    warning("No panels found in the plot.")
    return(g)
  }

  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])

  # 4. 强制修改宽度和高度
  #    注意：对于分面图 (facet)，有多个 panel，我们需要同时修改它们
  if (length(panel_index_w) == 1) {
    g$widths[panel_index_w] <- width
  } else {
    # 针对分面图，修改所有 panel 列的宽度
    g$widths[panel_index_w] <- rep(width, length(panel_index_w))
  }

  if (length(panel_index_h) == 1) {
    g$heights[panel_index_h] <- height
  } else {
    # 针对分面图，修改所有 panel 行的高度
    g$heights[panel_index_h] <- rep(height, length(panel_index_h))
  }

  # 5. (可选) 添加边距
  if (!is.null(margin)) {
    # 这里可以添加简单的 margin 处理，但在 gtable 层面上比较繁琐
    # 通常建议在 ggplot theme(plot.margin) 中设置
  }

  # 6. 为返回对象添加一个类，以便我们可以定义自定义的 print 方法
  class(g) <- c("fixed_panel_plot", class(g))
  return(g)
}

#' 定义 print 方法，让它能在 RStudio 中直接显示
#' @export
print.fixed_panel_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
}

p
p |> set_panel_size(width = 3, height = 5)


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
  return(invisible(g))
}

p + fix_panel(width = 4, height = 4)

p |> set_panel_size(width = 3, height = 5)


library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(patchwork) # 用于拼图

# --- 1. 必要的辅助符号 ---
# 定义空值合并运算符 (如果 a 为空则用 b)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- 2. 简化的配色函数 (如果您之前没定义 scale_fill_scp) ---
# 如果您已经有了之前的 my_themes.R，这部分可以跳过
if (!exists("palette_scp")) {
  palette_scp <- function(x = NULL, palette = "NPG", ...) {
    # 这里放一个简单的兜底，实际建议用您之前的 scale_fill_scp
    cols <- c(
      "#E64B35",
      "#4DBBD5",
      "#00A087",
      "#3C5488",
      "#F39B7F",
      "#8491B4",
      "#91D1C2"
    )
    if (is.null(x)) {
      return(cols)
    }
    colorRampPalette(cols)(length(unique(x)))
  }
}

#' 通用统计绘图函数 (源自 SCP::StatPlot)
#'
#' @param data 一个普通的数据框 (data.frame)
#' @param stat.by 需要统计的列名 (例如 "CellType" 或 "Gender")，用作填充色
#' @param group.by 分组列名 (例如 "Sample")，用作 X 轴
#' @param split.by 分面列名 (可选)
#' @param plot_type 图表类型: "bar", "pie", "ring", "rose", "area"
#' @param stat_type 统计方式: "percent" (百分比) 或 "count" (计数)
#' @param position 堆叠方式: "stack" (堆叠) 或 "dodge" (并排)
#' @param palette 配色方案名称 (传递给 palette_scp)
#' @param label 是否显示数值标签
#' @param ... 传递给 theme() 的其他参数
#' @export
Universal_StatPlot <- function(
  data,
  stat.by,
  group.by = NULL,
  split.by = NULL,
  plot_type = c("bar", "rose", "ring", "pie", "area"),
  stat_type = c("percent", "count"),
  position = c("stack", "dodge"),
  palette = "NPG",
  label = FALSE,
  title = NULL,
  xlab = NULL,
  ylab = NULL
) {
  # --- 参数校验与预处理 ---
  plot_type <- match.arg(plot_type)
  stat_type <- match.arg(stat_type)
  position <- match.arg(position)

  if (is.null(group.by)) {
    group.by <- "All_Groups"
    data[[group.by]] <- factor("All")
    xlab <- ""
  }
  if (is.null(split.by)) {
    split.by <- "All_Split"
    data[[split.by]] <- factor("All")
  }

  # 确保指定的列存在
  cols_needed <- c(stat.by, group.by, split.by)
  if (!all(cols_needed %in% colnames(data))) {
    stop(
      "指定的列名不存在于数据框中: ",
      paste(setdiff(cols_needed, colnames(data)), collapse = ", ")
    )
  }

  # 确保列为因子 (Factor) 以保证绘图顺序
  for (col in cols_needed) {
    if (!is.factor(data[[col]])) {
      data[[col]] <- factor(data[[col]], levels = unique(data[[col]]))
    }
  }

  # --- 核心计算逻辑 (dplyr) ---
  # 这里完全复刻了 SCP 的计算逻辑：先计数，再算比例

  # 1. 拆分数据 (如果有 split.by)
  dat_split <- split(data, data[[split.by]])
  plist <- list()

  # 2. 循环绘图
  for (sp in names(dat_split)) {
    sub_data <- dat_split[[sp]]

    # 聚合统计 (Count & Percent)
    # 类似于 table() 但输出为 dataframe
    df_stat <- sub_data %>%
      group_by(across(all_of(c(stat.by, group.by)))) %>%
      summarise(Freq = n(), .groups = "drop") %>%
      group_by(across(all_of(group.by))) %>%
      mutate(Total = sum(Freq), Percent = Freq / Total) %>%
      ungroup()

    # 根据用户选择决定 Y 轴是 Count 还是 Percent
    if (stat_type == "percent") {
      df_stat$Value <- df_stat$Percent
      y_label_func <- scales::percent
      ylab_default <- "Percentage"
    } else {
      df_stat$Value <- df_stat$Freq
      y_label_func <- scales::number
      ylab_default <- "Count"
    }

    # --- 绘图构建 (ggplot2) ---
    # 定义基础映射
    p <- ggplot(
      df_stat,
      aes(x = .data[[group.by]], y = Value, fill = .data[[stat.by]])
    )

    # 2.1 几何图层选择
    if (plot_type == "bar") {
      p <- p +
        geom_col(
          position = if (position == "stack") {
            "stack"
          } else {
            position_dodge2(preserve = "single")
          },
          width = 0.8,
          color = "black",
          size = 0.2
        )
    } else if (plot_type == "area") {
      p <- p +
        geom_area(
          aes(group = .data[[stat.by]]),
          position = "stack",
          alpha = 0.8,
          color = "black"
        )
    } else if (plot_type %in% c("pie", "ring", "rose")) {
      p <- p + geom_col(width = 1, color = "white", size = 0.2)
    }

    # 2.2 坐标系变换 (针对圆环/玫瑰图)
    if (plot_type == "rose") {
      p <- p + coord_polar(theta = "x")
    } else if (plot_type == "pie") {
      p <- p + coord_polar(theta = "y") + theme_void() # 饼图通常需要极简主题
    } else if (plot_type == "ring") {
      p <- p + coord_polar(theta = "y") + xlim(c(0.5, 2.5)) + theme_void() # 环形图通过挖空圆心实现
      # 对于环形图，x 轴需要设为常数才能形成环
      # 这里简化处理：通常 ring 图不适合展示多组 group，建议单组使用
    }

    # 2.3 标尺与标签
    p <- p +
      scale_y_continuous(
        labels = y_label_func,
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        title = if (split.by != "All_Split") sp else title,
        x = xlab %||% group.by,
        y = ylab %||% ylab_default,
        fill = stat.by
      )

    # 2.4 标签显示 (ggrepel)
    if (label) {
      p <- p +
        geom_text_repel(
          aes(
            label = if (stat_type == "percent") {
              scales::percent(Value, accuracy = 0.1)
            } else {
              Value
            }
          ),
          position = if (position == "stack") {
            position_stack(vjust = 0.5)
          } else {
            position_dodge(width = 0.8)
          },
          size = 3,
          show.legend = FALSE
        )
    }

    # 2.5 应用您之前的主题 (如果有)
    if (exists("theme_my_stat")) {
      p <- p + theme_my_stat()
    } else {
      p <- p + theme_classic() # 兜底主题
    }

    # 2.6 应用配色 (这里调用 palette_scp 逻辑，或者直接用 ggplot 默认)
    # 为了演示，我们用简单的 scale_fill_manual
    # 实际使用建议换成您的 scale_fill_scp()
    if (exists("scale_fill_scp")) {
      p <- p + scale_fill_scp(palette = palette)
    }

    plist[[sp]] <- p
  }

  # 3. 输出结果
  if (length(plist) == 1) {
    return(plist[[1]])
  }

  # 如果有多个分面，用 patchwork 拼起来
  return(wrap_plots(plist))
}


# --- 准备数据 ---
df <- mtcars
df$cyl <- as.factor(df$cyl) # 气缸数 (假设为细胞类型)
df$gear <- as.factor(df$gear) # 档位数 (假设为样本分组)
df$am <- factor(df$am, labels = c("Auto", "Manual")) # 变速箱 (假设为分面)

# --- 场景 1: 基础堆叠柱状图 (算百分比) ---
# 自动计算每个 gear 中不同 cyl 的比例
Universal_StatPlot(
  df,
  stat.by = "cyl", # 填充颜色 (CellType)
  group.by = "gear", # X轴分组 (Sample)
  stat_type = "percent",
  plot_type = "bar",
  label = TRUE
) # 显示百分比数值

# --- 场景 2: 玫瑰图 (Rose Plot) ---
Universal_StatPlot(df, stat.by = "cyl", group.by = "gear", plot_type = "rose")

# --- 场景 3: 自动分面 (Split by) ---
# 按照 am (Auto/Manual) 自动把图分开画
Universal_StatPlot(
  df,
  stat.by = "cyl",
  group.by = "gear",
  split.by = "am", # <--- 这里
  position = "dodge", # 并排柱状图
  stat_type = "count"
)


library(ggplot2)
library(dplyr)

#' 添加背景分块染色图层 (仿 SCP bg.by)
#'
#' @description
#' 根据 X 轴的分组，在背景绘制矩形色块。
#' 必须在 geom_col/geom_point 之前添加此图层，以防遮挡数据。
#'
#' @param data 绘图所用的数据框
#' @param x_var X 轴的列名 (必须是 Factor 或 Character)
#' @param bg_var 背景分组的列名 (bg.by)
#' @param palette 使用的色板名称 (如 "Pastel1", "Set2" 或自定义颜色向量)
#' @param alpha 透明度
#'
#' @return geom_rect 图层
#' @export
add_shading_layer <- function(
  data,
  x_var,
  bg_var,
  palette = "Pastel1",
  alpha = 0.2
) {
  # 1. 确保 X 轴是因子，并转换为数字坐标 (1, 2, 3...)
  if (!is.factor(data[[x_var]])) {
    data[[x_var]] <- factor(data[[x_var]], levels = unique(data[[x_var]]))
  }

  # 2. 提取坐标映射关系
  #    计算每个 bg_group 包含哪些 x_group，以及它们的 X 轴数值范围
  bg_data <- data %>%
    select(all_of(c(x_var, bg_var))) %>%
    distinct() %>%
    mutate(x_num = as.numeric(.data[[x_var]])) %>%
    group_by(.data[[bg_var]]) %>%
    summarise(
      xmin = min(x_num) - 0.5,
      xmax = max(x_num) + 0.5,
      .groups = "drop"
    )

  # 3. 处理颜色 (避免占用 ggplot 的 scale_fill)
  #    我们手动计算好颜色值，直接传给 fill 参数
  n_groups <- nrow(bg_data)

  if (
    length(palette) == 1 && palette %in% rownames(RColorBrewer::brewer.pal.info)
  ) {
    # 如果是 Brewer 色板名
    if (n_groups <= 2) {
      cols <- RColorBrewer::brewer.pal(3, palette)[1:n_groups]
    } else {
      cols <- colorRampPalette(RColorBrewer::brewer.pal(
        min(n_groups, 8),
        palette
      ))(n_groups)
    }
  } else if (length(palette) >= n_groups) {
    # 如果是自定义颜色向量
    cols <- palette[1:n_groups]
  } else {
    # 兜底
    cols <- scales::hue_pal()(n_groups)
  }

  # 将颜色加入数据框
  bg_data$fill_color <- cols

  # 4. 返回 geom_rect 图层
  #    注意：inherit.aes = FALSE 防止受到主图 aes 的干扰
  geom_rect(
    data = bg_data,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = bg_data$fill_color, # 直接使用计算好的颜色码
    alpha = alpha,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}


library(ggplot2)
library(ggpubr) # 核心统计绘图包
library(dplyr)
library(reshape2) # 用于数据重塑
library(patchwork) # 用于拼图


#' 通用特征统计绘图函数 (源自 SCP::FeatureStatPlot)
#'
#' @param data 普通数据框 (data.frame)
#' @param value.var 数值列名 (例如 "Expression" 或 "Value")
#' @param group.by 分组列名 (例如 "Group" 或 "Treatment")，用于 X 轴
#' @param split.by 分面列名 (可选，用于分面展示)
#' @param plot_type 图表类型: "violin", "box", "bar", "dot" (均值点图)
#' @param add_box 是否在小提琴图里叠加箱线图 (默认 FALSE)
#' @param add_point 是否叠加散点 (默认 FALSE)
#' @param add_trend 是否添加趋势线 (默认 FALSE, 适合有序分组)
#' @param add_stat 添加统计检验: "none", "mean" (显示均值), "median", 或统计方法名 "t.test", "wilcox.test", "anova", "kruskal.test"
#' @param comparisons 指定比较组 (list), 例如 list(c("A", "B"), c("C", "D"))。如果为 NULL 且 add_stat 为检验方法，则自动做两两比较。
#' @param palette 配色方案 (传递给 scale_fill_scp)
#' @param title 标题
#' @param ... 传递给 theme() 的参数
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  split.by = NULL,
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
  theme_use = theme_my_stat, # 默认使用你的自定义主题
  ...
) {
  # --- 1. 参数预处理 ---
  plot_type <- match.arg(plot_type)
  if (is.null(xlab)) {
    xlab <- group.by
  }
  if (is.null(ylab)) {
    ylab <- value.var
  }

  # 确保列存在
  if (!value.var %in% colnames(data)) {
    stop("Value variable not found.")
  }
  if (!group.by %in% colnames(data)) {
    stop("Group variable not found.")
  }

  # 转换因子以保持顺序
  if (!is.factor(data[[group.by]])) {
    data[[group.by]] <- factor(
      data[[group.by]],
      levels = unique(data[[group.by]])
    )
  }

  # --- 2. 基础绘图构建 ---
  p <- ggplot(
    data,
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[group.by]])
  )

  # --- 3. 添加几何图层 (核心可视化) ---

  # A. 主图层 (Violin / Box / Bar / Dot)
  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3
      )
  } else if (plot_type == "box") {
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black"
      )
  } else if (plot_type == "bar") {
    # 柱状图通常展示均值 + 误差棒
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  }

  # B. 叠加层 (Box / Point / Trend)
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8
      )
  }

  if (isTRUE(add_point)) {
    p <- p +
      geom_jitter(
        width = 0.2,
        size = 1.5,
        alpha = 0.6,
        color = "grey30",
        shape = 16
      )
  }

  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(group = 1),
        color = "black",
        linewidth = 1,
        linetype = "dashed"
      )
  }

  # --- 4. 添加统计信息 (核心统计) ---

  # 简单统计标注 (均值/中位数)
  if (add_stat %in% c("mean", "median")) {
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
  }

  # 假设检验 (t.test, wilcox 等) - 依赖 ggpubr
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    # 如果没指定比较组，且是多组比较 (anova/kruskal)，显示全局 P 值
    if (is.null(comparisons) && length(unique(data[[group.by]])) > 2) {
      p <- p + stat_compare_means(method = add_stat, label.y.npc = "top")
    } else {
      # 如果指定了比较组，或者自动进行两两比较
      p <- p +
        stat_compare_means(
          comparisons = comparisons,
          method = add_stat,
          label = "p.signif"
        )
    }
  }

  # --- 5. 分面与外观 ---
  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  # 应用主题与配色
  if (exists("theme_my_stat")) {
    p <- p + theme_use(...) # 使用您自定义的主题
  } else {
    p <- p + theme_classic()
  }

  if (exists("scale_fill_scp")) {
    p <- p + scale_fill_scp(palette = palette)
  }

  # 标签
  p <- p + labs(title = title, x = xlab, y = ylab, fill = group.by)

  return(p)
}


# --- 1. 模拟数据 ---
df_mouse <- data.frame(
  Group = rep(c("Control", "Drug_A", "Drug_B"), each = 20),
  Weight = c(rnorm(20, 20, 2), rnorm(20, 25, 2), rnorm(20, 22, 2)),
  Gender = rep(rep(c("Male", "Female"), each = 10), 3)
)

# --- 2. 基础小提琴图 + 箱线图 + 抖动点 ---
Universal_FeatureStatPlot(
  data = df_mouse,
  value.var = "Weight",
  group.by = "Group",
  plot_type = "violin",
  add_box = TRUE, # 叠加箱线图
  add_point = TRUE, # 叠加散点
  title = "Mouse Weight Distribution"
)

# --- 3. 带有统计检验的箱线图 ---
# 自动进行 Wilcoxon 两两比较，并标注显著性 (*, **)
Universal_FeatureStatPlot(
  data = df_mouse,
  value.var = "Weight",
  group.by = "Group",
  plot_type = "box",
  add_stat = "wilcox.test", # 添加统计检验
  comparisons = list(c("Control", "Drug_A"), c("Control", "Drug_B")), # 指定比较对象
  palette = "AAAS" # 使用 Science 风格配色
)

# --- 4. 分面展示 + 柱状图 + 趋势线 ---
Universal_FeatureStatPlot(
  data = df_mouse,
  value.var = "Weight",
  group.by = "Group",
  split.by = "Gender", # 按性别分面
  plot_type = "bar", # 柱状图 (带误差棒)
  add_trend = TRUE, # 添加趋势线
  add_stat = "anova", # 显示 ANOVA P值
  theme_use = theme_my_stat # 使用您的自定义主题
)


library(ggplot2)
library(ggpubr)
library(dplyr)

#' 通用特征统计绘图函数 (集成斑马纹背景版)
#'
#' @param data 普通数据框
#' @param value.var Y轴数值列名
#' @param group.by X轴分组列名
#' @param split.by 分面列名
#' @param add_bg 是否添加交替的灰白背景 (斑马纹), 默认 FALSE
#' @param bg_color 背景斑马纹的颜色, 默认 "grey95"
#' @param plot_type 图表类型: "violin", "box", "bar", "dot"
#' @param add_box 是否叠加箱线图
#' @param add_point 是否叠加散点
#' @param add_trend 是否叠加趋势线
#' @param add_stat 统计方法: "none", "mean", "median", "t.test", "wilcox.test", "anova", "kruskal.test"
#' @param comparisons 指定比较组 list(c("A","B"), c("C","D"))
#' @param palette 配色方案
#' @param theme_use 使用的主题函数
#' @param ... 传递给 theme() 的其他参数
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  split.by = NULL,
  add_bg = FALSE, # <--- 新增: 开关
  bg_color = "grey95", # <--- 新增: 颜色
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
  theme_use = theme_my_stat,
  ...
) {
  # --- 1. 参数预处理 ---
  plot_type <- match.arg(plot_type)
  if (is.null(xlab)) {
    xlab <- group.by
  }
  if (is.null(ylab)) {
    ylab <- value.var
  }

  # 确保列存在
  if (!value.var %in% colnames(data)) {
    stop(paste0("Column '", value.var, "' not found."))
  }
  if (!group.by %in% colnames(data)) {
    stop(paste0("Column '", group.by, "' not found."))
  }

  # 转换因子以保持顺序 (非常重要，否则背景会对不齐)
  if (!is.factor(data[[group.by]])) {
    data[[group.by]] <- factor(
      data[[group.by]],
      levels = unique(data[[group.by]])
    )
  }

  # --- 2. 初始化 ggplot 对象 ---
  p <- ggplot(
    data,
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[group.by]])
  )

  # --- 3. [核心新增] 自动添加斑马纹背景 ---
  # 逻辑：必须在所有数据图层之前添加
  if (isTRUE(add_bg)) {
    # 获取 X 轴的分组数量
    n_groups <- nlevels(data[[group.by]])

    # 生成背景数据框：只选取奇数位置 (1, 3, 5...)
    # ggplot2 的离散轴坐标本质是 1, 2, 3...，每个宽度的范围是 [i-0.5, i+0.5]
    bg_df <- data.frame(
      x_pos = seq(1, n_groups, 2) # 选取 1, 3, 5...
    ) %>%
      mutate(
        xmin = x_pos - 0.5,
        xmax = x_pos + 0.5
      )

    # 添加 geom_rect
    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = bg_color,
        inherit.aes = FALSE, # 防止继承主图的 aes
        show.legend = FALSE
      )
  }

  # --- 4. 添加主数据图层 ---

  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3
      )
  } else if (plot_type == "box") {
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black"
      )
  } else if (plot_type == "bar") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  }

  # --- 5. 叠加层 (Box / Point / Trend) ---
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8
      )
  }

  if (isTRUE(add_point)) {
    p <- p +
      geom_jitter(
        width = 0.2,
        size = 1.5,
        alpha = 0.6,
        color = "grey30",
        shape = 16
      )
  }

  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(group = 1),
        color = "black",
        linewidth = 1,
        linetype = "dashed"
      )
  }

  # --- 6. 统计标注 ---
  if (add_stat %in% c("mean", "median")) {
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    if (is.null(comparisons) && length(unique(data[[group.by]])) > 2) {
      p <- p + stat_compare_means(method = add_stat, label.y.npc = "top")
    } else {
      p <- p +
        stat_compare_means(
          comparisons = comparisons,
          method = add_stat,
          label = "p.signif"
        )
    }
  }

  # --- 7. 收尾 (分面、主题、配色) ---
  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  if (exists("theme_my_stat")) {
    p <- p + theme_use(...)
  } else {
    p <- p + theme_classic()
  }

  if (exists("scale_fill_scp")) {
    p <- p + scale_fill_scp(palette = palette)
  }

  p <- p + labs(title = title, x = xlab, y = ylab, fill = group.by)

  return(p)
}


library(ggplot2)
library(ggpubr)
library(dplyr)

#' 通用特征统计绘图函数 (斑马纹背景 + 虚线网格共存版)
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "grey95",
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
  theme_use = theme_my_stat,
  ...
) {
  plot_type <- match.arg(plot_type)
  if (is.null(xlab)) {
    xlab <- group.by
  }
  if (is.null(ylab)) {
    ylab <- value.var
  }

  if (!value.var %in% colnames(data)) {
    stop(paste0("Column '", value.var, "' not found."))
  }
  if (!group.by %in% colnames(data)) {
    stop(paste0("Column '", group.by, "' not found."))
  }

  if (!is.factor(data[[group.by]])) {
    data[[group.by]] <- factor(
      data[[group.by]],
      levels = unique(data[[group.by]])
    )
  }

  # --- 初始化 ---
  p <- ggplot(
    data,
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[group.by]])
  )

  # --- 1. [关键修正] 先画斑马纹背景 ---
  if (isTRUE(add_bg)) {
    n_groups <- nlevels(data[[group.by]])
    bg_df <- data.frame(x_pos = seq(1, n_groups, 2)) %>%
      mutate(xmin = x_pos - 0.5, xmax = x_pos + 0.5)

    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = bg_color,
        inherit.aes = FALSE,
        show.legend = FALSE
      )

    # --- 2. [关键修正] 强制在背景之上重绘虚线网格 ---
    # 因为 geom_rect 会遮挡默认的主题网格线，所以我们需要手动画上去
    # 这里根据 Y 轴数据的范围生成网格线
    y_range <- range(data[[value.var]], na.rm = TRUE)
    # 使用 pretty() 函数生成美观的刻度点，和 ggplot 默认算法一致
    grid_breaks <- pretty(y_range, n = 5)

    p <- p +
      geom_hline(
        yintercept = grid_breaks,
        linetype = "dashed",
        color = "grey85", # 颜色要比背景深一点
        size = 0.5
      )
  }

  # --- 3. 添加数据图层 (和之前一样) ---
  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3
      )
  } else if (plot_type == "box") {
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black"
      )
  } else if (plot_type == "bar") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  }

  # --- 4. 叠加层 ---
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8
      )
  }
  if (isTRUE(add_point)) {
    p <- p +
      geom_jitter(
        width = 0.2,
        size = 1.5,
        alpha = 0.6,
        color = "grey30",
        shape = 16
      )
  }
  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(group = 1),
        color = "black",
        linewidth = 1,
        linetype = "dashed"
      )
  }

  # --- 5. 统计与美化 ---
  if (add_stat %in% c("mean", "median")) {
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    if (is.null(comparisons) && length(unique(data[[group.by]])) > 2) {
      p <- p + stat_compare_means(method = add_stat, label.y.npc = "top")
    } else {
      p <- p +
        stat_compare_means(
          comparisons = comparisons,
          method = add_stat,
          label = "p.signif"
        )
    }
  }

  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  if (exists("theme_my_stat")) {
    # 如果开启了 add_bg，我们已经在上面手动加了网格线
    # 为了避免双重网格线（可能造成重影），这里可以临时关掉主题里的网格
    if (add_bg) {
      # 先应用主题，再覆盖掉默认网格
      p <- p + theme_use(...) + theme(panel.grid.major.y = element_blank())
    } else {
      p <- p + theme_use(...)
    }
  } else {
    p <- p + theme_classic()
  }

  if (exists("scale_fill_scp")) {
    p <- p + scale_fill_scp(palette = palette)
  }

  p <- p + labs(title = title, x = xlab, y = ylab, fill = group.by)
  return(p)
}


library(ggplot2)
library(ggpubr)
library(dplyr)

#' 通用特征统计绘图函数 (斑马纹背景 + 仅横向虚线)
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "grey95",
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
  theme_use = theme_my_stat,
  ...
) {
  plot_type <- match.arg(plot_type)
  if (is.null(xlab)) {
    xlab <- group.by
  }
  if (is.null(ylab)) {
    ylab <- value.var
  }

  # 数据校验
  if (!value.var %in% colnames(data)) {
    stop(paste0("Column '", value.var, "' not found."))
  }
  if (!group.by %in% colnames(data)) {
    stop(paste0("Column '", group.by, "' not found."))
  }
  if (!is.factor(data[[group.by]])) {
    data[[group.by]] <- factor(
      data[[group.by]],
      levels = unique(data[[group.by]])
    )
  }

  # --- 初始化 ---
  p <- ggplot(
    data,
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[group.by]])
  )

  # --- 1. [背景层] 斑马纹 + 手动横向网格 ---
  if (isTRUE(add_bg)) {
    # A. 画背景块
    n_groups <- nlevels(data[[group.by]])
    bg_df <- data.frame(x_pos = seq(1, n_groups, 2)) %>%
      mutate(xmin = x_pos - 0.5, xmax = x_pos + 0.5)

    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = bg_color,
        inherit.aes = FALSE,
        show.legend = FALSE
      )

    # B. [关键] 手动添加横向虚线 (浮在背景上)
    # 计算优化的刻度位置 (和 ggplot 默认的一致)
    y_range <- range(data[[value.var]], na.rm = TRUE)
    grid_breaks <- pretty(y_range, n = 5)

    p <- p +
      geom_hline(
        yintercept = grid_breaks,
        linetype = "dashed",
        color = "white", # 建议：在灰背景上用白色网格(反白)会更高级，或者用 "grey80"
        size = 0.5
      )
  }

  # --- 2. [数据层] ---
  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3
      )
  } else if (plot_type == "box") {
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black"
      )
  } else if (plot_type == "bar") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  }

  # --- 3. [叠加层] ---
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8
      )
  }
  if (isTRUE(add_point)) {
    p <- p +
      geom_jitter(
        width = 0.2,
        size = 1.5,
        alpha = 0.6,
        color = "grey30",
        shape = 16
      )
  }
  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(group = 1),
        color = "black",
        linewidth = 1,
        linetype = "dashed"
      )
  }

  # --- 4. [统计层] ---
  if (add_stat %in% c("mean", "median")) {
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    if (is.null(comparisons) && length(unique(data[[group.by]])) > 2) {
      p <- p + stat_compare_means(method = add_stat, label.y.npc = "top")
    } else {
      p <- p +
        stat_compare_means(
          comparisons = comparisons,
          method = add_stat,
          label = "p.signif"
        )
    }
  }

  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  # --- 5. [主题层] 最终清洗 ---
  if (exists("theme_my_stat")) {
    if (add_bg) {
      # 核心修改：如果加了背景，就应用主题但强制清空所有默认网格
      # 这样只保留我们刚才 geom_hline 画的横线
      p <- p + theme_use(...) + theme(panel.grid = element_blank())
    } else {
      p <- p + theme_use(...)
    }
  } else {
    p <- p + theme_classic()
  }

  if (exists("scale_fill_scp")) {
    p <- p + scale_fill_scp(palette = palette)
  }

  p <- p + labs(title = title, x = xlab, y = ylab, fill = group.by)
  return(p)
}


library(ggplot2)

# --- 1. 模拟数据 ---
df_test <- data.frame(
  Group = factor(rep(
    c("Ctrl", "Treat1", "Treat2", "Treat3", "Treat4"),
    each = 20
  )),
  Value = c(
    rnorm(20, 10),
    rnorm(20, 12),
    rnorm(20, 15),
    rnorm(20, 11),
    rnorm(20, 18)
  )
)

# --- 2. 基础效果：自动交替背景 ---
Universal_FeatureStatPlot(
  data = df_test,
  value.var = "Value",
  group.by = "Group",
  plot_type = "violin",
  add_box = TRUE,
  # 开启斑马纹背景
  add_bg = TRUE
)

# --- 3. 进阶效果：自定义背景色 + 统计 + 散点 ---
p <- Universal_FeatureStatPlot(
  data = df_test,
  value.var = "Value",
  group.by = "Group",
  plot_type = "box",
  add_point = TRUE,
  add_stat = "anova",
  # 开启背景并加深一点颜色
  add_bg = TRUE,
  bg_color = "grey95",
  title = "Zebra Striping Background Demo"
)

p + theme(panel.grid.major.x = element_blank())


library(ggplot2)
# 模拟 12 个分组的数据
df <- data.frame(
  Group = factor(rep(LETTERS[1:12], each = 10)),
  Value = rnorm(120)
)

# 使用 Paired 色板 (适合多组，最多12色)
Universal_FeatureStatPlot(
  data = df,
  value.var = "Value",
  group.by = "Group",
  plot_type = "box",
  add_bg = TRUE,
  add_trend = T,
  bg_color = "grey95",
  palette = "Paired" # <--- 直接写名字
)

# 使用 NPG 色板 (Nature 风格)
Universal_FeatureStatPlot(
  data = df,
  value.var = "Value",
  group.by = "Group",
  plot_type = "violin",
  add_box = T,
  palette = "NPG" # <--- 直接写名字
)


library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

#' 通用特征统计绘图函数 (实线趋势 + Box均值白点)
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "grey95",
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
  theme_use = theme_my_stat,
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

  if (!value.var %in% colnames(data)) {
    stop(paste0("Column '", value.var, "' not found."))
  }
  if (!group.by %in% colnames(data)) {
    stop(paste0("Column '", group.by, "' not found."))
  }
  if (!is.factor(data[[group.by]])) {
    data[[group.by]] <- factor(
      data[[group.by]],
      levels = unique(data[[group.by]])
    )
  }

  # --- 1. 颜色解析逻辑 ---
  n_groups <- nlevels(data[[group.by]])
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
      warning(paste("Palette", palette, "not found. Using NPG."))
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
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[group.by]])
  )

  # --- 3. 背景层 (斑马纹) ---
  if (isTRUE(add_bg)) {
    bg_df <- data.frame(x_pos = seq(1, n_groups, 2)) %>%
      mutate(xmin = x_pos - 0.5, xmax = x_pos + 0.5)

    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = bg_color,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    # 手动添加横向虚线
    y_range <- range(data[[value.var]], na.rm = TRUE)
    p <- p +
      geom_hline(
        yintercept = pretty(y_range, n = 5),
        linetype = "dashed",
        color = "white",
        size = 0.5
      )
  }

  # --- 4. 数据图层 ---

  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3
      )
  } else if (plot_type == "box") {
    # 主 Boxplot
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black"
      ) +
      # [新增] 均值小白点 (Shape 21)
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2.5,
        fill = "white",
        color = "black"
      )
  } else if (plot_type == "bar") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  }

  # --- 5. 叠加层 ---

  # [新增] 叠加 Box 时的均值小白点
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8
      ) +
      # 在叠加的 box 上也加上均值点
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2,
        fill = "white",
        color = "black"
      )
  }

  if (isTRUE(add_point)) {
    p <- p +
      geom_jitter(
        width = 0.2,
        size = 1.5,
        alpha = 0.6,
        color = "grey30",
        shape = 16
      )
  }

  # [修改] Trend 趋势线改为实线 (solid)
  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(group = 1),
        color = "black",
        linewidth = 1,
        linetype = "solid"
      ) # <--- 改为 solid
  }

  # --- 6. 统计层 ---
  if (add_stat %in% c("mean", "median")) {
    # 如果已经在 box 里加了点，这里可能会重复，但为了兼容 mean 模式，保留逻辑
    # 通常 box 自带点后，add_stat 建议设为 "none" 或检验方法
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    label_format <- if (is.null(comparisons) && n_groups > 2) {
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

  # --- 7. 收尾 ---
  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  if (exists("theme_my_stat")) {
    if (add_bg) {
      p <- p + theme_use(...) + theme(panel.grid = element_blank())
    } else {
      p <- p + theme_use(...)
    }
  } else {
    p <- p + theme_classic()
  }

  p <- p + scale_fill_manual(values = final_colors)
  p <- p + labs(title = title, x = xlab, y = ylab, fill = group.by)
  return(p)
}


library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

#' 通用特征统计绘图函数 (最终版：点在趋势线之上)
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "grey95",
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
  theme_use = theme_my_stat,
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

  if (!value.var %in% colnames(data)) {
    stop(paste0("Column '", value.var, "' not found."))
  }
  if (!group.by %in% colnames(data)) {
    stop(paste0("Column '", group.by, "' not found."))
  }
  if (!is.factor(data[[group.by]])) {
    data[[group.by]] <- factor(
      data[[group.by]],
      levels = unique(data[[group.by]])
    )
  }

  # --- 1. 颜色解析 ---
  n_groups <- nlevels(data[[group.by]])
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
      warning(paste("Palette", palette, "not found. Using NPG."))
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
    aes(x = .data[[group.by]], y = .data[[value.var]], fill = .data[[group.by]])
  )

  # --- 3. 背景层 ---
  if (isTRUE(add_bg)) {
    bg_df <- data.frame(x_pos = seq(1, n_groups, 2)) %>%
      mutate(xmin = x_pos - 0.5, xmax = x_pos + 0.5)

    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = bg_color,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    y_range <- range(data[[value.var]], na.rm = TRUE)
    p <- p +
      geom_hline(
        yintercept = pretty(y_range, n = 5),
        linetype = "dashed",
        color = "white",
        size = 0.5
      )
  }

  # --- 4. 主数据图层 ---
  if (plot_type == "violin") {
    p <- p +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.8,
        color = "black",
        size = 0.3
      )
  } else if (plot_type == "box") {
    p <- p +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = NA,
        color = "black"
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2.5,
        fill = "white",
        color = "black"
      )
  } else if (plot_type == "bar") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.7,
        alpha = 0.8,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else if (plot_type == "dot") {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 21,
        color = "black"
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  }

  # --- 5. 叠加层 (调整了顺序) ---

  # [顺序调整 1] 先画 Trend (这样它在最下面)
  if (isTRUE(add_trend)) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(group = 1),
        color = "black",
        linewidth = 1,
        linetype = "solid"
      )
  }

  # [顺序调整 2] 再画叠加 Box (如果在 Trend 上面)
  if (isTRUE(add_box) && plot_type == "violin") {
    p <- p +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA,
        color = "black",
        alpha = 0.8
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2,
        fill = "white",
        color = "black"
      )
  }

  # [顺序调整 3] 最后画 Point (这样点就在 Trend 之上)
  if (isTRUE(add_point)) {
    p <- p +
      geom_jitter(
        width = 0.2,
        size = 1.5,
        alpha = 0.6,
        color = "grey30",
        shape = 16
      )
  }

  # --- 6. 统计层 ---
  if (add_stat %in% c("mean", "median")) {
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    label_format <- if (is.null(comparisons) && n_groups > 2) {
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

  # --- 7. 收尾 ---
  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  if (exists("theme_my_stat")) {
    if (add_bg) {
      p <- p + theme_use(...) + theme(panel.grid = element_blank())
    } else {
      p <- p + theme_use(...)
    }
  } else {
    p <- p + theme_classic()
  }

  p <- p + scale_fill_manual(values = final_colors)
  p <- p + labs(title = title, x = xlab, y = ylab, fill = group.by)
  return(p)
}


Universal_FeatureStatPlot(
  data = df,
  value.var = "Value",
  group.by = "Group",
  plot_type = "violin",
  add_box = T,
  comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")),
  palette = "NPG" # <--- 直接写名字
)


library(ggplot2)
library(ggpubr)
library(dplyr)

# 设置随机种子以便复现
set.seed(2024)

# 构造数据
df_split_demo <- data.frame(
  # X 轴分组：对照组 vs 治疗组
  Group = rep(c("Control", "Treat"), times = 40),

  # 分面分组：雄性 vs 雌性
  Gender = rep(c("Male", "Female"), each = 40),

  # 数值：故意制造差异
  # Male: Treat 比 Control 高很多 (15 vs 10)
  # Female: Treat 和 Control 差不多 (10.5 vs 10)
  Value = c(
    rnorm(20, mean = 10, sd = 1),
    rnorm(20, mean = 15, sd = 1), # Male 的数据
    rnorm(20, mean = 10, sd = 1),
    rnorm(20, mean = 10.5, sd = 1) # Female 的数据
  )
)

# 确保因子顺序 (Control 在前)
df_split_demo$Group <- factor(
  df_split_demo$Group,
  levels = c("Control", "Treat")
)


Universal_FeatureStatPlot(
  data = df_split_demo,
  value.var = "Value",
  group.by = "Group", # X轴：比较 Control vs Treat
  split.by = "Gender", # 分面：按性别分开画

  # --- 统计设置 ---
  add_stat = "t.test", # 使用 T 检验
  # 指定比较对象：比较 X 轴上的 "Control" 和 "Treat"
  comparisons = list(c("Control", "Treat")),

  # --- 绘图风格 ---
  plot_type = "box", # 箱线图
  add_point = TRUE, # 加散点
  add_bg = F, # 加斑马纹背景
  palette = "NPG", # 使用 NPG 配色
  title = "Independent T-test by Gender"
)


library(ggplot2)
library(ggpubr)

# 1. 构造模拟数据 (3个组)
set.seed(2024)
df_one_way <- data.frame(
  Group = factor(
    rep(c("Control", "Treat1", "Treat2"), each = 30),
    levels = c("Control", "Treat1", "Treat2")
  ),
  Value = c(rnorm(30, 10), rnorm(30, 12), rnorm(30, 15)) # 故意制造差异
)

# 2. 进行组间比较 (不使用 split.by)
p <- Universal_FeatureStatPlot(
  data = df_one_way,
  value.var = "Value",
  group.by = "Group",

  # --- 核心设置 ---
  add_stat = "t.test", # 选择统计方法
  comparisons = list(
    # 指定比较对
    c("Control", "Treat1"),
    c("Control", "Treat2"),
    c("Treat1", "Treat2")
  ),

  # --- 绘图风格 ---
  plot_type = "violin",
  add_box = TRUE,
  add_point = TRUE,
  add_bg = TRUE, # 斑马纹背景
  palette = "NPG",
  theme_use = theme_bw,
  title = "Pairwise Comparisons (No Split)"
)


library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

#' 通用特征统计绘图函数 (纯净斑马纹背景版)
#' @export
Universal_FeatureStatPlot <- function(
  data,
  value.var,
  group.by,
  fill.by = NULL,
  split.by = NULL,
  add_bg = FALSE,
  bg_color = "#00000008",
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
  theme_use = theme_my_stat,
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
    stop("Column not found.")
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

  # --- 3. [纯净背景层] 只画斑马纹 ---
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

  # --- 4. 数据图层 ---
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
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2.5,
        fill = "white",
        color = "black",
        position = pos
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
      ) # 实线趋势
    if (!is_grouped) p <- p + aes(color = NULL)
  }

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
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2,
        fill = "white",
        color = "black",
        position = pos
      ) # 均值白点
  }

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
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black",
        position = pos
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

  if (exists("theme_my_stat")) {
    if (add_bg) {
      # 如果加了斑马纹，为了美观，我们通常会把主题自带的纵向网格也去掉
      # 只保留斑马纹作为区分，非常干净
      p <- p + theme_use(...) + theme(panel.grid = element_blank())
    } else {
      p <- p + theme_use(...)
    }
  } else {
    p <- p + theme_classic()
  }

  p <- p + scale_fill_manual(values = final_colors)
  p <- p + labs(title = title, x = xlab, y = ylab, fill = fill_var)

  if (is_grouped && add_trend) {
    p <- p + scale_color_manual(values = rep("black", n_groups), guide = "none")
  }

  return(p)
}


source('E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\my_themes.R')
p + fix_panel(5, 4)


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
  bg_color = "#0000000D", # 默认 5% 透明度黑色
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
  theme_use = theme_classic, # 默认主题
  ...
) {
  # --- 0. 参数校验与预处理 ---
  plot_type <- match.arg(plot_type)
  if (is.null(xlab)) {
    xlab <- group.by
  }
  if (is.null(ylab)) {
    ylab <- value.var
  }

  # 确定填充变量
  fill_var <- if (is.null(fill.by)) group.by else fill.by

  # 检查列是否存在
  cols_needed <- c(value.var, group.by, fill_var)
  if (!all(cols_needed %in% colnames(data))) {
    stop("指定的列名在数据框中不存在，请检查拼写。")
  }

  # 转换为因子 (Factor) 以保证绘图顺序
  for (col in unique(c(group.by, fill_var))) {
    if (!is.factor(data[[col]])) {
      data[[col]] <- factor(data[[col]], levels = unique(data[[col]]))
    }
  }

  # --- 1. 智能颜色解析 ---
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

  # 决定颜色向量
  if (length(palette) == 1 && is.character(palette)) {
    if (palette %in% names(sci_palettes)) {
      base_cols <- sci_palettes[[palette]]
    } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      base_cols <- RColorBrewer::brewer.pal(min(n_groups, max_n), palette)
    } else {
      warning(paste("未找到色板", palette, "。已回退到 NPG 配色。"))
      base_cols <- sci_palettes[["NPG"]]
    }
  } else {
    base_cols <- palette # 用户自定义向量
  }

  # 自动插值 (如果分组数 > 颜色数)
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

  # 判断是否为 Grouped (Dodge) 模式
  is_grouped <- (fill_var != group.by)
  dodge_width <- 0.8

  # --- 3. 背景层 (斑马纹) ---
  if (isTRUE(add_bg)) {
    n_x_groups <- nlevels(data[[group.by]])
    # 生成奇数位置的背景矩形
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
  # 设置 Dodge 参数
  pos <- if (is_grouped) {
    position_dodge(width = dodge_width)
  } else {
    position_dodge(width = 0)
  }

  # 绘制主图
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
      # Box 自带均值白点
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2.5,
        fill = "white",
        color = "black",
        position = pos
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

  # --- 5. 叠加层 (注意绘制顺序: Trend -> Box -> Point) ---

  # (A) 趋势线 (实线)
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

  # (B) 叠加 Box (带均值白点)
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
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 21,
        size = 2,
        fill = "white",
        color = "black",
        position = pos
      )
  }

  # (C) 散点 (最顶层)
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
    p <- p +
      stat_summary(
        fun = add_stat,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black",
        position = pos
      )
  }
  if (add_stat %in% c("t.test", "wilcox.test", "anova", "kruskal.test")) {
    if (is_grouped) {
      # 组内比较
      p <- p +
        stat_compare_means(
          aes(group = .data[[fill_var]]),
          method = add_stat,
          label = "p.signif"
        )
    } else {
      # 组间比较
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

  # --- 7. 收尾 (主题与配色) ---
  if (!is.null(split.by)) {
    p <- p + facet_wrap(as.formula(paste("~", split.by)), scales = "free")
  }

  # 应用主题函数
  if (is.function(theme_use)) {
    p <- p + theme_use(...)
  } else {
    warning("theme_use 不是一个有效的函数，使用默认 theme_classic")
    p <- p + theme_classic(...)
  }

  # 如果开启了斑马纹背景，强制去除面板网格线，保持干净
  if (add_bg) {
    p <- p + theme(panel.grid = element_blank())
  }

  # 应用颜色
  p <- p + scale_fill_manual(values = final_colors)
  p <- p + labs(title = title, x = xlab, y = ylab, fill = fill_var)

  # 修正 grouped trend 可能产生的额外图例
  if (is_grouped && add_trend) {
    p <- p + scale_color_manual(values = rep("black", n_groups), guide = "none")
  }

  return(p)
}


df <- data.frame(
  Group = factor(rep(c("Ctrl", "Treat"), each = 20)),
  Value = rnorm(40)
)

compare_plot(
  data = df,
  value.var = "Value",
  group.by = "Group",
  plot_type = "box",
  add_point = TRUE,
  add_stat = "t.test",
  add_bg = TRUE, # 开启斑马纹
  theme_use = theme_classic
)

compare_plot(
  data = df,
  value.var = "Value",
  group.by = "Group",
  theme_use = theme_my_stat, # <--- 传入 theme_bw
  base_size = 14 # <--- 传递参数调整字号
)


library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

#' 组间比较可视化函数 (已修复组内比较时的均值点错位问题)
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


df_complex <- data.frame(
  Cluster = factor(rep(c("C1", "C2", "C3"), each = 20)),
  Treatment = factor(rep(c("Ctrl", "Exp"), times = 30)),
  Expression = rnorm(60)
)

compare_plot(
  data = df_complex,
  value.var = "Expression",
  group.by = "Cluster", # X轴
  fill.by = "Treatment", # 填充色 (实现组内并排)
  plot_type = "box",
  add_stat = "wilcox.test", # 自动计算组内差异
  palette = "AAAS", # 使用杂志配色
  title = "Intra-cluster Treatment Comparison"
)


df_multi <- data.frame(Group = factor(rep(1:4, each = 10)), Value = rnorm(40))

compare_plot(
  data = df_multi,
  value.var = "Value",
  group.by = "Group",
  add_stat = "anova", # 全局 ANOVA
  plot_type = "violin",
  add_box = TRUE,
  add_bg = TRUE
)


library(ggplot2)
library(dplyr)

set.seed(2026)
# 模拟单细胞降维与表达数据
sc_data <- data.frame(
  UMAP_1 = rnorm(500),
  UMAP_2 = rnorm(500),
  CellType = sample(
    c("B cells", "T cells", "Macrophages", "Fibroblasts", "Endothelial"),
    500,
    replace = TRUE
  ),
  Treatment = sample(c("Ctrl", "Treated"), 500, replace = TRUE),
  Gene_Expr = runif(500, min = 0, max = 10) # 连续表达量
)

# 模拟空间转录组的网格热图数据
spatial_grid <- expand.grid(X = 1:15, Y = 1:15)
spatial_grid$Enrichment_Score <- rnorm(225)

# --- 定义全局字典 ---
# 严格命名向量 (锁定细胞颜色)
set_palette(
  "MyCells",
  c(
    "B cells" = "#3C5488",
    "T cells" = "#E64B35",
    "Macrophages" = "#00A087",
    "Fibroblasts" = "#F39B7F",
    "Endothelial" = "#4DBBD5"
  )
)

# 故意定义一个只有 3 种颜色的短色板 (用于测试自动插值)
set_palette("ShortPal", c("#FF0000", "#00FF00", "#0000FF"))


# 故意剔除两种细胞，制造缺失情况
sub_data <- sc_data %>%
  filter(CellType %in% c("B cells", "T cells", "Fibroblasts"))

ggplot(sub_data, aes(x = CellType, y = Gene_Expr, fill = CellType)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  # 调用全局字典
  scale_fill_univ_d("MyCells", alpha = 0.9) +
  labs(
    title = "Test 1: Strict Mapping & Subset Immunity",
    subtitle = "Colors remain locked to CellType despite missing categories."
  )

# 验证效果：虽然只有3个箱子，但 B cells 依然是深蓝，T cells 依然是红色。没有发生颜色顺延！

ggplot(sc_data, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 2) +
  theme_minimal() +
  # 5 个细胞组，但硬塞给它一个只有 3 色的字典
  scale_color_univ_d(palette = "NineteenEightyR::miami2") +
  labs(
    title = "Test 2: Auto-Interpolation for Discrete Scale",
    subtitle = "5 groups flawlessly mapped using a 3-color palette."
  )

# 验证效果：不会报错。系统会自动把 红-绿-蓝 平滑切分成 5 种离散颜色分配给点。

ggplot(sc_data, aes(x = UMAP_1, y = UMAP_2, color = Gene_Expr)) +
  geom_point(size = 3) +
  # theme_dark() + # 暗色背景更容易看出高亮颜色
  # 调用内置 Viridis (适合表达量)，并设置透明度
  scale_color_univ_c(palette = "viridis::turbo") +
  labs(
    title = "Test 3: Continuous Color Mapping with Alpha",
    subtitle = "Viridis gradient with 60% opacity."
  )

# 验证效果：完美的渐变色点图，点与点重叠的地方透明度叠加自然，没有变成黑块。

ggplot(spatial_grid, aes(x = X, y = Y, fill = Enrichment_Score)) +
  geom_tile(color = "white") +
  theme_void() +
  # 调用内置 RWB，并反转颜色顺序 (变成 Blue-White-Red)
  scale_fill_univ_c("RWB", reverse = TRUE) +
  labs(
    title = "Test 4: Continuous Fill Reversal",
    subtitle = "Spatial grid using reversed RWB (Red is high)."
  )

# 验证效果：生成一张网格热图。正值（高分）映射为深红色，负值映射为深蓝色，0 附近为白色。

get_colors("ShortPal")


# 1. 在一个专门的地方定义好你的超级大字典 (List)
my_project_master_palettes <- list(
  # 肿瘤微环境细胞专属画板
  "My_TME_Cells" = c(
    "STC1+ CAFs" = "#DC0000",
    "Normal Fibroblasts" = "#4DBBD5",
    "Macrophages" = "#00A087"
  ),

  # 组织取样来源画板
  "Tissue_Source" = c(
    "Tumor_Core" = "#E64B35",
    "Invasive_Margin" = "#F39B7F",
    "Normal_Adjacent" = "#3C5488"
  ),

  # 连续型测序深度渐变色板 (不需要名字的纯向量)
  "Seq_Depth_Grad" = c("#e0f3db", "#a8ddb5", "#43a2ca")
)

# 2. 一键批量打入系统！
set_palette(my_project_master_palettes)


source('配色系统调试.R')

else
if (
  requireNamespace("paletteer", quietly = TRUE) &&
    grepl("::", palette)
) {
  tryCatch(
    {
      # 第 1 步：优先尝试作为离散型画板提取 (绝大多数情况)
      cols <- as.character(paletteer::paletteer_d(palette))
    },
    error = function(e1) {
      # 第 2 步：【核心升级】如果离散型报错，尝试作为连续型画板提取
      tryCatch(
        {
          # 提取连续渐变带上的 100 个颜色节点作为基准色池
          cols <- as.character(paletteer::paletteer_c(palette, n = 100))
        },
        error = function(e2) {
          # 如果两边都找不到，再弹出警告
          warning(sprintf(
            "无法从 paletteer 加载画板 '%s'，请检查拼写。",
            palette
          ))
          cols <<- NULL # 使用 <<- 赋值给外层环境
        }
      )
    }
  )
}


install.packages("gggenes")
library(gggenes)

library(ggplot2)
library(gggenes)
data(example_genes)

ggplot(
  example_genes,
  aes(xmin = start, xmax = end, y = molecule, fill = gene)
) +
  geom_gene_arrow() +
  facet_wrap(~molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")


# 确保已安装并加载所需包
library(tidyverse)
library(scales)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(eulerr)
library(grid)

plot_circos_heatmap <- function(
  data,
  group_col = "Group", # 分组列名
  feature_col = "Gene", # 特征列名(如Gene, OTU, 代谢物等)
  value_cols, # 表达量/丰度列名向量
  qval_col = "qval", # 显著性列名
  group_levels = NULL, # 分组的顺序
  center_title = "Shared\nFeatures" # 中心韦恩图下方的文本
) {
  # 1. 处理分组因子级别
  if (!is.null(group_levels)) {
    data[[group_col]] <- factor(data[[group_col]], levels = group_levels)
  } else {
    data[[group_col]] <- as.factor(data[[group_col]])
    group_levels <- levels(data[[group_col]])
  }

  # 2. 动态生成调色板 (可以根据需要提取到函数参数中让用户自定义)
  green_pink <- colorRamp2(
    seq(0, 0.5, length.out = 100),
    colorRampPalette(rev(brewer.pal(n = 5, name = "PiYG")))(100)
  )
  red_black <- colorRamp2(
    seq(0, 1, length.out = 6),
    c("#fff1e7", "#f99c70", "#dd4139", "#97044d", "#4f2043", "black")
  )

  # 为动态数量的分组生成颜色
  group_colors <- brewer.pal(max(3, length(group_levels)), "Set2")[
    1:length(group_levels)
  ]
  names(group_colors) <- group_levels

  # 3. 数据预处理与标准化
  data_normalized <- data %>%
    group_by(!!sym(group_col)) %>%
    arrange(!!sym(qval_col), .by_group = TRUE) %>%
    mutate(across(all_of(value_cols), ~ rescale(., to = c(0, 0.5)))) %>%
    ungroup()

  # 转为矩阵
  data_matrix <- data_normalized %>%
    select(all_of(value_cols)) %>%
    as.matrix()
  rownames(data_matrix) <- data_normalized[[feature_col]]

  # 4. 初始化 Circos 图
  circos.clear()
  circos.par(
    start.degree = 60,
    gap.after = c(rep(5, length(group_levels) - 1), 30),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  # 5. 绘制核心热图
  circos.heatmap(
    data_matrix,
    split = data_normalized[[group_col]],
    cluster = FALSE,
    bg.border = "black",
    bg.lwd = 1,
    cell.border = "white",
    cell.lwd = 0.5,
    rownames.side = "outside",
    rownames.cex = 1.2,
    col = green_pink,
    track.height = 0.25
  )

  # 6. 添加列名轨道
  circos.track(
    track.index = get.current.track.index(),
    bg.border = NA,
    panel.fun = function(x, y) {
      if (CELL_META$sector.numeric.index == length(group_levels)) {
        cn <- colnames(data_matrix)
        n <- length(cn)
        cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
        y_coords <- seq(
          CELL_META$cell.ylim[1] + cell_height / 2,
          CELL_META$cell.ylim[2] - cell_height / 2,
          length.out = n
        )
        for (i in 1:n) {
          circos.lines(
            c(
              CELL_META$cell.xlim[2],
              CELL_META$cell.xlim[2] + convert_x(1, "mm")
            ),
            c(y_coords[i], y_coords[i]),
            col = "black",
            lwd = 2
          )
        }
        circos.text(
          rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
          y_coords,
          cn,
          cex = 1.4,
          adj = c(0, 0.5),
          facing = "inside"
        )
      }
    }
  )

  # 7. 添加 qval 散点图轨道
  circos.track(
    ylim = c(0, 1),
    track.height = 0.05,
    bg.border = NA,
    panel.fun = function(x, y) {
      # 动态过滤当前扇区的数据
      sector_data <- data_normalized[
        data_normalized[[group_col]] == CELL_META$sector.index,
      ]

      for (i in 1:nrow(sector_data)) {
        circos.points(
          CELL_META$xlim[1] +
            (CELL_META$xlim[2] - CELL_META$xlim[1]) *
              (i - 0.5) /
              nrow(sector_data),
          0.5,
          pch = 18,
          cex = 1.8,
          col = red_black(sector_data[[qval_col]][i] * 100)
        )
      }

      if (CELL_META$sector.numeric.index == length(group_levels)) {
        circos.lines(
          c(
            CELL_META$cell.xlim[2],
            CELL_META$cell.xlim[2] + convert_x(1, "mm")
          ),
          c(0.5, 0.5),
          col = "black",
          lwd = 2
        )
        circos.text(
          CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
          0.5,
          "qval",
          cex = 1.4,
          adj = c(0, 0.5),
          facing = "inside"
        )
      }
    }
  )

  # 8. 添加分组标签轨道
  circos.track(
    ylim = c(0, 1),
    track.height = 0.065,
    bg.col = adjustcolor(
      group_colors[levels(data_normalized[[group_col]])],
      alpha.f = 0.3
    ),
    panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[2] - 0.75,
        CELL_META$sector.index,
        facing = "bending.inside",
        cex = 1.5,
        adj = c(0.5, 0)
      )
    }
  )

  # 9. 添加图例
  heatmap_legend <- Legend(
    title = "Expression\nfraction",
    col_fun = green_pink,
    at = seq(0, 0.5, length.out = 6),
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 14),
    labels_gp = gpar(fontsize = 14)
  )
  qval_legend <- Legend(
    title = "qval",
    col_fun = red_black,
    at = seq(0, 1, length.out = 6),
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 14),
    labels_gp = gpar(fontsize = 14)
  )

  draw(
    heatmap_legend,
    x = unit(0.95, "npc") - unit(5, "mm"),
    y = unit(0.9, "npc") - unit(5, "mm"),
    just = c("right", "top")
  )
  draw(
    qval_legend,
    x = unit(0.93, "npc") - unit(5, "mm"),
    y = unit(0.25, "npc") - unit(5, "mm"),
    just = c("right", "top")
  )

  # 10. 动态生成韦恩图数据并插入中心
  # 使用 split 自动根据组别提取特征列表，无需再写死基因集提取逻辑
  feature_groups <- split(data[[feature_col]], data[[group_col]])
  fit <- euler(feature_groups)

  pushViewport(viewport(width = 0.24, height = 0.24))
  circos.track(
    ylim = c(0, 1),
    track.height = 0.3,
    bg.border = NA,
    panel.fun = function(x, y) {
      grid.draw(plot(
        fit,
        fills = list(fill = group_colors, alpha = 0.2), # 颜色与组别对应
        edges = list(col = "black", lty = 2, lwd = 2),
        labels = FALSE,
        quantities = TRUE
      ))
      grid.text(
        center_title,
        x = unit(0.5, "npc"),
        y = unit(-0.01, "npc"),
        gp = gpar(fontsize = 25, lineheight = 0.8)
      )
    }
  )
  popViewport() # 养成好习惯，画完退出视口
}


library(dplyr)

# 设置随机种子以保证结果可重复
set.seed(2026)

# 定义一些特征基因
target_genes <- c(
  "STC1",
  "PA28gamma",
  "C1QBP",
  "FAP",
  "ACTA2",
  "COL1A1",
  "MMP9",
  "VEGFA",
  "CXCL12",
  "IL6",
  "TGFB1",
  "CD3E",
  "CD8A",
  "FOXP3",
  "CD68",
  "CD163",
  "PDCD1",
  "CTLA4",
  "MKI67",
  "TP53",
  "EGFR",
  "MYC",
  "SNAI1",
  "TWIST1",
  "VIM"
)

# 定义3个感兴趣的分组
tme_groups <- c("STC1_High_CAFs", "STC1_Low_CAFs", "Normal_Fibroblasts")

# 构造初始数据框 (随机抽取生成60行数据)
simulated_data <- data.frame(
  Gene = sample(target_genes, 60, replace = TRUE),
  Group = sample(tme_groups, 60, replace = TRUE),
  # 模拟5个不同空间位置的表达量分布 (对应原图的5个数据列)
  Tumor_Core = runif(60, min = 10, max = 500),
  Invasive_Margin = runif(60, min = 10, max = 500),
  Stroma = runif(60, min = 10, max = 500),
  Adjacent_Normal = runif(60, min = 0, max = 200),
  Lymph_Node = runif(60, min = 0, max = 300),
  # 模拟显著性 q-value (使用 beta 分布让数据偏向于较小的显著数值)
  qval = rbeta(60, shape1 = 1, shape2 = 5)
)

# 去重：确保同一个基因在同一个Group中只出现一次，避免热图矩阵行名冲突
test_df <- simulated_data %>%
  distinct(Gene, Group, .keep_all = TRUE)

# 查看生成的数据结构
head(test_df)


# 运行绘图函数
plot_circos_heatmap(
  data = test_df,
  group_col = "Group",
  feature_col = "Gene",
  # 传入我们模拟的5个空间位置列名
  value_cols = c(
    "Tumor_Core",
    "Invasive_Margin",
    "Stroma",
    "Adjacent_Normal",
    "Lymph_Node"
  ),
  qval_col = "qval",
  # 指定分组在环图上的逆时针排列顺序
  group_levels = c("STC1_High_CAFs", "STC1_Low_CAFs", "Normal_Fibroblasts"),
  center_title = "Shared\nMarkers"
)


# 1. 模拟数据 (增加 logFC 列)
set.seed(42)
genes <- paste0("Gene_", 1:45)
groups <- c(
  rep("STC1+ CAFs", 15),
  rep("Normal Fibroblasts", 15),
  rep("Other TME Cells", 15)
)

test_data <- data.frame(
  Gene = genes,
  Group = groups,
  Sample1 = runif(45, 10, 100),
  Sample2 = runif(45, 10, 100),
  Sample3 = runif(45, 10, 100),
  Sample4 = runif(45, 10, 100),
  qval = rbeta(45, 1, 5), # 偏小的 p/q 值
  logFC = rnorm(45, mean = 0, sd = 2) # 有正有负的 logFC
)

# 2. 定义指定对象和具体颜色的字典 (直接接收您定义的 list)
my_tme_palette <- c(
  "STC1+ CAFs" = "#DC0000",
  "Normal Fibroblasts" = "#4DBBD5",
  "Other TME Cells" = "#00A087"
)

# 3. 定义动态轨道的配置
# 可以无限往列表中增加元素，函数会自动由内向外叠放轨道
my_tracks <- list(
  # 轨道1: 散点图展示显著性
  list(
    column = "qval",
    type = "points",
    label = "q-value",
    height = 0.06,
    color = colorRamp2(
      seq(0, 1, length.out = 3),
      c("#dd4139", "grey50", "black")
    )
  ),

  # 轨道2: 柱状图展示 logFC
  list(
    column = "logFC",
    type = "bars",
    label = "log2FC",
    height = 0.1,
    color = "#3C5488"
  )
)

# 4. 执行绘图
plot_multi_track_circos(
  data = test_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = c("Sample1", "Sample2", "Sample3", "Sample4"),
  group_palette = my_tme_palette, # 传入自定义颜色板
  track_configs = my_tracks # 传入轨道配置列表
)


# 1. 模拟更大规模的数据
set.seed(2026)
genes <- paste0("Gene_", 1:50) # 50个基因
groups <- c(
  rep("STC1+ CAFs", 17),
  rep("Normal Fibroblasts", 17),
  rep("Other TME Cells", 16)
)

# 动态生成 20 个样本的列名
sample_cols <- paste0("Sample_", 1:20)

# 构造基础数据框
test_data <- data.frame(
  Gene = genes,
  Group = groups
)

# 批量添加 20 个样本的表达量数据 (随机生成 10 到 100 的丰度)
for (col in sample_cols) {
  test_data[[col]] <- runif(50, 10, 100)
}

# 添加 4 种额外的轨道维度数据
test_data$qval <- rbeta(50, 1, 5) # 轨道1：显著性 q-value (0~1)
test_data$logFC <- rnorm(50, mean = 0, sd = 2) # 轨道2：差异倍数 log2FC (有正有负)
test_data$CNV_score <- rnorm(50, mean = 5, sd = 2) # 轨道3：拷贝数变异打分 (连续数值)
test_data$Methylation <- runif(50, 0, 1) # 轨道4：甲基化比例 (0~1)


# 2. 预设分组的调色板
my_tme_palette <- c(
  "STC1+ CAFs" = "#DC0000",
  "Normal Fibroblasts" = "#4DBBD5",
  "Other TME Cells" = "#00A087"
)


# 3. 配置 4 层外围轨道 (高度加起来约为 0.23，加上热图的 0.2，完全在安全范围内)
my_tracks <- list(
  # 轨道 1: q-value 散点图 (越接近 0 颜色越红)
  list(
    column = "qval",
    type = "points",
    label = "q-value",
    height = 0.05,
    color = colorRamp2(
      seq(0, 1, length.out = 3),
      c("#dd4139", "grey50", "black")
    )
  ),

  # 轨道 2: logFC 柱状图 (✨ 核心技巧：利用匿名函数实现条件着色，大于0为红，小于0为蓝)
  list(
    column = "logFC",
    type = "bars",
    label = "log2FC",
    height = 0.08,
    color = function(x) ifelse(x > 0, "#DC0000", "#3C5488")
  ),

  # 轨道 3: CNV 变异折线图
  list(
    column = "CNV_score",
    type = "lines",
    label = "CNV",
    height = 0.05,
    color = "#F39B7F"
  ),

  # 轨道 4: 甲基化程度 散点图 (从白到紫的渐变映射)
  list(
    column = "Methylation",
    type = "points",
    label = "Methyl",
    height = 0.05,
    color = colorRamp2(c(0, 1), c("white", "#4f2043"))
  )
)

# 4. 执行绘图
plot_multi_track_circos(
  data = test_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = sample_cols, # 直接传入包含 20 个样本名称的向量
  group_palette = my_tme_palette, # 传入自定义颜色板
  track_configs = my_tracks # 传入 4 层轨道配置
)


# 制作一个经典的“蓝-白-红” Z-score 风格渐变色带
my_bwr_palette <- colorRamp2(
  seq(0, 0.5, length.out = 100),
  colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100)
)

# 画图时传入
plot_multi_track_circos(
  data = test_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = sample_cols,
  group_palette = my_tme_palette,
  heatmap_col = my_bwr_palette, # ✨ 在这里传入你定义的色带
  track_configs = my_tracks
)

source('环状热图.R')
# 1. 模拟含极端异常值的数据
set.seed(42)
test_data <- data.frame(
  Gene = paste0("Gene_", 1:45),
  Group = c(
    rep("STC1+ CAFs", 15),
    rep("Normal Fibroblasts", 15),
    rep("Other", 15)
  ),
  Sample1 = runif(45, 10, 100),
  Sample2 = runif(45, 10, 100),
  Sample3 = runif(45, 10, 100),
  # 模拟原始丰度：第10个基因故意设为 5000（极值）
  Raw_Abundance = c(runif(9, 1, 50), 5000, runif(35, 1, 50)),
  # 模拟免疫得分：数值起伏大
  Immune_Score = rnorm(45, mean = 10, sd = 15)
)

my_palette <- c(
  "STC1+ CAFs" = "#DC0000",
  "Normal Fibroblasts" = "#4DBBD5",
  "Other" = "#00A087"
)

# 2. 轨道配置：激活限制引擎
my_tracks <- list(
  # 轨道1：棒棒糖图展示原始丰度，启用盖帽引擎
  # 若不加 cap_outliers = TRUE，其他 44 个基因会被压缩成一条贴底的平线
  list(
    column = "Raw_Abundance",
    type = "lollipop",
    label = "Abundance",
    height = 0.1,
    color = "#E64B35",
    cap_outliers = TRUE
  ),

  # 轨道2：面积填充图展示免疫得分，同时执行 Z-score 标准化
  list(
    column = "Immune_Score",
    type = "area",
    label = "Immune\nZ-score",
    height = 0.08,
    color = "#3C5488",
    scale = "zscore"
  )
)

# 3. 运行函数：关闭行名和列名展示
plot_advanced_circos(
  data = test_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = c("Sample1", "Sample2", "Sample3"),
  group_palette = my_palette,
  track_configs = my_tracks,
  heatmap_col = my_bwr_palette,
  show_rownames = FALSE, # ✨ 隐藏基因名，让内圈更干净
  show_colnames = FALSE # ✨ 隐藏样本列名
)


# 1. 构造基础数据
set.seed(123)
genes <- paste0("Gene_", 1:60)
groups <- c(
  rep("STC1+ CAFs", 20),
  rep("Normal Fibroblasts", 20),
  rep("Other TME Cells", 20)
)
sample_cols <- paste0("Sample_", 1:10)

test_data <- data.frame(Gene = genes, Group = groups)
for (col in sample_cols) {
  test_data[[col]] <- runif(60, 0, 100)
}
test_data$logFC <- rnorm(60, 0, 2)

# 2. 准备轨道配置和调色板
my_palette <- c(
  "STC1+ CAFs" = "#DC0000",
  "Normal Fibroblasts" = "#4DBBD5",
  "Other TME Cells" = "#00A087"
)
my_tracks <- list(
  list(
    column = "logFC",
    type = "bars",
    label = "logFC",
    height = 0.1,
    color = function(x) ifelse(x > 0, "#DC0000", "#3C5488")
  )
)

# ✨ 3. 构造连线网络数据 (From -> To)
# 模拟 STC1+ CAFs 细胞中的基因，去调控其他细胞群的基因
interaction_links <- data.frame(
  from = c("Gene_2", "Gene_5", "Gene_8", "Gene_15", "Gene_15"), # 起点：大部分都在 CAFs 扇区
  to = c("Gene_50", "Gene_25", "Gene_45", "Gene_55", "Gene_30"), # 终点：指向其他扇区
  # 为不同的连线赋予不同的颜色（比如依据互作强度或通路分配颜色）
  color = c("#DC000080", "#00A08780", "#DC000080", "#4DBBD580", "#3C548880")
)

# 4. 运行绘图
# 注意：要把前面补全好的完整 plot_advanced_circos 函数先运行加载进内存
plot_advanced_circos(
  data = test_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = sample_cols,
  group_palette = my_palette,
  track_configs = my_tracks,
  links_data = interaction_links, # ✨ 将网络关系网传入
  show_rownames = FALSE,
  show_colnames = FALSE
)


# ==========================================
# 1. 构造多维度模拟数据 (TME 场景)
# ==========================================
set.seed(2026)
genes <- paste0("Gene_", 1:60)
# 设定三个经典的 TME 细胞类群
groups <- c(
  rep("STC1+ CAFs", 20),
  rep("Normal Fibroblasts", 20),
  rep("Tumor Cells", 20)
)
sample_cols <- paste0("Sample_", 1:5)

test_data <- data.frame(Gene = genes, Group = groups)

# 模拟 5 个样本的表达量矩阵 (用于最内圈热图)
for (col in sample_cols) {
  test_data[[col]] <- runif(60, 0, 100)
}

# 模拟附加轨道的数据
# 轨道数据 A: 差异倍数 (基因水平，每个基因一个值)
test_data$logFC <- rnorm(60, mean = 0, sd = 2)

# 轨道数据 B: 总体表达分布特征 (用于绘制箱线图)
# 我们让 STC1+ CAFs 的整体表达均值偏高，Normal Fibroblasts 偏低
test_data$Expr_Dist <- c(rnorm(20, 80, 10), rnorm(20, 30, 8), rnorm(20, 50, 15))


# ==========================================
# 2. 预设调色板与网络连线数据
# ==========================================
# 严格按照您预先定义的列表，为指定的细胞类群分配精准的颜色
My_TME_Cells <- c(
  "STC1+ CAFs" = "#DC0000",
  "Normal Fibroblasts" = "#4DBBD5",
  "Tumor Cells" = "#00A087"
)

# 模拟细胞间通讯网络 (从 STC1+ CAFs 的靶点指向 Tumor Cells)
interaction_links <- data.frame(
  from = c("Gene_3", "Gene_8", "Gene_15", "Gene_12"),
  to = c("Gene_45", "Gene_52", "Gene_45", "Gene_58"),
  color = c("#DC000080", "#DC000080", "#4DBBD580", "#00A08780") # 带透明度的线条
)


# ==========================================
# 3. 配置外围动态轨道
# ==========================================
my_tracks <- list(
  # 第 1 圈 (靠内): logFC 差异倍数柱状图 (动态着色：大于0标红，小于0标蓝)
  list(
    column = "logFC",
    type = "bars",
    label = "logFC",
    height = 0.08,
    color = function(x) ifelse(x > 0, "#DC0000", "#3C5488")
  ),

  # 第 2 圈 (靠外): 箱线图 (展示该类群整体的 Expr_Dist 分布特征)
  # 注意：由于它是箱线图，每个扇区中心只会生成一个代表该组数据的“箱子”
  list(
    column = "Expr_Dist",
    type = "boxplot",
    label = "Dist.\n(Box)",
    height = 0.12, # 箱线图建议给足高度
    color = function(x) "#F39B7F"
  ) # 给箱体统一一个柔和的橙色填充
)


# ==========================================
# 4. 执行绘图
# ==========================================
plot_advanced_circos(
  data = test_data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = sample_cols,
  group_palette = My_TME_Cells, # 传入您指定的颜色字典
  track_configs = my_tracks, # 传入轨道配置 (包含 boxplot)
  links_data = interaction_links, # 传入跨圆心连线数据
  show_rownames = FALSE, # 隐藏基因名防重叠
  show_colnames = TRUE, # 保留热图右侧的样本标签
  gap_degree = 90,
  start_degree = 270
)

source('环状热图.R')


# 1. 构造 4 个分组的模拟数据
set.seed(2026)
genes <- paste0("Gene_", 1:80) # 增加基因数到 80
groups <- c(
  rep("STC1+ CAFs", 20),
  rep("Normal Fibroblasts", 20),
  rep("Tumor Cells", 20),
  rep("Immune Cells", 20)
) # ✨ 新增第 4 组

test_data_4groups <- data.frame(Gene = genes, Group = groups)
sample_cols <- paste0("Sample_", 1:5)
for (col in sample_cols) {
  test_data_4groups[[col]] <- runif(80, 0, 100)
}

# 增加轨道数据
test_data_4groups$logFC <- rnorm(80, 0, 2)

# 2. 定义 4 个组的颜色字典
my_4group_palette <- c(
  "STC1+ CAFs" = "#DC0000",
  "Normal Fibroblasts" = "#4DBBD5",
  "Tumor Cells" = "#00A087",
  "Immune Cells" = "#E64B35" # ✨ 新增颜色
)

# 3. 轨道配置 (无需修改，逻辑通用)
my_tracks <- list(
  list(
    column = "logFC",
    type = "bars",
    label = "logFC",
    height = 0.1,
    color = function(x) ifelse(x > 0, "#DC0000", "#3C5488")
  )
)

# 4. 运行函数 (自动适配)
plot_advanced_circos(
  data = test_data_4groups,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = sample_cols,
  group_palette = my_4group_palette, # 传入 4 色字典
  track_configs = my_tracks,
  gap_degree = 90 # 依然画 270 度环形
)


# ==============================================================================
# 1. 构造多组、多样本、多维度的模拟数据
# ==============================================================================
set.seed(2026)

# --- A. 基础信息 ---
# 5 个细胞分组 (满足 >4 组的需求)
group_levels <- c("B_Cells", "CD4_T", "CD8_T", "NK_Cells", "Myeloid")
n_genes <- 100 # 总共 100 个基因

# 生成主数据框
df <- data.frame(
  Gene = paste0("Gene_", 1:n_genes),
  Group = factor(
    sample(group_levels, n_genes, replace = TRUE),
    levels = group_levels
  )
)
# 排序以便观察（非必须，因为函数会自动处理，但为了模拟真实数据的有序性）
df <- df[order(df$Group), ]

# --- B. 模拟 6 个样本的表达量热图数据 ---
# 假设前3个样本是对照组，后3个是治疗组
sample_cols <- paste0("Sample_", 1:6)
for (col in sample_cols) {
  # 模拟不同组有不同的表达模式
  df[[col]] <- runif(n_genes, 0, 10) + ifelse(df$Group == "B_Cells", 5, 0)
}

# --- C. 模拟多维度轨道数据 ---
# 1. 差异倍数 (logFC): 有正有负
df$logFC <- rnorm(n_genes, mean = 0, sd = 2)

# 2. 甲基化水平 (Methylation): 0-1 之间
df$Methylation <- runif(n_genes, 0, 1)

# 3. 免疫检查点评分 (Score): 连续变量
df$Immune_Score <- rnorm(n_genes, mean = 50, sd = 15)

# 4. 总体分布数据 (用于箱线图): 模拟每个组的样本内变异
# 注意：箱线图是按组画的，这里我们模拟一个代表组内离散度的指标
df$Variance_Metric <- abs(rnorm(n_genes, mean = 10, sd = 5))


# ==============================================================================
# 2. 定义颜色与配置
# ==============================================================================

# --- A. 分组颜色 (5组配色) ---
my_palette <- c(
  "B_Cells" = "#E64B35", # 红
  "CD4_T" = "#4DBBD5", # 蓝
  "CD8_T" = "#00A087", # 绿
  "NK_Cells" = "#3C5488", # 深蓝
  "Myeloid" = "#F39B7F" # 橙
)

# --- B. 热图颜色 (蓝-白-红 经典配色) ---
heatmap_col_fun <- colorRamp2(c(0, 0.5, 1), c("#2166AC", "white", "#B2182B"))

# --- C. 轨道配置 (4层豪华轨道) ---
my_tracks <- list(
  # 第1层: logFC 柱状图 (红蓝双色)
  list(
    column = "logFC",
    type = "bars",
    height = 0.08,
    label = "logFC",
    color = function(x) ifelse(x > 0, "#E64B35", "#3C5488")
  ),

  # 第2层: 甲基化水平 折线图 (面积填充)
  list(
    column = "Methylation",
    type = "area",
    height = 0.06,
    label = "Methy",
    color = "#00A087",
    scale = "minmax" # 归一化到 0-1
  ),

  # 第3层: 免疫评分 散点图 (点的大小固定，颜色映射)
  list(
    column = "Immune_Score",
    type = "points",
    height = 0.06,
    label = "Score",
    color = colorRamp2(c(20, 80), c("grey", "purple")),
    scale = "zscore" # Z-score 标准化
  ),

  # 第4层: 组内变异 箱线图 (展示每个细胞群的整体离散度)
  list(
    column = "Variance_Metric",
    type = "boxplot",
    height = 0.1,
    label = "Var\n(Box)",
    color = "#8491B4" # 统一颜色
  )
)

# --- D. 跨组连线数据 (模拟细胞通讯) ---
# 随机选取一些基因作为通讯对
set.seed(123)
interaction_data <- data.frame(
  from = sample(df$Gene[df$Group == "Myeloid"], 5), # 髓系细胞发出的信号
  to = sample(df$Gene[df$Group == "CD8_T"], 5), # 作用于 CD8 T 细胞
  color = "#E64B3580" # 半透明红色连线
)


# ==============================================================================
# 3. 执行绘图 (开启聚类测试)
# ==============================================================================

# 运行函数
plot_advanced_circos(
  data = df,
  group_col = "Group",
  feature_col = "Gene",
  value_cols = sample_cols,
  group_palette = my_palette,
  heatmap_col = heatmap_col_fun,
  track_configs = my_tracks,
  links_data = interaction_data,

  # --- 这里的参数展示了函数的灵活性 ---
  show_rownames = T, # 基因太多，隐藏行名
  show_colnames = TRUE, # 显示样本名
  gap_degree = 45, # 开口留 45 度 (画 315 度的大圆环)
  start_degree = 90, # 从 12 点钟方向开始画

  # --- 传递给 circlize::circos.heatmap 的额外参数 ---
  cluster = F, # ✨ 开启聚类！测试连线和轨道是否对齐
  dend.side = "inside", # 聚类树画在内侧
  colnames.cex = 0.5,
  dend.track.height = 0.15 # 树高
)


library(circlize)
library(RColorBrewer)
library(dplyr)
library(scales)
library(ComplexHeatmap)

# ==============================================================================
# 测试数据生成：多组、多维、多样本的单细胞转录组场景
# ==============================================================================

#' String Repeat Operator
#' @keywords internal
`%R%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}

#' Generate Comprehensive Test Data for Circos Plot
#'
#' Simulates a single-cell RNA-seq dataset with multiple cell types,
#' marker genes, and various metadata
#'
#' @param n_celltypes Number of cell types (groups)
#' @param n_genes_per_type Number of marker genes per cell type
#' @param n_samples Number of samples/conditions
#' @param seed Random seed for reproducibility
#' @return Data frame ready for plot_advanced_circos()
#' @export
generate_circos_test_data <- function(
  n_celltypes = 6,
  n_genes_per_type = 15,
  n_samples = 12,
  seed = 42
) {
  set.seed(seed)

  # Define realistic cell types
  default_celltypes <- c(
    "CD4_T",
    "CD8_T",
    "B_cell",
    "Monocyte",
    "NK_cell",
    "Dendritic",
    "Macrophage",
    "Neutrophil"
  )
  celltypes <- default_celltypes[seq_len(min(
    n_celltypes,
    length(default_celltypes)
  ))]

  # If user wants more than 8, add generic names
  if (n_celltypes > length(default_celltypes)) {
    extra <- paste0(
      "CellType_",
      seq(length(default_celltypes) + 1, n_celltypes)
    )
    celltypes <- c(celltypes, extra)
  }

  # Generate gene names
  all_genes <- character()
  gene_celltype_map <- character()

  for (ct in celltypes) {
    genes <- paste0(ct, "_marker_", seq_len(n_genes_per_type))
    all_genes <- c(all_genes, genes)
    gene_celltype_map <- c(gene_celltype_map, rep(ct, n_genes_per_type))
  }

  n_total_genes <- length(all_genes)

  # Create sample metadata
  sample_metadata <- data.frame(
    Sample = paste0("Sample_", sprintf("%02d", seq_len(n_samples))),
    Condition = rep(c("Control", "Treatment"), length.out = n_samples),
    Timepoint = rep(c("0h", "6h", "12h", "24h"), length.out = n_samples),
    stringsAsFactors = FALSE
  )

  # Generate expression matrix
  # Each cell type has high expression in its own markers
  expr_matrix <- matrix(0, nrow = n_total_genes, ncol = n_samples)

  for (i in seq_len(n_total_genes)) {
    celltype <- gene_celltype_map[i]

    # Base expression (low)
    base_expr <- rnorm(n_samples, mean = 2, sd = 0.5)

    # Add cell-type specific signal
    # Simulate biological variation across samples
    signal <- rnorm(n_samples, mean = 8, sd = 1.5)

    # Add treatment effect (some genes respond to treatment)
    treatment_effect <- ifelse(
      sample_metadata$Condition == "Treatment",
      rnorm(n_samples, mean = 2, sd = 0.5),
      0
    )

    # Combine
    expr_matrix[i, ] <- pmax(base_expr + signal + treatment_effect, 0)

    # Add noise to non-marker genes (simulate off-target expression)
    noise_genes <- gene_celltype_map != celltype
    expr_matrix[noise_genes, ] <- expr_matrix[noise_genes, ] *
      runif(sum(noise_genes), 0.1, 0.3)
  }

  # Create main data frame
  data <- data.frame(
    Gene = all_genes,
    CellType = factor(gene_celltype_map, levels = celltypes),
    stringsAsFactors = FALSE
  )

  # Add expression values
  colnames(expr_matrix) <- sample_metadata$Sample
  data <- cbind(data, as.data.frame(expr_matrix))

  # Add gene-level metadata
  data$LogFC <- rnorm(n_total_genes, mean = 0, sd = 2)
  data$Pvalue <- runif(n_total_genes, min = 0.0001, max = 0.1)
  data$NegLog10P <- -log10(data$Pvalue)

  # Add gene scores (e.g., importance scores)
  data$Importance <- abs(rnorm(n_total_genes, mean = 0.5, sd = 0.2))

  # Add pathway enrichment scores (simulate)
  data$Pathway_Score <- runif(n_total_genes, min = 0, max = 1)

  # Add conservation score
  data$Conservation <- rbeta(n_total_genes, shape1 = 2, shape2 = 5)

  # Add chromatin accessibility (ATAC-seq like)
  data$Chromatin_Access <- rgamma(n_total_genes, shape = 2, rate = 2)

  return(list(
    data = data,
    sample_metadata = sample_metadata,
    celltypes = celltypes
  ))
}

#' Generate Realistic Link Data (Gene-Gene Interactions)
#'
#' @param data Main data frame
#' @param n_links Number of links to generate
#' @param inter_celltype_prob Probability of inter-cell-type links
generate_link_data <- function(data, n_links = 20, inter_celltype_prob = 0.3) {
  genes <- data$Gene
  celltypes <- data$CellType

  links <- data.frame(
    from = character(n_links),
    to = character(n_links),
    color = character(n_links),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_links)) {
    # Decide if link is within or between cell types
    inter_type <- runif(1) < inter_celltype_prob

    if (inter_type) {
      # Between different cell types
      from_gene <- sample(genes, 1)
      from_ct <- celltypes[genes == from_gene]
      to_gene <- sample(genes[celltypes != from_ct], 1)
      link_color <- "#FF000040" # Red for inter-type
    } else {
      # Within same cell type
      ct <- sample(unique(celltypes), 1)
      ct_genes <- genes[celltypes == ct]
      if (length(ct_genes) < 2) {
        next
      }
      sampled <- sample(ct_genes, 2, replace = FALSE)
      from_gene <- sampled[1]
      to_gene <- sampled[2]
      link_color <- "#0000FF40" # Blue for intra-type
    }

    links$from[i] <- from_gene
    links$to[i] <- to_gene
    links$color[i] <- link_color
  }

  # Remove any incomplete rows
  links <- links[complete.cases(links), ]

  return(links)
}

# ==============================================================================
# Example 1: Basic Multi-Group Plot
# ==============================================================================

test_example_1_basic <- function() {
  cat("=" %R% 80, "\n")
  cat("Example 1: Basic 6 Cell Types with 12 Samples\n")
  cat("=" %R% 80, "\n\n")

  # Generate data
  test_data <- generate_circos_test_data(
    n_celltypes = 6,
    n_genes_per_type = 12,
    n_samples = 12
  )

  data <- test_data$data
  celltypes <- test_data$celltypes

  # Define custom color palette
  celltype_colors <- setNames(
    RColorBrewer::brewer.pal(length(celltypes), "Set2"),
    celltypes
  )

  # Sample columns
  sample_cols <- grep("^Sample_", colnames(data), value = TRUE)

  cat(
    "Data dimensions:",
    nrow(data),
    "genes x",
    length(sample_cols),
    "samples\n"
  )
  cat("Cell types:", paste(celltypes, collapse = ", "), "\n\n")

  # Plot
  result <- plot_advanced_circos(
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = sample_cols,
    group_palette = celltype_colors,
    show_rownames = FALSE, # Too many genes
    show_colnames = TRUE,
    gap_degree = 60,
    start_degree = 90,
    show_legend = TRUE
  )

  cat("✓ Plot completed successfully\n")
  cat("Number of groups:", result$n_groups, "\n")
  cat("Clustering enabled:", result$clustered, "\n\n")

  return(invisible(result))
}

# ==============================================================================
# Example 2: With Clustering and External Tracks
# ==============================================================================

test_example_2_tracks <- function() {
  cat("=" %R% 80, "\n")
  cat("Example 2: With Clustering + Multiple External Tracks\n")
  cat("=" %R% 80, "\n\n")

  test_data <- generate_circos_test_data(
    n_celltypes = 5,
    n_genes_per_type = 20,
    n_samples = 10
  )

  data <- test_data$data
  celltypes <- test_data$celltypes
  sample_cols <- grep("^Sample_", colnames(data), value = TRUE)

  # Custom colors using the universal palette system
  celltype_colors <- setNames(
    c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")[seq_along(
      celltypes
    )],
    celltypes
  )

  # Define complex track configurations
  track_configs <- list(
    # Track 1: Log Fold Change (bars)
    list(
      column = "LogFC",
      type = "bars",
      height = 0.12,
      color = colorRamp2(c(-3, 0, 3), c("#2166AC", "#F7F7F7", "#B2182B")),
      scale = NULL, # Already in log scale
      cap_outliers = TRUE,
      label = "Log2 FC",
      bar_width = 0.8
    ),

    # Track 2: P-value significance (lollipop)
    list(
      column = "NegLog10P",
      type = "lollipop",
      height = 0.10,
      color = function(x) ifelse(x > 2, "#D62728", "#7F7F7F"), # Red if p < 0.01
      scale = NULL,
      label = "-log10(P)",
      lwd = 1.2
    ),

    # Track 3: Gene importance (area)
    list(
      column = "Importance",
      type = "area",
      height = 0.08,
      color = "#9467BD80",
      scale = "minmax",
      label = "Importance"
    ),

    # Track 4: Chromatin accessibility (lines)
    list(
      column = "Chromatin_Access",
      type = "lines",
      height = 0.08,
      color = "#FF7F0E",
      scale = "minmax",
      label = "ATAC",
      lwd = 2
    )
  )

  cat("Tracks configured:", length(track_configs), "\n")
  for (i in seq_along(track_configs)) {
    cat(sprintf(
      "  Track %d: %s (%s)\n",
      i,
      track_configs[[i]]$label,
      track_configs[[i]]$type
    ))
  }
  cat("\n")

  # Plot with clustering
  result <- plot_advanced_circos(
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = sample_cols,
    group_palette = celltype_colors,
    track_configs = track_configs,
    show_rownames = TRUE,
    show_colnames = TRUE,
    gap_degree = 45,
    start_degree = 90,
    cluster = TRUE, # Enable clustering!
    show_legend = TRUE,
    heatmap_col = c("#24bb7fff", "#F7F7F7", "#B2182B"),
    legend_params = list(
      title = "Expression",
      x = 0.88,
      y = 0.85
    ),
    plot_params = list(
      heatmap_rownames_cex = 0.5,
      track_label_cex = 0.9,
      group_label_cex = 1.3
    )
  )

  cat("✓ Plot with tracks completed\n")
  cat("Clustering:", result$clustered, "\n\n")

  return(invisible(result))
}

# ==============================================================================
# Example 3: Large Scale with Links (>8 Cell Types)
# ==============================================================================

test_example_3_large_scale <- function() {
  cat("=" %R% 80, "\n")
  cat("Example 3: Large Scale (8 Cell Types) with Gene-Gene Links\n")
  cat("=" %R% 80, "\n\n")

  test_data <- generate_circos_test_data(
    n_celltypes = 8,
    n_genes_per_type = 15,
    n_samples = 16
  )

  data <- test_data$data
  celltypes <- test_data$celltypes
  sample_cols <- grep("^Sample_", colnames(data), value = TRUE)

  # Generate links
  links_data <- generate_link_data(
    data,
    n_links = 25,
    inter_celltype_prob = 0.4
  )

  cat("Generated", nrow(links_data), "gene-gene links\n")
  cat("  - Red links: inter-cell-type interactions\n")
  cat("  - Blue links: intra-cell-type interactions\n\n")

  # Use more colors for 8 groups
  celltype_colors <- setNames(
    c(
      "#E64B35",
      "#4DBBD5",
      "#00A087",
      "#3C5488",
      "#F39B7F",
      "#8491B4",
      "#91D1C2",
      "#DC0000"
    ),
    celltypes
  )

  # Simpler track config for large plot
  track_configs <- list(
    list(
      column = "LogFC",
      type = "bars",
      height = 0.10,
      color = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
      label = "Log2FC"
    ),
    list(
      column = "NegLog10P",
      type = "points",
      height = 0.08,
      color = function(x) ifelse(x > 2, "red", "grey60"),
      label = "Significance",
      pch = 16,
      cex = 0.8
    )
  )

  # Plot
  result <- plot_advanced_circos(
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = sample_cols,
    group_palette = celltype_colors,
    track_configs = track_configs,
    links_data = links_data, # Add links!
    show_rownames = FALSE, # Too crowded with 120 genes
    show_colnames = FALSE,
    gap_degree = 30,
    gap_between = 3,
    start_degree = 90,
    cluster = TRUE,
    show_legend = TRUE,
    plot_params = list(
      group_label_cex = 1.0,
      track_point_cex = 0.7
    )
  )

  cat("✓ Large-scale plot with links completed\n")
  cat("Total genes:", nrow(data), "\n")
  cat("Cell types:", result$n_groups, "\n\n")

  # Verify link positions
  cat("Verifying feature positions after clustering...\n")
  feature_positions <- preview_feature_positions(result)

  return(invisible(result))
}

# ==============================================================================
# Example 4: Extreme Case - Maximum Complexity
# ==============================================================================

test_example_4_extreme <- function() {
  cat("=" %R% 80, "\n")
  cat("Example 4: Extreme Complexity Test\n")
  cat("=" %R% 80, "\n\n")

  test_data <- generate_circos_test_data(
    n_celltypes = 7,
    n_genes_per_type = 25,
    n_samples = 20
  )

  data <- test_data$data
  celltypes <- test_data$celltypes
  sample_cols <- grep("^Sample_", colnames(data), value = TRUE)

  links_data <- generate_link_data(
    data,
    n_links = 30,
    inter_celltype_prob = 0.5
  )

  celltype_colors <- setNames(
    RColorBrewer::brewer.pal(length(celltypes), "Dark2"),
    celltypes
  )

  # Maximum track configuration
  track_configs <- list(
    list(
      column = "LogFC",
      type = "bars",
      height = 0.10,
      color = colorRamp2(c(-3, 0, 3), c("#053061", "#F7F7F7", "#67001F")),
      scale = NULL,
      cap_outliers = TRUE,
      label = "Log2FC"
    ),
    list(
      column = "NegLog10P",
      type = "lollipop",
      height = 0.09,
      color = function(x) {
        ifelse(x > 3, "#B2182B", ifelse(x > 2, "#EF8A62", "#CCCCCC"))
      },
      label = "-log10P",
      lwd = 1
    ),
    list(
      column = "Importance",
      type = "area",
      height = 0.07,
      color = "#4575B480",
      scale = "minmax",
      label = "Importance"
    ),
    list(
      column = "Conservation",
      type = "lines",
      height = 0.07,
      color = "#1A9850",
      scale = "minmax",
      label = "Conservation",
      lwd = 1.5
    ),
    list(
      column = "Chromatin_Access",
      type = "points",
      height = 0.06,
      color = "#D73027",
      scale = "zscore",
      cap_outliers = TRUE,
      label = "ATAC",
      pch = 20,
      cex = 0.6
    )
  )

  cat("Configuration:\n")
  cat("  - Genes:", nrow(data), "\n")
  cat("  - Cell types:", length(celltypes), "\n")
  cat("  - Samples:", length(sample_cols), "\n")
  cat("  - Tracks:", length(track_configs), "\n")
  cat("  - Links:", nrow(links_data), "\n")
  cat("  - Clustering: ENABLED\n\n")

  # Plot
  result <- plot_advanced_circos(
    data = data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = sample_cols,
    group_palette = celltype_colors,
    heatmap_col = colorRamp2(
      c(0, 0.25, 0.5),
      c("#313695", "#FFFFBF", "#A50026")
    ),
    track_configs = track_configs,
    links_data = links_data,
    show_rownames = FALSE,
    show_colnames = FALSE,
    gap_degree = 40,
    gap_between = 2,
    start_degree = 85,
    value_range = c(0, 0.5),
    cluster = TRUE,
    show_legend = TRUE,
    legend_params = list(
      title = "Norm. Expr.",
      x = 0.87,
      y = 0.88
    ),
    plot_params = list(
      heatmap_cell_lwd = 0.3,
      heatmap_bg_lwd = 1.5,
      group_label_cex = 1.1,
      track_label_cex = 0.85,
      group_bg_alpha = 0.6
    )
  )

  cat("✓ Extreme complexity plot completed successfully!\n\n")

  return(invisible(result))
}

# ==============================================================================
# Run All Tests
# ==============================================================================

#' Run All Test Examples
#'
#' @param which Integer vector of which examples to run (1-4)
#' @param save_plots Logical. Save plots to PDF?
#' @param output_dir Directory for saving plots
#' @export
run_all_circos_tests <- function(
  which = 1:4,
  save_plots = FALSE,
  output_dir = "circos_test_plots"
) {
  `%R%` <- function(x, n) paste(rep(x, n), collapse = "")

  if (save_plots) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  results <- list()

  if (1 %in% which) {
    if (save_plots) {
      pdf(file.path(output_dir, "test_1_basic.pdf"), width = 10, height = 10)
    }
    results[[1]] <- test_example_1_basic()
    if (save_plots) dev.off()
  }

  if (2 %in% which) {
    if (save_plots) {
      pdf(file.path(output_dir, "test_2_tracks.pdf"), width = 12, height = 12)
    }
    results[[2]] <- test_example_2_tracks()
    if (save_plots) dev.off()
  }

  if (3 %in% which) {
    if (save_plots) {
      pdf(file.path(output_dir, "test_3_large.pdf"), width = 14, height = 14)
    }
    results[[3]] <- test_example_3_large_scale()
    if (save_plots) dev.off()
  }

  if (4 %in% which) {
    if (save_plots) {
      pdf(file.path(output_dir, "test_4_extreme.pdf"), width = 16, height = 16)
    }
    results[[4]] <- test_example_4_extreme()
    if (save_plots) dev.off()
  }

  cat("\n")
  cat("=" %R% 80, "\n")
  cat("All Tests Completed!\n")
  cat("=" %R% 80, "\n\n")

  if (save_plots) {
    cat("Plots saved to:", output_dir, "\n")
  }

  return(invisible(results))
}

# ==============================================================================
# Quick Start
# ==============================================================================

if (FALSE) {
  # Run individual tests
  result1 <- test_example_1_basic()
  result2 <- test_example_2_tracks()
  result3 <- test_example_3_large_scale()
  result4 <- test_example_4_extreme()

  # Run all tests and save to PDF
  all_results <- run_all_circos_tests(save_plots = TRUE)

  # Run only specific tests
  run_all_circos_tests(which = c(2, 4), save_plots = TRUE)

  # Custom test
  my_data <- generate_circos_test_data(
    n_celltypes = 6,
    n_genes_per_type = 20,
    n_samples = 15,
    seed = 123
  )

  my_result <- plot_advanced_circos(
    data = my_data$data,
    group_col = "CellType",
    feature_col = "Gene",
    value_cols = grep("^Sample_", colnames(my_data$data), value = TRUE),
    cluster = TRUE,
    show_legend = TRUE
  )
}

source('环状热图.R')
# 1. 基础测试（6组，12样本）
result1 <- test_example_1_basic()

# 2. 带轨道的测试（5组，10样本，4个外围轨道，聚类）
result2 <- test_example_2_tracks()

# 3. 大规模测试（8组，16样本，25条连线，聚类）
result3 <- test_example_3_large_scale()

# 4. 极限复杂度（7组，20样本，5个轨道，30条连线，聚类）
result4 <- test_example_4_extreme()

# 5. 运行所有测试并保存为PDF
all_results <- run_all_circos_tests(
  save_plots = TRUE,
  output_dir = "my_test_plots"
)

# 6. 只运行特定测试
run_all_circos_tests(which = c(2, 4), save_plots = TRUE)


library(ggplot2)
library(grid)
library(gtable)

# ==============================================================================
# 1. 定义用户接口函数
# ==============================================================================
fix_panel <- function(width = 4, height = 4, margin = NULL) {
  structure(
    list(width = width, height = height, margin = margin),
    class = "fix_panel_params"
  )
}

# ==============================================================================
# 2. 定义 + 号的加法逻辑
# ==============================================================================
#' @method ggplot_add fix_panel_params
#' @export
ggplot_add.fix_panel_params <- function(object, plot, object_name) {
  plot$fix_panel_data <- object
  class(plot) <- c("fixed_panel_plot", class(plot))
  return(plot)
}

# ==============================================================================
# 3. 简化的 Print 方法（直接处理所有逻辑）
# ==============================================================================
#' @method print fixed_panel_plot
#' @export
print.fixed_panel_plot <- function(x, newpage = TRUE, ...) {
  # 提取参数
  params <- x$fix_panel_data
  w <- unit(params$width, "cm")
  h <- unit(params$height, "cm")

  # 临时移除自定义类
  class(x) <- setdiff(class(x), "fixed_panel_plot")

  # 生成 gtable
  g <- ggplotGrob(x)

  # 修改面板尺寸
  panels <- grep("panel", g$layout$name)
  if (length(panels) > 0) {
    panel_index_w <- unique(g$layout$l[panels])
    panel_index_h <- unique(g$layout$t[panels])

    if (length(panel_index_w) == 1) {
      g$widths[panel_index_w] <- w
    } else {
      g$widths[panel_index_w] <- rep(w, length(panel_index_w))
    }

    if (length(panel_index_h) == 1) {
      g$heights[panel_index_h] <- h
    } else {
      g$heights[panel_index_h] <- rep(h, length(panel_index_h))
    }
  }

  # 绘制
  if (newpage) {
    grid.newpage()
  }
  grid.draw(g)

  # 恢复类
  class(x) <- c("fixed_panel_plot", class(x))

  # 返回原始对象
  invisible(x)
}


library(ggplot2)
library(grid)
library(gtable)

# ==============================================================================
# 1. 定义用户接口函数
# ==============================================================================
fix_panel <- function(width = 4, height = 4, margin = NULL) {
  structure(
    list(width = width, height = height, margin = margin),
    class = "fix_panel_params"
  )
}

# ==============================================================================
# 2. 定义 + 号的加法逻辑
# ==============================================================================
#' @method ggplot_add fix_panel_params
#' @export
ggplot_add.fix_panel_params <- function(object, plot, object_name) {
  plot$fix_panel_data <- object
  class(plot) <- c("fixed_panel_plot", class(plot))
  return(plot)
}

# ==============================================================================
# 3. 关键修复：添加 ggplot_gtable 方法
# ==============================================================================
#' @export
ggplot_gtable.fixed_panel_plot <- function(data) {
  # 提取固定尺寸参数
  params <- data$fix_panel_data
  w <- unit(params$width, "cm")
  h <- unit(params$height, "cm")

  # 临时移除自定义类，避免递归
  class(data) <- setdiff(class(data), "fixed_panel_plot")

  # 构建 gtable（使用正确的流程）
  built <- ggplot_build(data)
  g <- ggplot_gtable(built)

  # 恢复类
  class(data) <- c("fixed_panel_plot", class(data))

  # 修改面板尺寸
  panels <- grep("panel", g$layout$name)
  if (length(panels) > 0) {
    panel_index_w <- unique(g$layout$l[panels])
    panel_index_h <- unique(g$layout$t[panels])

    if (length(panel_index_w) == 1) {
      g$widths[panel_index_w] <- w
    } else {
      g$widths[panel_index_w] <- rep(w, length(panel_index_w))
    }

    if (length(panel_index_h) == 1) {
      g$heights[panel_index_h] <- h
    } else {
      g$heights[panel_index_h] <- rep(h, length(panel_index_h))
    }
  }

  return(g)
}

# ==============================================================================
# 4. Print 方法（用于交互查看）
# ==============================================================================
#' @export
print.fixed_panel_plot <- function(x, newpage = is.null(vp), vp = NULL, ...) {
  # 通过 ggplot_gtable 生成修改后的 gtable
  g <- ggplot_gtable(x)

  # 绘制
  if (newpage) {
    grid.newpage()
  }

  if (is.null(vp)) {
    grid.draw(g)
  } else {
    if (is.character(vp)) {
      seekViewport(vp)
    } else {
      pushViewport(vp)
    }
    grid.draw(g)
    upViewport()
  }

  invisible(x)
}
