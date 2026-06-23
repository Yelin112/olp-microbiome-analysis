library(tidyverse)
library(scales)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

# plot_advanced_circos <- function(
#   data,
#   group_col = "Group",
#   feature_col = "Gene",
#   value_cols,
#   group_palette = NULL,
#   heatmap_col = NULL,
#   track_configs = list(),
#   links_data = NULL,
#   show_rownames = TRUE, # ✨ 新增：是否显示行名(基因名)
#   show_colnames = TRUE # ✨ 新增：是否显示列名(样本名)
# ) {
#   # 1. 处理分组因子与颜色映射
#   if (!is.null(group_palette)) {
#     data[[group_col]] <- factor(
#       data[[group_col]],
#       levels = names(group_palette)
#     )
#     group_colors <- group_palette
#   } else {
#     data[[group_col]] <- as.factor(data[[group_col]])
#     group_levels <- levels(data[[group_col]])
#     group_colors <- brewer.pal(max(3, length(group_levels)), "Set2")[
#       1:length(group_levels)
#     ]
#     names(group_colors) <- group_levels
#   }

#   if (is.null(heatmap_col)) {
#     heatmap_col <- colorRamp2(
#       seq(0, 0.5, length.out = 100),
#       colorRampPalette(rev(brewer.pal(n = 5, name = "PiYG")))(100)
#     )
#   }

#   # 2. 数据预处理
#   data_normalized <- data %>%
#     group_by(!!sym(group_col)) %>%
#     mutate(across(all_of(value_cols), ~ rescale(., to = c(0, 0.5)))) %>%
#     ungroup()

#   data_matrix <- data_normalized %>%
#     select(all_of(value_cols)) %>%
#     as.matrix()
#   rownames(data_matrix) <- data_normalized[[feature_col]]

#   # 3. 初始化画布
#   circos.clear()
#   circos.par(
#     start.degree = 90,
#     gap.after = c(rep(5, length(names(group_colors)) - 1), 30),
#     track.margin = c(0.01, 0.01),
#     cell.padding = c(0, 0, 0, 0)
#   )

#   # 4. 绘制核心热图
#   rn_side <- ifelse(show_rownames, "outside", "none") # 控制行名

#   circos.heatmap(
#     data_matrix,
#     split = data_normalized[[group_col]],
#     cluster = FALSE,
#     bg.border = "black",
#     bg.lwd = 1,
#     cell.border = "white",
#     cell.lwd = 0.5,
#     rownames.side = rn_side,
#     rownames.cex = 0.8,
#     col = heatmap_col,
#     track.height = 0.2
#   )

#   # 5. 可选：绘制列名(样本名)轨道
#   if (show_colnames) {
#     circos.track(
#       track.index = get.current.track.index(),
#       bg.border = NA,
#       panel.fun = function(x, y) {
#         if (CELL_META$sector.numeric.index == length(names(group_colors))) {
#           cn <- colnames(data_matrix)
#           n <- length(cn)
#           cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
#           y_coords <- seq(
#             CELL_META$cell.ylim[1] + cell_height / 2,
#             CELL_META$cell.ylim[2] - cell_height / 2,
#             length.out = n
#           )
#           for (i in 1:n) {
#             circos.lines(
#               c(
#                 CELL_META$cell.xlim[2],
#                 CELL_META$cell.xlim[2] + convert_x(1, "mm")
#               ),
#               c(y_coords[i], y_coords[i]),
#               col = "black",
#               lwd = 1
#             )
#           }
#           circos.text(
#             rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
#             y_coords,
#             cn,
#             cex = 1,
#             adj = c(0, 0.5),
#             facing = "inside"
#           )
#         }
#       }
#     )
#   }

#   # 6. 动态生成外围附加轨道
#   if (length(track_configs) > 0) {
#     lapply(track_configs, function(config) {
#       # ✨ 限制引擎：处理极端值与标准化
#       track_data <- data_normalized[[config$column]]

#       # 1) 标准化转换
#       if (!is.null(config$scale)) {
#         if (config$scale == "zscore") {
#           track_data <- as.numeric(scale(track_data))
#         }
#         if (config$scale == "minmax") {
#           track_data <- rescale(track_data, to = c(0, 1))
#         }
#       }

#       # 2) 异常值盖帽 (去除顶端/底端 1% 的极值干扰)
#       if (isTRUE(config$cap_outliers)) {
#         q_low <- quantile(track_data, 0.01, na.rm = TRUE)
#         q_high <- quantile(track_data, 0.99, na.rm = TRUE)
#         track_data[track_data < q_low] <- q_low
#         track_data[track_data > q_high] <- q_high
#       }

#       data_normalized[[config$column]] <- track_data
#       val_range <- range(track_data, na.rm = TRUE)
#       if (config$type %in% c("bars", "lollipop", "area") && val_range[1] > 0) {
#         val_range[1] <- 0
#       }

#       circos.track(
#         ylim = val_range,
#         track.height = config$height,
#         bg.border = "grey80",
#         panel.fun = function(x, y) {
#           sector_data <- data_normalized[
#             data_normalized[[group_col]] == CELL_META$sector.index,
#           ]
#           col_data <- sector_data[[config$column]]
#           x_pos <- CELL_META$xlim[1] +
#             (CELL_META$xlim[2] - CELL_META$xlim[1]) *
#               (1:nrow(sector_data) - 0.5) /
#               nrow(sector_data)
#           point_cols <- if (is.function(config$color)) {
#             config$color(col_data)
#           } else {
#             config$color
#           }

#           # ✨ 新增图形分支
#           if (config$type == "points") {
#             circos.points(
#               x_pos,
#               col_data,
#               pch = 16,
#               cex = 1.2,
#               col = point_cols
#             )
#           } else if (config$type == "bars") {
#             circos.barplot(
#               col_data,
#               pos = x_pos,
#               col = point_cols,
#               bar_width = 0.6
#             )
#           } else if (config$type == "lines") {
#             circos.lines(x_pos, col_data, col = point_cols, lwd = 2)
#           } else if (config$type == "area") {
#             # 面积填充图
#             circos.lines(
#               x_pos,
#               col_data,
#               col = point_cols,
#               area = TRUE,
#               baseline = val_range[1],
#               border = point_cols
#             )
#           } else if (config$type == "rect") {
#             # 矩形区块图
#             circos.rect(
#               xleft = x_pos - 0.4,
#               ybottom = val_range[1],
#               xright = x_pos + 0.4,
#               ytop = col_data,
#               col = point_cols,
#               border = NA
#             )
#           } else if (config$type == "lollipop") {
#             # 棒棒糖图 (线 + 顶部的点)
#             circos.segments(
#               x0 = x_pos,
#               y0 = val_range[1],
#               x1 = x_pos,
#               y1 = col_data,
#               col = point_cols,
#               lwd = 1.5
#             )
#             circos.points(
#               x_pos,
#               col_data,
#               pch = 16,
#               cex = 1.2,
#               col = point_cols
#             )
#           } else if (config$type == "boxplot") {
#             # 注意坐标：使用 CELL_META$xcenter 定位到扇区正中心
#             # col_data 是当前扇区所有基因的该列数值向量
#             # color 取 point_cols 的第一个值，保持整个扇区颜色一致
#             circos.boxplot(
#               col_data,
#               pos = CELL_META$xcenter,
#               col = point_cols[1],
#               border = "black"
#             )
#           } else if (config$type == "violin") {
#             # 同理，绘制小提琴图展示密度分布
#             circos.violin(
#               col_data,
#               pos = CELL_META$xcenter,
#               col = point_cols[1],
#               border = "black"
#             )
#           }

#           # 轨道标签
#           if (CELL_META$sector.numeric.index == length(names(group_colors))) {
#             mid_y <- mean(val_range)
#             circos.lines(
#               c(
#                 CELL_META$cell.xlim[2],
#                 CELL_META$cell.xlim[2] + convert_x(1, "mm")
#               ),
#               c(mid_y, mid_y),
#               col = "black",
#               lwd = 1
#             )
#             circos.text(
#               CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
#               mid_y,
#               config$label,
#               cex = 1,
#               adj = c(0, 0.5),
#               facing = "inside"
#             )
#           }
#         }
#       )
#     })
#   }

#   # 7. 外围分组标签轨道
#   circos.track(
#     ylim = c(0, 1),
#     track.height = 0.08,
#     bg.border = NA,
#     bg.col = adjustcolor(
#       group_colors[levels(data_normalized[[group_col]])],
#       alpha.f = 0.5
#     ),
#     panel.fun = function(x, y) {
#       circos.text(
#         CELL_META$xcenter,
#         0.5,
#         CELL_META$sector.index,
#         facing = "bending.inside",
#         cex = 1.2,
#         font = 2,
#         col = "white",
#         adj = c(0.5, 0.5)
#       )
#     }
#   )

#   # ✨ 8. 核心新增：解析并在中心绘制网络连线
#   if (!is.null(links_data) && nrow(links_data) > 0) {
#     # 构建坐标映射字典：计算每个基因在哪个扇区，以及它的 X 轴绝对坐标
#     feature_mapping <- data_normalized %>%
#       group_by(!!sym(group_col)) %>%
#       mutate(x_pos = row_number() - 0.5) %>% # 计算几何中心坐标
#       ungroup() %>%
#       select(feature = !!sym(feature_col), sector = !!sym(group_col), x_pos)

#     # 遍历外部传入的连线数据并逐条绘制
#     for (i in 1:nrow(links_data)) {
#       from_meta <- feature_mapping[
#         feature_mapping$feature == links_data$from[i],
#       ]
#       to_meta <- feature_mapping[feature_mapping$feature == links_data$to[i], ]

#       # 如果起点和终点都在图上，则画线
#       if (nrow(from_meta) == 1 && nrow(to_meta) == 1) {
#         # 支持自定义连线颜色，如果没有则默认使用带透明度的灰色
#         link_col <- if ("color" %in% colnames(links_data)) {
#           links_data$color[i]
#         } else {
#           "#00000050"
#         }

#         circos.link(
#           sector.index1 = as.character(from_meta$sector), # 起点扇区
#           point1 = from_meta$x_pos, # 起点X坐标
#           sector.index2 = as.character(to_meta$sector), # 终点扇区
#           point2 = to_meta$x_pos, # 终点X坐标
#           col = link_col, # 颜色
#           lwd = 2, # 线宽
#           directional = 1, # ✨ 1表示画箭头，展示调控方向
#           arr.length = 0.2 # 箭头大小
#         )
#       }
#     }
#   }

#   circos.clear()
# }

# # ==============================================================================
# # 依赖包加载
# # ==============================================================================
# library(tidyverse)
# library(scales)
# library(RColorBrewer)
# library(circlize)
# library(ComplexHeatmap)

# # ==============================================================================
# # 核心绘图函数定义
# # ==============================================================================
# # plot_advanced_circos <- function(
# #   data, # 主数据框 (Long/Wide format)
# #   group_col = "Group", # 分组列名
# #   feature_col = "Gene", # 特征列名 (如 Gene Symbol)
# #   value_cols, # 热图数值列名向量 (如 c("S1", "S2"))
# #   group_palette = NULL, # 分组颜色命名向量 (Name -> Color)
# #   heatmap_col = NULL, # 热图颜色映射函数 (colorRamp2)
# #   track_configs = list(), # 外围轨道配置列表
# #   links_data = NULL, # 中心连线数据框 (from, to, color)
# #   show_rownames = TRUE, # 是否显示基因名
# #   show_colnames = TRUE, # 是否显示样本名
# #   gap_degree = 90, # 缺口大小 (默认90度，即绘制270度图)
# #   start_degree = 270 # 起始角度 (配合缺口调整开口方向)
# # ) {
# #   # --------------------------------------------------------------------------
# #   # 1. 预处理：分组因子与颜色校验
# #   # --------------------------------------------------------------------------
# #   if (!is.null(group_palette)) {
# #     # 按照调色板顺序强制设置因子级别
# #     data[[group_col]] <- factor(
# #       data[[group_col]],
# #       levels = names(group_palette)
# #     )
# #     group_colors <- group_palette
# #   } else {
# #     data[[group_col]] <- as.factor(data[[group_col]])
# #     group_levels <- levels(data[[group_col]])
# #     # 默认生成 Set2 配色
# #     group_colors <- brewer.pal(max(3, length(group_levels)), "Set2")[
# #       1:length(group_levels)
# #     ]
# #     names(group_colors) <- group_levels
# #   }

# #   # 默认热图颜色 (粉-绿)
# #   if (is.null(heatmap_col)) {
# #     heatmap_col <- colorRamp2(
# #       seq(0, 0.5, length.out = 100),
# #       colorRampPalette(rev(brewer.pal(n = 5, name = "PiYG")))(100)
# #     )
# #   }

# #   # --------------------------------------------------------------------------
# #   # 2. 预处理：热图数据标准化 (0-0.5)
# #   # --------------------------------------------------------------------------
# #   # 注意：这里只标准化热图数据，外围轨道数据在后面单独处理
# #   data_normalized <- data %>%
# #     group_by(!!sym(group_col)) %>%
# #     mutate(across(all_of(value_cols), ~ rescale(., to = c(0, 0.5)))) %>%
# #     ungroup()

# #   # 转换为矩阵用于 circos.heatmap
# #   data_matrix <- data_normalized %>%
# #     select(all_of(value_cols)) %>%
# #     as.matrix()
# #   rownames(data_matrix) <- data_normalized[[feature_col]]

# #   # --------------------------------------------------------------------------
# #   # 3. 初始化画布：设置角度与间隙
# #   # --------------------------------------------------------------------------
# #   circos.clear() # 清除旧状态

# #   n_groups <- length(names(group_colors))

# #   # 逻辑：组间保留 2 度微小缝隙，最后保留 gap_degree 度大开口
# #   if (n_groups > 1) {
# #     gaps_vector <- c(rep(2, n_groups - 1), gap_degree)
# #   } else {
# #     gaps_vector <- gap_degree
# #   }

# #   circos.par(
# #     start.degree = start_degree,
# #     gap.after = gaps_vector,
# #     track.margin = c(0.01, 0.01), # 轨道间距
# #     cell.padding = c(0, 0, 0, 0) # 单元格内边距
# #   )

# #   # --------------------------------------------------------------------------
# #   # 4. 绘制核心热图 (Inner Track)
# #   # --------------------------------------------------------------------------
# #   rn_side <- ifelse(show_rownames, "outside", "none")

# #   circos.heatmap(
# #     data_matrix,
# #     split = data_normalized[[group_col]],
# #     cluster = FALSE, # 关闭自动聚类，保持数据框顺序
# #     bg.border = "black",
# #     bg.lwd = 1,
# #     cell.border = "white",
# #     cell.lwd = 0.5,
# #     rownames.side = rn_side,
# #     rownames.cex = 0.7,
# #     col = heatmap_col,
# #     track.height = 0.2
# #   )

# #   # --------------------------------------------------------------------------
# #   # 5. 绘制列名/样本名 (Optional)
# #   # --------------------------------------------------------------------------
# #   if (show_colnames) {
# #     circos.track(
# #       track.index = get.current.track.index(),
# #       bg.border = NA,
# #       panel.fun = function(x, y) {
# #         # 只在最后一个扇区绘制列名
# #         if (CELL_META$sector.numeric.index == n_groups) {
# #           cn <- colnames(data_matrix)
# #           n <- length(cn)
# #           # 计算每个样本名的 y 坐标
# #           cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
# #           y_coords <- seq(
# #             CELL_META$cell.ylim[1] + cell_height / 2,
# #             CELL_META$cell.ylim[2] - cell_height / 2,
# #             length.out = n
# #           )

# #           # 画引导线和文字
# #           for (i in 1:n) {
# #             circos.lines(
# #               c(
# #                 CELL_META$cell.xlim[2],
# #                 CELL_META$cell.xlim[2] + convert_x(1, "mm")
# #               ),
# #               c(y_coords[i], y_coords[i]),
# #               col = "black",
# #               lwd = 1
# #             )
# #           }
# #           circos.text(
# #             rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
# #             y_coords,
# #             cn,
# #             cex = 0.8,
# #             adj = c(0, 0.5),
# #             facing = "inside"
# #           )
# #         }
# #       }
# #     )
# #   }

# #   # --------------------------------------------------------------------------
# #   # 6. 动态绘制外围轨道 (Outer Tracks)
# #   # --------------------------------------------------------------------------
# #   if (length(track_configs) > 0) {
# #     lapply(track_configs, function(config) {
# #       # --- A. 数据提取与限制引擎 ---
# #       track_data <- data_normalized[[config$column]]

# #       # 缩放处理 (Scale)
# #       if (!is.null(config$scale)) {
# #         if (config$scale == "zscore") {
# #           track_data <- as.numeric(scale(track_data))
# #         }
# #         if (config$scale == "minmax") {
# #           track_data <- rescale(track_data, to = c(0, 1))
# #         }
# #       }

# #       # 异常值盖帽 (Cap Outliers)
# #       if (isTRUE(config$cap_outliers)) {
# #         q_low <- quantile(track_data, 0.01, na.rm = TRUE)
# #         q_high <- quantile(track_data, 0.99, na.rm = TRUE)
# #         track_data[track_data < q_low] <- q_low
# #         track_data[track_data > q_high] <- q_high
# #       }

# #       # 更新临时数据
# #       data_normalized[[config$column]] <- track_data

# #       # 计算 Y 轴范围
# #       val_range <- range(track_data, na.rm = TRUE)
# #       # 如果是柱状图或面积图，且数据都是正数，强行让底边从0开始
# #       if (config$type %in% c("bars", "lollipop", "area") && val_range[1] > 0) {
# #         val_range[1] <- 0
# #       }

# #       # --- B. 创建轨道与绘图 ---
# #       circos.track(
# #         ylim = val_range,
# #         track.height = config$height,
# #         bg.border = "grey90",
# #         panel.fun = function(x, y) {
# #           # 提取当前扇区数据
# #           sector_data <- data_normalized[
# #             data_normalized[[group_col]] == CELL_META$sector.index,
# #           ]
# #           col_data <- sector_data[[config$column]]

# #           # 计算 X 轴坐标 (基因水平)
# #           x_pos <- CELL_META$xlim[1] +
# #             (CELL_META$xlim[2] - CELL_META$xlim[1]) *
# #               (1:nrow(sector_data) - 0.5) /
# #               nrow(sector_data)

# #           # 颜色处理
# #           point_cols <- if (is.function(config$color)) {
# #             config$color(col_data)
# #           } else {
# #             config$color
# #           }

# #           # --- C. 图形分发逻辑 ---
# #           if (config$type == "points") {
# #             circos.points(x_pos, col_data, pch = 16, cex = 1, col = point_cols)
# #           } else if (config$type == "bars") {
# #             circos.barplot(
# #               col_data,
# #               pos = x_pos,
# #               col = point_cols,
# #               bar_width = 0.6,
# #               border = NA
# #             )
# #           } else if (config$type == "lines") {
# #             circos.lines(x_pos, col_data, col = point_cols, lwd = 1.5)
# #           } else if (config$type == "area") {
# #             circos.lines(
# #               x_pos,
# #               col_data,
# #               col = point_cols,
# #               area = TRUE,
# #               baseline = val_range[1],
# #               border = NA
# #             )
# #           } else if (config$type == "lollipop") {
# #             circos.segments(
# #               x0 = x_pos,
# #               y0 = val_range[1],
# #               x1 = x_pos,
# #               y1 = col_data,
# #               col = point_cols,
# #               lwd = 1
# #             )
# #             circos.points(x_pos, col_data, pch = 16, cex = 1, col = point_cols)
# #           } else if (config$type == "boxplot") {
# #             # 箱线图：绘制在扇区中心 (pos = CELL_META$xcenter)
# #             # 颜色取第一个值作为填充色
# #             circos.boxplot(
# #               col_data,
# #               pos = CELL_META$xcenter,
# #               col = point_cols[1],
# #               border = "black",
# #               outline = FALSE
# #             )
# #           } else if (config$type == "violin") {
# #             # 小提琴图
# #             circos.violin(
# #               col_data,
# #               pos = CELL_META$xcenter,
# #               col = point_cols[1],
# #               border = NA
# #             )
# #           }

# #           # --- D. 添加轨道标签 ---
# #           if (CELL_META$sector.numeric.index == n_groups) {
# #             mid_y <- mean(val_range)
# #             # 画引导线
# #             circos.lines(
# #               c(
# #                 CELL_META$cell.xlim[2],
# #                 CELL_META$cell.xlim[2] + convert_x(1, "mm")
# #               ),
# #               c(mid_y, mid_y),
# #               col = "black",
# #               lwd = 1
# #             )
# #             # 写标签
# #             circos.text(
# #               CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
# #               mid_y,
# #               config$label,
# #               cex = 0.8,
# #               adj = c(0, 0.5),
# #               facing = "inside"
# #             )
# #           }
# #         }
# #       )
# #     })
# #   }

# #   # --------------------------------------------------------------------------
# #   # 7. 绘制外围分组标签 (Labels)
# #   # --------------------------------------------------------------------------
# #   circos.track(
# #     ylim = c(0, 1),
# #     track.height = 0.08,
# #     bg.border = NA,
# #     bg.col = adjustcolor(
# #       group_colors[levels(data_normalized[[group_col]])],
# #       alpha.f = 0.5
# #     ),
# #     panel.fun = function(x, y) {
# #       circos.text(
# #         CELL_META$xcenter,
# #         0.5,
# #         CELL_META$sector.index,
# #         facing = "bending.inside",
# #         niceFacing = TRUE,
# #         cex = 1.2,
# #         font = 2,
# #         col = "white",
# #         adj = c(0.5, 0.5)
# #       )
# #     }
# #   )

# #   # --------------------------------------------------------------------------
# #   # 8. 绘制中心连线 (Links)
# #   # --------------------------------------------------------------------------
# #   if (!is.null(links_data) && nrow(links_data) > 0) {
# #     # 建立坐标映射表
# #     feature_mapping <- data_normalized %>%
# #       group_by(!!sym(group_col)) %>%
# #       mutate(x_pos = row_number() - 0.5) %>% # 计算基因在扇区内的相对位置
# #       ungroup() %>%
# #       select(feature = !!sym(feature_col), sector = !!sym(group_col), x_pos)

# #     # 遍历绘制每一条线
# #     for (i in 1:nrow(links_data)) {
# #       from_meta <- feature_mapping[
# #         feature_mapping$feature == links_data$from[i],
# #       ]
# #       to_meta <- feature_mapping[feature_mapping$feature == links_data$to[i], ]

# #       # 确保起点终点都存在于数据中
# #       if (nrow(from_meta) == 1 && nrow(to_meta) == 1) {
# #         link_col <- if ("color" %in% colnames(links_data)) {
# #           links_data$color[i]
# #         } else {
# #           "#00000050"
# #         }

# #         circos.link(
# #           sector.index1 = as.character(from_meta$sector),
# #           point1 = from_meta$x_pos,
# #           sector.index2 = as.character(to_meta$sector),
# #           point2 = to_meta$x_pos,
# #           col = link_col,
# #           lwd = 1.5,
# #           directional = 1, # 1=箭头指向终点
# #           arr.length = 0.2
# #         )
# #       }
# #     }
# #   }

# #   # 结束绘图，释放参数
# #   circos.clear()
# # }

# library(circlize)
# library(RColorBrewer)
# library(dplyr)
# library(scales)
# library(ComplexHeatmap)
# library(grid)

# # ==============================================================================
# # 1. Helper Functions (工具模块)
# # ==============================================================================

# #' Validate Input Data
# .validate_circos_input <- function(data, group_col, feature_col, value_cols) {
#   if (!is.data.frame(data)) {
#     stop("'data' must be a data.frame")
#   }
#   if (nrow(data) == 0) {
#     stop("'data' is empty")
#   }

#   missing_cols <- setdiff(c(group_col, feature_col, value_cols), colnames(data))
#   if (length(missing_cols) > 0) {
#     stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
#   }

#   # 检查特征列是否有重复 (环状图要求特征唯一)
#   if (anyDuplicated(data[[feature_col]])) {
#     stop(sprintf(
#       "Duplicate values found in '%s'. Please consolidate rows.",
#       feature_col
#     ))
#   }

#   invisible(TRUE)
# }

# #' Prepare Colors
# .prepare_colors <- function(data, group_col, group_palette) {
#   # 确保分组按照因子顺序排列，防止颜色错乱
#   if (!is.null(group_palette)) {
#     data[[group_col]] <- factor(
#       data[[group_col]],
#       levels = names(group_palette)
#     )
#   } else {
#     data[[group_col]] <- as.factor(data[[group_col]])
#     group_levels <- levels(data[[group_col]])
#     n <- length(group_levels)
#     pal <- if (n <= 8) {
#       brewer.pal(max(3, n), "Set2")[1:n]
#     } else {
#       scales::hue_pal()(n)
#     }
#     group_palette <- setNames(pal, group_levels)
#   }
#   return(list(data = data, palette = group_palette))
# }

# #' Generate Legends
# .create_legends <- function(heatmap_col, track_configs, group_palette) {
#   leg_list <- list()

#   # 1. 热图图例
#   leg_list[[1]] <- Legend(
#     title = "Expression",
#     col_fun = heatmap_col,
#     title_gp = gpar(fontsize = 10, fontface = "bold")
#   )

#   # 2. 轨道图例 (自动推断连续或离散)
#   for (cfg in track_configs) {
#     if (!is.null(cfg$label)) {
#       if (is.function(cfg$color)) {
#         # 连续型图例 (简化处理，假设是colorRamp2)
#         # 实际场景可能需要更复杂的判断
#         try(
#           {
#             # 尝试生成一个简单的连续图例
#             at_vals <- seq(0, 1, length.out = 5) # 假定归一化后的范围
#             leg_list[[length(leg_list) + 1]] <- Legend(
#               title = cfg$label,
#               col_fun = cfg$color,
#               at = c(0, 0.5, 1),
#               labels = c("Low", "Mid", "High")
#             )
#           },
#           silent = TRUE
#         )
#       }
#     }
#   }
#   return(leg_list)
# }

# # ==============================================================================
# # 2. Main Function (核心函数)
# # ==============================================================================

# #' Optimized Advanced Circular Heatmap
# #'
# #' @export
# plot_advanced_circos <- function(
#   data,
#   group_col = "Group",
#   feature_col = "Gene",
#   value_cols,
#   group_palette = NULL,
#   heatmap_col = NULL,
#   track_configs = list(),
#   links_data = NULL,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   gap_degree = 90,
#   start_degree = 135,
#   value_range = c(0, 0.5),
#   plot_params = list(),
#   draw_legends = TRUE, # 新增：控制是否绘制图例
#   ...
# ) {
#   # --- 1. 初始化与校验 ---
#   .validate_circos_input(data, group_col, feature_col, value_cols)

#   # 合并绘图参数
#   pp <- modifyList(
#     list(
#       heatmap_bg_border = "black",
#       heatmap_bg_lwd = 1,
#       heatmap_cell_border = "white",
#       heatmap_cell_lwd = 0.5,
#       track_height = 0.08,
#       label_cex = 0.8
#     ),
#     plot_params
#   )

#   # 准备分组颜色与因子水平
#   prep_res <- .prepare_colors(data, group_col, group_palette)
#   data <- prep_res$data
#   group_colors <- prep_res$palette

#   # --- 2. 数据标准化 ---
#   # 仅对热图数据进行标准化
#   data_norm <- data %>%
#     group_by(!!sym(group_col)) %>%
#     mutate(across(
#       all_of(value_cols),
#       ~ scales::rescale(., to = value_range)
#     )) %>%
#     ungroup()

#   mat <- as.matrix(data_norm[, value_cols])
#   rownames(mat) <- data_norm[[feature_col]]

#   # 准备热图颜色
#   if (is.null(heatmap_col)) {
#     heatmap_col <- colorRamp2(
#       seq(value_range[1], value_range[2], length.out = 100),
#       colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
#     )
#   }

#   # --- 3. 初始化画布 ---
#   circos.clear()
#   n_groups <- length(group_colors)
#   # 确保缺口逻辑严密：如果有多个组，前n-1个是小缝，最后一个是大缝
#   gaps <- if (n_groups > 1) c(rep(2, n_groups - 1), gap_degree) else gap_degree

#   circos.par(
#     start.degree = start_degree,
#     gap.after = gaps,
#     track.margin = c(0.01, 0.01),
#     cell.padding = c(0, 0, 0, 0)
#   )

#   # --- 4. 绘制热图 (Capture Row Order) ---
#   # 使用 tryCatch 确保出错也能清理画布
#   tryCatch(
#     {
#       # 获取用户传入的 heatmap 参数，合并默认值
#       args <- list(...)
#       hm_args <- modifyList(
#         list(
#           mat = mat,
#           split = data_norm[[group_col]],
#           col = heatmap_col,
#           track.height = 0.2,
#           bg.border = pp$heatmap_bg_border,
#           bg.lwd = pp$heatmap_bg_lwd,
#           cell.border = pp$heatmap_cell_border,
#           cell.lwd = pp$heatmap_cell_lwd,
#           rownames.side = if (show_rownames) "outside" else "none",
#           rownames.cex = 0.7,
#           cluster = FALSE # 默认为 FALSE，如果用户传入 TRUE 会被 modifyList 覆盖
#         ),
#         args
#       )

#       do.call(circos.heatmap, hm_args)

#       # ✨ 核心修复：捕获聚类后的基因顺序 ✨
#       # circos.heatmap 执行后，如果 cluster=TRUE，扇区内的基因顺序会改变。
#       # 我们必须遍历所有扇区，获取新的 row_order，用于后续的 Tracks 对齐和 Links 坐标。

#       clustered <- isTRUE(hm_args$cluster)
#       sector_levels <- levels(data_norm[[group_col]])

#       # 建立全局映射表：Feature -> (Sector, New_Index_in_Sector)
#       # 这对正确绘制连线至关重要
#       feature_pos_map <- list()

#       for (sec in sector_levels) {
#         # 获取该扇区在 data_norm 中的原始子集
#         sec_data <- data_norm[data_norm[[group_col]] == sec, ]
#         original_features <- sec_data[[feature_col]]

#         # 获取 circos 当前状态下的 row_order (这是相对于该扇区数据的索引)
#         # get.cell.meta.data("row_order") 返回的是当前扇区绘图时数据的排列索引
#         # 注意：需要在对应的 track index 下获取。circos.heatmap 通常占用 Track 1
#         current_ro <- get.cell.meta.data(
#           "row_order",
#           sector.index = sec,
#           track.index = 1
#         )

#         # 记录映射关系
#         # 真实的物理顺序 features：original_features[current_ro]
#         ordered_feats <- original_features[current_ro]

#         # 存入 map
#         for (i in seq_along(ordered_feats)) {
#           feature_pos_map[[ordered_feats[i]]] <- list(
#             sector = sec,
#             x_idx = i,
#             n_total = length(ordered_feats)
#           )
#         }
#       }

#       # --- 5. 绘制列名 ---
#       if (show_colnames) {
#         circos.track(
#           track.index = get.current.track.index(),
#           bg.border = NA,
#           panel.fun = function(x, y) {
#             if (CELL_META$sector.numeric.index == n_groups) {
#               cn <- colnames(mat)
#               n <- length(cn)
#               # 简单的列名绘制逻辑
#               y_pos <- seq(
#                 CELL_META$cell.ylim[1],
#                 CELL_META$cell.ylim[2],
#                 length.out = n
#               )
#               circos.text(
#                 rep(CELL_META$cell.xlim[2] + convert_x(2, "mm"), n),
#                 y_pos,
#                 cn,
#                 cex = 0.8,
#                 adj = c(0, 0.5),
#                 facing = "inside"
#               )
#             }
#           }
#         )
#       }

#       # --- 6. 绘制外部轨道 ---
#       if (length(track_configs) > 0) {
#         for (cfg in track_configs) {
#           # 数据准备
#           track_raw <- data_norm[[cfg$column]]
#           # 处理 Scale
#           if (!is.null(cfg$scale)) {
#             if (cfg$scale == "zscore") {
#               track_raw <- as.numeric(scale(track_raw))
#             }
#             if (cfg$scale == "minmax") {
#               track_raw <- scales::rescale(track_raw, to = c(0, 1))
#             }
#           }
#           # 处理 Outliers
#           if (isTRUE(cfg$cap_outliers)) {
#             q <- quantile(track_raw, c(0.01, 0.99), na.rm = TRUE)
#             track_raw[track_raw < q[1]] <- q[1]
#             track_raw[track_raw > q[2]] <- q[2]
#           }
#           data_norm[[cfg$column]] <- track_raw # 更新临时数据

#           # 确定 Y 轴范围
#           yrange <- range(track_raw, na.rm = TRUE)
#           if (cfg$type %in% c("bars", "area", "lollipop") && yrange[1] > 0) {
#             yrange[1] <- 0
#           }

#           circos.track(
#             ylim = yrange,
#             track.height = cfg$height %||% pp$track_height,
#             bg.border = "grey90",
#             panel.fun = function(x, y) {
#               sec_id <- CELL_META$sector.index
#               # 获取该扇区的原始数据
#               sec_df <- data_norm[data_norm[[group_col]] == sec_id, ]

#               # ✨ 关键：使用与热图一致的 row_order 对数据进行重排 ✨
#               # 无论是否聚类，get.cell.meta.data("row_order") 都是最保险的
#               ro <- get.cell.meta.data(
#                 "row_order",
#                 sector.index = sec_id,
#                 track.index = 1
#               )

#               # 提取并重排数据
#               y_data <- sec_df[[cfg$column]][ro]

#               # 计算 X 坐标 (居中)
#               n <- length(y_data)
#               x_center <- CELL_META$xlim[1] +
#                 (CELL_META$xlim[2] - CELL_META$xlim[1]) * (1:n - 0.5) / n

#               # 颜色映射
#               cols <- if (is.function(cfg$color)) {
#                 cfg$color(y_data)
#               } else {
#                 cfg$color
#               }

#               # 绘图逻辑 (Switch Case)
#               switch(
#                 cfg$type,
#                 points = circos.points(
#                   x_center,
#                   y_data,
#                   pch = 16,
#                   cex = 0.8,
#                   col = cols
#                 ),
#                 bars = circos.barplot(
#                   y_data,
#                   pos = x_center,
#                   col = cols,
#                   bar_width = 0.6,
#                   border = NA
#                 ),
#                 lines = circos.lines(x_center, y_data, col = cols, lwd = 1.5),
#                 area = circos.lines(
#                   x_center,
#                   y_data,
#                   col = cols,
#                   area = TRUE,
#                   baseline = yrange[1],
#                   border = NA
#                 ),
#                 boxplot = circos.boxplot(
#                   y_data,
#                   pos = CELL_META$xcenter,
#                   col = cols[1],
#                   outline = FALSE
#                 ),
#                 violin = circos.violin(
#                   y_data,
#                   pos = CELL_META$xcenter,
#                   col = cols[1],
#                   border = NA
#                 )
#               )

#               # 添加轨道标签 (仅在最后一个扇区)
#               if (CELL_META$sector.numeric.index == n_groups) {
#                 circos.text(
#                   CELL_META$cell.xlim[2] + convert_x(2, "mm"),
#                   mean(yrange),
#                   cfg$label,
#                   cex = pp$label_cex,
#                   adj = c(0, 0.5),
#                   facing = "inside"
#                 )
#               }
#             }
#           )
#         }
#       }

#       # --- 7. 分组标签 ---
#       circos.track(
#         ylim = c(0, 1),
#         track.height = 0.05,
#         bg.border = NA,
#         panel.fun = function(x, y) {
#           circos.text(
#             CELL_META$xcenter,
#             0.5,
#             CELL_META$sector.index,
#             facing = "bending.inside",
#             niceFacing = TRUE,
#             font = 2,
#             cex = 1.2
#           )
#         }
#       )

#       # --- 8. 绘制连线 (使用修正后的坐标) ---
#       if (!is.null(links_data) && nrow(links_data) > 0) {
#         for (i in 1:nrow(links_data)) {
#           u <- links_data$from[i]
#           v <- links_data$to[i]

#           # 从 map 中查找真实位置
#           pos_u <- feature_pos_map[[u]]
#           pos_v <- feature_pos_map[[v]]

#           if (!is.null(pos_u) && !is.null(pos_v)) {
#             # 计算物理坐标： (index - 0.5)
#             pt1 <- pos_u$x_idx - 0.5
#             pt2 <- pos_v$x_idx - 0.5

#             lcol <- if ("color" %in% names(links_data)) {
#               links_data$color[i]
#             } else {
#               "#00000040"
#             }

#             circos.link(
#               pos_u$sector,
#               pt1,
#               pos_v$sector,
#               pt2,
#               col = lcol,
#               lwd = 1.5,
#               directional = 1
#             )
#           }
#         }
#       }
#     },
#     error = function(e) {
#       message("Error occurred in plot_advanced_circos: ", e$message)
#     },
#     finally = {
#       # 这里的 clear 只能在所有绘图结束后调用，但如果我们要画 Legend，
#       # circos.clear() 会重置视口，导致 Legend 画不上去。
#       # 所以通常不在这里 clear，而是让用户画完手动 clear，或者我们只 clear 参数不 clear 画布
#       # 为了安全起见，我们仅在出错时完全重置，正常结束保留状态以便后续添加 grid 元素
#     }
#   )

#   # --- 9. 绘制图例 (使用 ComplexHeatmap Grid 系统) ---
#   if (draw_legends) {
#     lgd_list <- .create_legends(heatmap_col, track_configs, group_palette)
#     # 将图例打包
#     pd <- packLegend(list = lgd_list)
#     # 绘制到右侧
#     draw(pd, x = unit(1, "npc") - unit(2, "mm"), just = "right")
#   }

#   invisible(list(feature_map = feature_pos_map))
# }

# # 简化的空值合并算子
# `%||%` <- function(a, b) if (is.null(a)) b else a

# library(circlize)
# library(RColorBrewer)
# library(dplyr)
# library(scales)

# # ==============================================================================
# # Helper Functions
# # ==============================================================================

# #' Validate Input Data
# #' @keywords internal
# .validate_circos_input <- function(data, group_col, feature_col, value_cols) {
#   if (!is.data.frame(data)) {
#     stop("'data' must be a data.frame")
#   }
#   if (nrow(data) == 0) {
#     stop("'data' is empty")
#   }

#   missing_cols <- setdiff(c(group_col, feature_col, value_cols), colnames(data))
#   if (length(missing_cols) > 0) {
#     stop(sprintf(
#       "Missing columns in data: %s",
#       paste(missing_cols, collapse = ", ")
#     ))
#   }

#   if (anyDuplicated(data[[feature_col]])) {
#     stop(sprintf("Duplicate values found in '%s' column", feature_col))
#   }

#   non_numeric <- value_cols[!sapply(data[value_cols], is.numeric)]
#   if (length(non_numeric) > 0) {
#     stop(sprintf(
#       "Non-numeric value columns: %s",
#       paste(non_numeric, collapse = ", ")
#     ))
#   }

#   invisible(TRUE)
# }

# #' Prepare Color Palette for Groups
# #' @keywords internal
# .prepare_group_colors <- function(data, group_col, group_palette = NULL) {
#   if (!is.null(group_palette)) {
#     if (!is.character(group_palette) || is.null(names(group_palette))) {
#       stop("'group_palette' must be a named character vector")
#     }

#     # 确保因子水平顺序与颜色向量顺序一致
#     data[[group_col]] <- factor(
#       data[[group_col]],
#       levels = names(group_palette)
#     )
#     missing_levels <- setdiff(unique(data[[group_col]]), names(group_palette))
#     if (length(missing_levels) > 0) {
#       warning(sprintf(
#         "Groups without colors in palette: %s",
#         paste(missing_levels, collapse = ", ")
#       ))
#     }
#     return(list(data = data, colors = group_palette))
#   } else {
#     data[[group_col]] <- as.factor(data[[group_col]])
#     group_levels <- levels(data[[group_col]])
#     n_groups <- length(group_levels)

#     if (n_groups <= 12) {
#       colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set3")[seq_len(
#         n_groups
#       )]
#     } else {
#       colors <- scales::hue_pal()(n_groups)
#     }

#     group_colors <- setNames(colors, group_levels)
#     return(list(data = data, colors = group_colors))
#   }
# }

# #' Prepare Heatmap Color Function
# #' @keywords internal
# .prepare_heatmap_colors <- function(
#   heatmap_col = NULL,
#   value_range = c(0, 0.5)
# ) {
#   if (is.null(heatmap_col)) {
#     colors <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
#     heatmap_col <- circlize::colorRamp2(
#       seq(value_range[1], value_range[2], length.out = length(colors)),
#       colors
#     )
#   }

#   if (!is.function(heatmap_col)) {
#     stop(
#       "'heatmap_col' must be a color mapping function (e.g., from colorRamp2)"
#     )
#   }

#   return(heatmap_col)
# }

# #' Process Track Data with Scaling
# #' @keywords internal
# .process_track_data <- function(data, config) {
#   track_data <- data[[config$column]]

#   if (!is.null(config$scale)) {
#     track_data <- switch(
#       config$scale,
#       zscore = as.numeric(scale(track_data)),
#       minmax = scales::rescale(track_data, to = c(0, 1)),
#       log2 = log2(pmax(track_data, 0) + 1),
#       log10 = log10(pmax(track_data, 0) + 1),
#       track_data
#     )
#   }

#   if (isTRUE(config$cap_outliers)) {
#     q_low <- quantile(track_data, 0.01, na.rm = TRUE)
#     q_high <- quantile(track_data, 0.99, na.rm = TRUE)
#     track_data <- pmin(pmax(track_data, q_low), q_high)
#   }

#   return(track_data)
# }

# #' Extract Global Row Order After Clustering
# #'
# #' This is critical for maintaining consistency across tracks and links
# #' @keywords internal
# .extract_row_order <- function(data, group_col, clustered) {
#   if (!clustered) {
#     # No clustering: maintain original order per group
#     row_order_list <- data %>%
#       dplyr::group_by(!!sym(group_col)) %>%
#       dplyr::mutate(
#         sector_index = dplyr::cur_group_id(),
#         row_index = dplyr::row_number()
#       ) %>%
#       dplyr::ungroup() %>%
#       dplyr::select(!!sym(group_col), sector_index, row_index)

#     return(row_order_list)
#   } else {
#     # Clustering: extract actual order from circlize
#     # This must be done immediately after circos.heatmap()
#     row_order_list <- list()

#     sector_names <- get.all.sector.index()

#     for (sector in sector_names) {
#       # Get the actual row order for this sector
#       tryCatch(
#         {
#           # Access the internal structure to get row order
#           # Note: This is somewhat fragile and depends on circlize internals
#           cell_data <- get.cell.meta.data(
#             "row_order",
#             sector.index = sector,
#             track.index = 1
#           )

#           row_order_list[[sector]] <- cell_data
#         },
#         error = function(e) {
#           warning(sprintf(
#             "Could not extract row order for sector '%s': %s",
#             sector,
#             e$message
#           ))
#           # Fallback to sequential order
#           n_rows <- sum(data[[group_col]] == sector)
#           row_order_list[[sector]] <- seq_len(n_rows)
#         }
#       )
#     }

#     # Convert to data frame with global indices
#     order_df <- data %>%
#       dplyr::group_by(!!sym(group_col)) %>%
#       dplyr::mutate(
#         sector = as.character(!!sym(group_col)),
#         original_index = dplyr::row_number()
#       ) %>%
#       dplyr::ungroup()

#     # Map the reordered indices
#     for (sector in names(row_order_list)) {
#       sector_mask <- order_df$sector == sector
#       order_df$new_index[sector_mask] <- row_order_list[[sector]]
#     }

#     return(order_df)
#   }
# }

# #' Calculate Feature Positions After Clustering
# #'
# #' Critical for accurate link drawing
# #' @keywords internal
# .calculate_feature_positions <- function(
#   data,
#   group_col,
#   feature_col,
#   row_order_df
# ) {
#   feature_mapping <- data %>%
#     dplyr::mutate(
#       original_index = dplyr::row_number(),
#       sector = as.character(!!sym(group_col))
#     ) %>%
#     dplyr::select(
#       feature = !!sym(feature_col),
#       sector,
#       original_index
#     )

#   # If we have clustering info, use it
#   if ("new_index" %in% colnames(row_order_df)) {
#     feature_mapping <- feature_mapping %>%
#       dplyr::left_join(
#         row_order_df %>%
#           dplyr::select(sector, original_index, new_index),
#         by = c("sector", "original_index")
#       ) %>%
#       dplyr::group_by(sector) %>%
#       dplyr::arrange(new_index) %>%
#       dplyr::mutate(x_pos = dplyr::row_number() - 0.5) %>%
#       dplyr::ungroup()
#   } else {
#     # No clustering: use original order
#     feature_mapping <- feature_mapping %>%
#       dplyr::group_by(sector) %>%
#       dplyr::mutate(x_pos = dplyr::row_number() - 0.5) %>%
#       dplyr::ungroup()
#   }

#   return(feature_mapping)
# }

# #' Draw Single Track
# #' @keywords internal
# .draw_track <- function(
#   sector_data,
#   config,
#   group_col,
#   n_groups,
#   row_order_info,
#   plot_params
# ) {
#   # Extract current sector's row order
#   current_sector <- CELL_META$sector.index

#   if (!is.null(row_order_info) && "new_index" %in% colnames(row_order_info)) {
#     # Get the reordered indices for this sector
#     sector_order_info <- row_order_info %>%
#       dplyr::filter(sector == current_sector) %>%
#       dplyr::arrange(new_index)

#     # Reorder the data
#     col_data <- sector_data[[config$column]][sector_order_info$original_index]
#   } else {
#     col_data <- sector_data[[config$column]]
#   }

#   n_points <- length(col_data)
#   x_pos <- CELL_META$xlim[1] +
#     (CELL_META$xlim[2] - CELL_META$xlim[1]) *
#       (seq_len(n_points) - 0.5) /
#       n_points

#   point_cols <- if (is.function(config$color)) {
#     config$color(col_data)
#   } else {
#     rep(config$color, length.out = n_points)
#   }

#   val_range <- range(col_data, na.rm = TRUE)

#   switch(
#     config$type,
#     points = {
#       circos.points(
#         x_pos,
#         col_data,
#         pch = config$pch %||% 16,
#         cex = config$cex %||% plot_params$track_point_cex,
#         col = point_cols
#       )
#     },
#     bars = {
#       circos.barplot(
#         col_data,
#         pos = x_pos,
#         col = point_cols,
#         bar_width = config$bar_width %||% 0.6,
#         border = config$border %||% NA
#       )
#     },
#     lines = {
#       circos.lines(
#         x_pos,
#         col_data,
#         col = point_cols,
#         lwd = config$lwd %||% plot_params$track_line_lwd
#       )
#     },
#     area = {
#       baseline <- config$baseline %||% min(val_range[1], 0)
#       circos.lines(
#         x_pos,
#         col_data,
#         col = point_cols,
#         area = TRUE,
#         baseline = baseline,
#         border = NA
#       )
#     },
#     lollipop = {
#       baseline <- config$baseline %||% min(val_range[1], 0)
#       circos.segments(
#         x0 = x_pos,
#         y0 = baseline,
#         x1 = x_pos,
#         y1 = col_data,
#         col = point_cols,
#         lwd = config$lwd %||% 1
#       )
#       circos.points(
#         x_pos,
#         col_data,
#         pch = config$pch %||% 16,
#         cex = config$cex %||% plot_params$track_point_cex,
#         col = point_cols
#       )
#     },
#     boxplot = {
#       circos.boxplot(
#         col_data,
#         pos = CELL_META$xcenter,
#         col = point_cols[1],
#         border = "black",
#         outline = FALSE
#       )
#     },
#     violin = {
#       circos.violin(
#         col_data,
#         pos = CELL_META$xcenter,
#         col = point_cols[1],
#         border = NA
#       )
#     },
#     warning(sprintf("Unknown track type: %s", config$type))
#   )

#   if (CELL_META$sector.numeric.index == n_groups) {
#     mid_y <- mean(val_range)
#     circos.lines(
#       c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
#       c(mid_y, mid_y),
#       col = "black",
#       lwd = 1
#     )
#     circos.text(
#       CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
#       mid_y,
#       config$label,
#       cex = plot_params$track_label_cex,
#       adj = c(0, 0.5),
#       facing = "inside"
#     )
#   }
# }

# #' Add Color Legend
# #' @keywords internal
# .add_heatmap_legend <- function(
#   heatmap_col,
#   value_range,
#   title = "Value",
#   x = 0.85,
#   y = 0.85,
#   width = 0.08,
#   height = 0.5
# ) {
#   # Create legend
#   lgd <- ComplexHeatmap::Legend(
#     col_fun = heatmap_col,
#     title = title,
#     at = seq(value_range[1], value_range[2], length.out = 5),
#     labels_gp = grid::gpar(fontsize = 8),
#     title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
#     grid_height = grid::unit(4, "mm"),
#     grid_width = grid::unit(4, "mm")
#   )

#   # Draw legend
#   ComplexHeatmap::draw(
#     lgd,
#     x = grid::unit(x, "npc"),
#     y = grid::unit(y, "npc"),
#     just = c("left", "top")
#   )
# }

# `%||%` <- function(a, b) if (is.null(a)) b else a

# # ==============================================================================
# # Main Function
# # ==============================================================================

# #' Advanced Circular Heatmap with Tracks
# #'
# #' @param data Data frame with features (rows) and groups
# #' @param group_col Column name for grouping variable
# #' @param feature_col Column name for feature identifiers (must be unique)
# #' @param value_cols Character vector of numeric columns for heatmap
# #' @param group_palette Named character vector of colors for groups (important: order matters!)
# #' @param heatmap_col Color mapping function from colorRamp2()
# #' @param track_configs List of track configurations
# #' @param links_data Data frame with 'from', 'to', and optional 'color' columns
# #' @param show_rownames Logical. Show feature names?
# #' @param show_colnames Logical. Show column names?
# #' @param gap_degree Numeric. Gap after last sector (in degrees)
# #' @param gap_between Numeric. Gap between other sectors (in degrees)
# #' @param start_degree Numeric. Starting angle
# #' @param value_range Numeric vector. Range for data normalization
# #' @param plot_params List of plotting parameters
# #' @param show_legend Logical. Show heatmap legend?
# #' @param legend_params List of legend parameters (x, y, width, height, title)
# #' @param ... Additional arguments passed to circos.heatmap()
# #'
# #' @return Invisible list with plotting information
# #' @export
# plot_advanced_circos <- function(
#   data,
#   group_col = "Group",
#   feature_col = "Gene",
#   value_cols,
#   group_palette = NULL,
#   heatmap_col = NULL,
#   track_configs = list(),
#   links_data = NULL,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   gap_degree = 90,
#   gap_between = 2,
#   start_degree = 135,
#   value_range = c(0, 0.5),
#   plot_params = list(),
#   show_legend = TRUE,
#   legend_params = list(),
#   ...
# ) {
#   # Ensure clean state
#   tryCatch(circos.clear(), error = function(e) NULL)

#   # Use tryCatch wrapper for entire function
#   tryCatch(
#     {
#       # ========================================================================
#       # 1. Input Validation
#       # ========================================================================
#       .validate_circos_input(data, group_col, feature_col, value_cols)

#       # Merge plot parameters
#       default_plot_params <- list(
#         heatmap_rownames_cex = 0.7,
#         heatmap_cell_border = "white",
#         heatmap_cell_lwd = 0.5,
#         heatmap_bg_border = "black",
#         heatmap_bg_lwd = 1,
#         colnames_cex = 0.8,
#         track_point_cex = 1.0,
#         track_line_lwd = 1.5,
#         track_label_cex = 0.8,
#         group_label_cex = 1.2,
#         group_bg_alpha = 0.5
#       )
#       plot_params <- utils::modifyList(default_plot_params, plot_params)

#       # Merge legend parameters
#       default_legend_params <- list(
#         x = 0.85,
#         y = 0.85,
#         width = 0.08,
#         height = 0.5,
#         title = "Value"
#       )
#       legend_params <- utils::modifyList(default_legend_params, legend_params)

#       # ========================================================================
#       # 2. Data Preparation
#       # ========================================================================
#       color_result <- .prepare_group_colors(data, group_col, group_palette)
#       data <- color_result$data # Use the factor-ordered data
#       group_colors <- color_result$colors

#       heatmap_col <- .prepare_heatmap_colors(heatmap_col, value_range)

#       # Normalize data
#       data_normalized <- data %>%
#         dplyr::group_by(!!sym(group_col)) %>%
#         dplyr::mutate(dplyr::across(
#           dplyr::all_of(value_cols),
#           ~ scales::rescale(., to = value_range)
#         )) %>%
#         dplyr::ungroup()

#       # Convert to matrix
#       data_matrix <- data_normalized %>%
#         dplyr::select(dplyr::all_of(value_cols)) %>%
#         as.matrix()
#       rownames(data_matrix) <- data_normalized[[feature_col]]

#       # ========================================================================
#       # 3. Initialize Circos Canvas
#       # ========================================================================
#       circos.clear()

#       # 【修复】确保 gap.after 与因子水平顺序一致
#       group_levels <- levels(data_normalized[[group_col]])
#       n_groups <- length(group_levels)

#       gaps_vector <- rep(gap_between, n_groups)
#       gaps_vector[n_groups] <- gap_degree # Last gap

#       circos.par(
#         start.degree = start_degree,
#         gap.after = gaps_vector,
#         track.margin = c(0.01, 0.01),
#         cell.padding = c(0, 0, 0, 0)
#       )

#       # ========================================================================
#       # 4. Draw Main Heatmap
#       # ========================================================================
#       user_args <- list(...)

#       default_heatmap_args <- list(
#         mat = data_matrix,
#         split = data_normalized[[group_col]], # This is already a properly ordered factor
#         col = heatmap_col,
#         track.height = 0.2,
#         cluster = FALSE,
#         bg.border = plot_params$heatmap_bg_border,
#         bg.lwd = plot_params$heatmap_bg_lwd,
#         cell.border = plot_params$heatmap_cell_border,
#         cell.lwd = plot_params$heatmap_cell_lwd,
#         rownames.side = ifelse(show_rownames, "outside", "none"),
#         rownames.cex = plot_params$heatmap_rownames_cex
#       )

#       final_heatmap_args <- utils::modifyList(default_heatmap_args, user_args)

#       do.call(circos.heatmap, final_heatmap_args)

#       clustered <- isTRUE(final_heatmap_args$cluster)

#       # ========================================================================
#       # 5. Extract Row Order (Critical!)
#       # ========================================================================
#       # This must happen immediately after circos.heatmap()
#       row_order_info <- .extract_row_order(
#         data_normalized,
#         group_col,
#         clustered
#       )

#       # ========================================================================
#       # 6. Add Column Names Track
#       # ========================================================================
#       if (show_colnames) {
#         circos.track(
#           track.index = get.current.track.index(),
#           bg.border = NA,
#           panel.fun = function(x, y) {
#             if (CELL_META$sector.numeric.index == n_groups) {
#               cn <- colnames(data_matrix)
#               n <- length(cn)
#               cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) /
#                 n
#               y_coords <- seq(
#                 CELL_META$cell.ylim[1] + cell_height / 2,
#                 CELL_META$cell.ylim[2] - cell_height / 2,
#                 length.out = n
#               )

#               for (i in seq_len(n)) {
#                 circos.lines(
#                   c(
#                     CELL_META$cell.xlim[2],
#                     CELL_META$cell.xlim[2] + convert_x(1, "mm")
#                   ),
#                   c(y_coords[i], y_coords[i]),
#                   col = "black",
#                   lwd = 1
#                 )
#               }

#               circos.text(
#                 rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
#                 y_coords,
#                 cn,
#                 cex = plot_params$colnames_cex,
#                 adj = c(0, 0.5),
#                 facing = "inside"
#               )
#             }
#           }
#         )
#       }

#       # ========================================================================
#       # 7. Add External Tracks (with correct ordering)
#       # ========================================================================
#       if (length(track_configs) > 0) {
#         for (i in seq_along(track_configs)) {
#           config <- track_configs[[i]]

#           required_fields <- c("column", "type", "height", "label")
#           missing_fields <- setdiff(required_fields, names(config))
#           if (length(missing_fields) > 0) {
#             warning(sprintf(
#               "Track %d missing fields: %s. Skipping.",
#               i,
#               paste(missing_fields, collapse = ", ")
#             ))
#             next
#           }

#           if (!config$column %in% colnames(data_normalized)) {
#             warning(sprintf(
#               "Track column '%s' not found. Skipping.",
#               config$column
#             ))
#             next
#           }

#           track_data <- .process_track_data(data_normalized, config)
#           data_normalized[[config$column]] <- track_data

#           val_range <- range(track_data, na.rm = TRUE)
#           if (
#             config$type %in% c("bars", "lollipop", "area") && val_range[1] > 0
#           ) {
#             val_range[1] <- 0
#           }

#           circos.track(
#             ylim = val_range,
#             track.height = config$height,
#             bg.border = "grey90",
#             panel.fun = function(x, y) {
#               sector_data <- data_normalized[
#                 data_normalized[[group_col]] == CELL_META$sector.index,
#               ]

#               # 【修复】使用全局 row_order_info
#               .draw_track(
#                 sector_data,
#                 config,
#                 group_col,
#                 n_groups,
#                 row_order_info,
#                 plot_params
#               )
#             }
#           )
#         }
#       }

#       # ========================================================================
#       # 8. Add Group Labels Track
#       # ========================================================================
#       circos.track(
#         ylim = c(0, 1),
#         track.height = 0.08,
#         bg.border = NA,
#         bg.col = adjustcolor(
#           group_colors[group_levels],
#           alpha.f = plot_params$group_bg_alpha
#         ),
#         panel.fun = function(x, y) {
#           circos.text(
#             CELL_META$xcenter,
#             0.5,
#             CELL_META$sector.index,
#             facing = "bending.inside",
#             niceFacing = TRUE,
#             cex = plot_params$group_label_cex,
#             font = 2,
#             col = "white",
#             adj = c(0.5, 0.5)
#           )
#         }
#       )

#       # ========================================================================
#       # 9. Add Links (with correct positions!)
#       # ========================================================================
#       if (!is.null(links_data) && nrow(links_data) > 0) {
#         required_link_cols <- c("from", "to")
#         if (!all(required_link_cols %in% colnames(links_data))) {
#           warning(
#             "links_data must contain 'from' and 'to' columns. Skipping links."
#           )
#         } else {
#           # 【修复】使用考虑聚类的位置计算
#           feature_mapping <- .calculate_feature_positions(
#             data_normalized,
#             group_col,
#             feature_col,
#             row_order_info
#           )

#           links_drawn <- 0
#           links_failed <- character()

#           for (i in seq_len(nrow(links_data))) {
#             from_meta <- feature_mapping[
#               feature_mapping$feature == links_data$from[i],
#             ]
#             to_meta <- feature_mapping[
#               feature_mapping$feature == links_data$to[i],
#             ]

#             if (nrow(from_meta) == 1 && nrow(to_meta) == 1) {
#               link_col <- if ("color" %in% colnames(links_data)) {
#                 links_data$color[i]
#               } else {
#                 "#00000050"
#               }

#               tryCatch(
#                 {
#                   circos.link(
#                     sector.index1 = as.character(from_meta$sector),
#                     point1 = from_meta$x_pos,
#                     sector.index2 = as.character(to_meta$sector),
#                     point2 = to_meta$x_pos,
#                     col = link_col,
#                     lwd = 1.5,
#                     directional = 1,
#                     arr.length = 0.2
#                   )
#                   links_drawn <- links_drawn + 1
#                 },
#                 error = function(e) {
#                   links_failed <<- c(
#                     links_failed,
#                     sprintf("%s->%s", links_data$from[i], links_data$to[i])
#                   )
#                 }
#               )
#             } else {
#               links_failed <- c(
#                 links_failed,
#                 sprintf(
#                   "%s->%s (not found)",
#                   links_data$from[i],
#                   links_data$to[i]
#                 )
#               )
#             }
#           }

#           message(sprintf(
#             "Drew %d of %d requested links",
#             links_drawn,
#             nrow(links_data)
#           ))
#           if (length(links_failed) > 0) {
#             warning(sprintf(
#               "Failed to draw links: %s",
#               paste(links_failed, collapse = ", ")
#             ))
#           }
#         }
#       }

#       # ========================================================================
#       # 10. Add Legend
#       # ========================================================================
#       if (show_legend && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#         tryCatch(
#           {
#             .add_heatmap_legend(
#               heatmap_col,
#               value_range,
#               title = legend_params$title,
#               x = legend_params$x,
#               y = legend_params$y,
#               width = legend_params$width,
#               height = legend_params$height
#             )
#           },
#           error = function(e) {
#             warning(sprintf("Failed to add legend: %s", e$message))
#           }
#         )
#       }

#       # ========================================================================
#       # 11. Return Information
#       # ========================================================================
#       result <- list(
#         data_matrix = data_matrix,
#         data_normalized = data_normalized,
#         group_colors = group_colors,
#         group_levels = group_levels,
#         n_groups = n_groups,
#         clustered = clustered,
#         row_order_info = row_order_info,
#         feature_mapping = if (!is.null(links_data)) {
#           .calculate_feature_positions(
#             data_normalized,
#             group_col,
#             feature_col,
#             row_order_info
#           )
#         } else {
#           NULL
#         },
#         plot_params = plot_params
#       )

#       return(invisible(result))
#     },
#     error = function(e) {
#       # Cleanup on error
#       tryCatch(circos.clear(), error = function(e2) NULL)
#       stop(sprintf("Error in plot_advanced_circos: %s", e$message))
#     },
#     finally = {
#       # Always clear at the end
#       circos.clear()
#     }
#   )
# }

# # ==============================================================================
# # Utility Function: Preview Feature Positions
# # ==============================================================================

# #' Preview Feature Positions (Debugging Helper)
# #'
# #' Useful for verifying that clustering hasn't broken your links
# #'
# #' @param result Output from plot_advanced_circos()
# #' @export
# preview_feature_positions <- function(result) {
#   if (is.null(result$feature_mapping)) {
#     stop("No feature mapping available. Did you provide links_data?")
#   }

#   cat("Feature Positions:\n")
#   cat("==================\n\n")

#   if (
#     !is.null(result$row_order_info) &&
#       "new_index" %in% colnames(result$row_order_info)
#   ) {
#     cat("Note: Clustering was enabled. Positions reflect reordered layout.\n\n")
#   }

#   print(
#     result$feature_mapping %>%
#       dplyr::arrange(sector, x_pos) %>%
#       dplyr::select(feature, sector, x_pos)
#   )

#   invisible(result$feature_mapping)
# }

# library(circlize)
# library(RColorBrewer)
# library(dplyr)
# library(scales)

# # ==============================================================================
# # Helper Functions (模块化)
# # ==============================================================================
# #' Validate and Adjust Track Heights
# #'
# #' 确保所有轨道加上热图不会超过可用空间
# #' @keywords internal
# .validate_track_heights <- function(
#   track_configs,
#   heatmap_height = 0.2,
#   show_colnames = TRUE,
#   show_group_labels = TRUE
# ) {
#   # 计算固定占用的空间
#   fixed_space <- heatmap_height # 热图
#   if (show_colnames) {
#     fixed_space <- fixed_space + 0.05
#   } # 列名轨道（估算）
#   if (show_group_labels) {
#     fixed_space <- fixed_space + 0.08
#   } # 分组标签轨道

#   # 计算轨道总高度
#   track_heights <- sapply(track_configs, function(x) x$height %||% 0.1)
#   total_track_height <- sum(track_heights)

#   # 可用空间（保守估计，留20%缓冲）
#   available_space <- 0.8 - fixed_space

#   if (total_track_height > available_space) {
#     warning(sprintf(
#       "Track heights (%.2f) exceed available space (%.2f). Auto-adjusting...",
#       total_track_height,
#       available_space
#     ))

#     # 按比例缩小所有轨道
#     scale_factor <- available_space / total_track_height * 0.95 # 再留5%安全边际

#     for (i in seq_along(track_configs)) {
#       original_height <- track_configs[[i]]$height
#       new_height <- original_height * scale_factor
#       track_configs[[i]]$height <- new_height

#       message(sprintf(
#         "  Track %d (%s): %.3f -> %.3f",
#         i,
#         track_configs[[i]]$label,
#         original_height,
#         new_height
#       ))
#     }
#   }

#   return(track_configs)
# }

# #' Calculate Optimal Track Heights
# #'
# #' 根据轨道数量智能分配高度
# #' @keywords internal
# .suggest_track_heights <- function(n_tracks, total_available = 0.5) {
#   if (n_tracks == 0) {
#     return(numeric(0))
#   }

#   # 根据轨道数量使用不同策略
#   if (n_tracks <= 3) {
#     # 少量轨道：可以给较大空间
#     base_height <- 0.12
#   } else if (n_tracks <= 5) {
#     # 中等数量：均匀分配
#     base_height <- total_available / n_tracks * 0.9
#   } else {
#     # 大量轨道：压缩高度
#     base_height <- total_available / n_tracks * 0.85
#   }

#   # 最小高度限制
#   min_height <- 0.04
#   base_height <- max(base_height, min_height)

#   return(rep(base_height, n_tracks))
# }

# #' Validate Input Data
# #' @keywords internal
# .validate_circos_input <- function(data, group_col, feature_col, value_cols) {
#   # Check data frame
#   if (!is.data.frame(data)) {
#     stop("'data' must be a data.frame")
#   }
#   if (nrow(data) == 0) {
#     stop("'data' is empty")
#   }

#   # Check required columns
#   missing_cols <- setdiff(c(group_col, feature_col, value_cols), colnames(data))
#   if (length(missing_cols) > 0) {
#     stop(sprintf(
#       "Missing columns in data: %s",
#       paste(missing_cols, collapse = ", ")
#     ))
#   }

#   # Check for duplicates
#   if (anyDuplicated(data[[feature_col]])) {
#     stop(sprintf("Duplicate values found in '%s' column", feature_col))
#   }

#   # Check value columns are numeric
#   non_numeric <- value_cols[!sapply(data[value_cols], is.numeric)]
#   if (length(non_numeric) > 0) {
#     stop(sprintf(
#       "Non-numeric value columns: %s",
#       paste(non_numeric, collapse = ", ")
#     ))
#   }

#   invisible(TRUE)
# }

# #' Prepare Color Palette for Groups
# #' @keywords internal
# .prepare_group_colors <- function(data, group_col, group_palette = NULL) {
#   if (!is.null(group_palette)) {
#     # Validate palette
#     if (!is.character(group_palette) || is.null(names(group_palette))) {
#       stop("'group_palette' must be a named character vector")
#     }

#     # 确保因子水平顺序与颜色向量顺序一致
#     data[[group_col]] <- factor(
#       data[[group_col]],
#       levels = names(group_palette)
#     )
#     missing_levels <- setdiff(unique(data[[group_col]]), names(group_palette))
#     if (length(missing_levels) > 0) {
#       warning(sprintf(
#         "Groups without colors in palette: %s",
#         paste(missing_levels, collapse = ", ")
#       ))
#     }
#     return(list(data = data, colors = group_palette))
#   } else {
#     data[[group_col]] <- as.factor(data[[group_col]])
#     group_levels <- levels(data[[group_col]])
#     n_groups <- length(group_levels)

#     if (n_groups <= 12) {
#       colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set3")[seq_len(
#         n_groups
#       )]
#     } else {
#       colors <- scales::hue_pal()(n_groups)
#     }

#     group_colors <- setNames(colors, group_levels)
#     return(list(data = data, colors = group_colors))
#   }
# }

# #' Prepare Heatmap Color Function
# #' @keywords internal
# .prepare_heatmap_colors <- function(
#   heatmap_col = NULL,
#   value_range = c(0, 0.5)
# ) {
#   if (is.null(heatmap_col)) {
#     colors <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
#     heatmap_col <- circlize::colorRamp2(
#       seq(value_range[1], value_range[2], length.out = length(colors)),
#       colors
#     )
#   }

#   if (!is.function(heatmap_col)) {
#     stop(
#       "'heatmap_col' must be a color mapping function (e.g., from colorRamp2)"
#     )
#   }

#   return(heatmap_col)
# }

# #' Process Track Data with Scaling (修复版)
# #' @keywords internal
# .process_track_data <- function(data, config) {
#   track_data <- data[[config$column]]

#   # 首先处理异常值（在缩放之前）
#   if (isTRUE(config$cap_outliers)) {
#     q_low <- quantile(track_data, 0.01, na.rm = TRUE)
#     q_high <- quantile(track_data, 0.99, na.rm = TRUE)
#     track_data <- pmin(pmax(track_data, q_low), q_high)
#   }

#   # 然后应用缩放
#   if (!is.null(config$scale)) {
#     track_data <- switch(
#       config$scale,
#       zscore = as.numeric(scale(track_data)),
#       minmax = scales::rescale(track_data, to = c(0, 1)),
#       log2 = log2(pmax(track_data, 0) + 1),
#       log10 = log10(pmax(track_data, 0) + 1),
#       track_data
#     )
#   }

#   return(track_data)
# }

# #' Extract Global Row Order After Clustering (修复版)
# #' @keywords internal
# .extract_row_order <- function(data, group_col, clustered) {
#   if (!clustered) {
#     # 非聚类模式：返回简单的顺序信息，不包含 new_index
#     row_order_list <- data %>%
#       dplyr::group_by(!!sym(group_col)) %>%
#       dplyr::mutate(
#         sector = as.character(!!sym(group_col)),
#         original_index = dplyr::row_number()
#       ) %>%
#       dplyr::ungroup() %>%
#       dplyr::select(sector, original_index)

#     return(row_order_list)
#   } else {
#     # 聚类模式：需要提取实际的重排序信息
#     order_df <- data %>%
#       dplyr::group_by(!!sym(group_col)) %>%
#       dplyr::mutate(
#         sector = as.character(!!sym(group_col)),
#         original_index = dplyr::row_number()
#       ) %>%
#       dplyr::ungroup()

#     # 尝试从 circlize 获取聚类后的顺序
#     sector_names <- get.all.sector.index()

#     for (sector in sector_names) {
#       tryCatch(
#         {
#           cell_order <- get.cell.meta.data(
#             "row_order",
#             sector.index = sector,
#             track.index = 1
#           )

#           sector_mask <- order_df$sector == sector
#           order_df$new_index[sector_mask] <- cell_order
#         },
#         error = function(e) {
#           warning(sprintf(
#             "Could not extract row order for sector '%s'",
#             sector
#           ))
#           # Fallback: 使用原始顺序
#           sector_mask <- order_df$sector == sector
#           order_df$new_index[sector_mask] <- order_df$original_index[
#             sector_mask
#           ]
#         }
#       )
#     }

#     return(order_df)
#   }
# }

# #' Calculate Feature Positions After Clustering (修复版)
# #' @keywords internal
# .calculate_feature_positions <- function(
#   data,
#   group_col,
#   feature_col,
#   row_order_df
# ) {
#   feature_mapping <- data %>%
#     dplyr::mutate(
#       original_index = dplyr::row_number(),
#       sector = as.character(!!sym(group_col))
#     ) %>%
#     dplyr::select(
#       feature = !!sym(feature_col),
#       sector,
#       original_index
#     )

#   # 如果有聚类信息，使用它
#   if ("new_index" %in% colnames(row_order_df)) {
#     feature_mapping <- feature_mapping %>%
#       dplyr::left_join(
#         row_order_df %>%
#           dplyr::select(sector, original_index, new_index),
#         by = c("sector", "original_index")
#       ) %>%
#       dplyr::group_by(sector) %>%
#       dplyr::arrange(new_index) %>%
#       dplyr::mutate(x_pos = dplyr::row_number() - 0.5) %>%
#       dplyr::ungroup()
#   } else {
#     # 无聚类：使用原始顺序
#     feature_mapping <- feature_mapping %>%
#       dplyr::group_by(sector) %>%
#       dplyr::mutate(x_pos = dplyr::row_number() - 0.5) %>%
#       dplyr::ungroup()
#   }

#   return(feature_mapping)
# }

# #' Draw Single Track (修复版)
# #' @keywords internal
# .draw_track <- function(
#   sector_data,
#   config,
#   group_col,
#   n_groups,
#   row_order_info,
#   plot_params
# ) {
#   current_sector <- CELL_META$sector.index

#   # 检查是否有 new_index 列
#   if (!is.null(row_order_info) && "new_index" %in% colnames(row_order_info)) {
#     # 聚类模式：重新排序数据
#     sector_order_info <- row_order_info %>%
#       dplyr::filter(sector == current_sector) %>%
#       dplyr::arrange(new_index)

#     col_data <- sector_data[[config$column]][sector_order_info$original_index]
#   } else {
#     # 非聚类模式：使用原始顺序
#     col_data <- sector_data[[config$column]]
#   }

#   n_points <- length(col_data)
#   x_pos <- CELL_META$xlim[1] +
#     (CELL_META$xlim[2] - CELL_META$xlim[1]) *
#       (seq_len(n_points) - 0.5) /
#       n_points

#   point_cols <- if (is.function(config$color)) {
#     config$color(col_data)
#   } else {
#     rep(config$color, length.out = n_points)
#   }

#   # 绘图
#   switch(
#     config$type,
#     points = {
#       circos.points(
#         x_pos,
#         col_data,
#         pch = config$pch %||% 16,
#         cex = config$cex %||% plot_params$track_point_cex,
#         col = point_cols
#       )
#     },
#     bars = {
#       circos.barplot(
#         col_data,
#         pos = x_pos,
#         col = point_cols,
#         bar_width = config$bar_width %||% 0.6,
#         border = config$border %||% NA
#       )
#     },
#     lines = {
#       circos.lines(
#         x_pos,
#         col_data,
#         col = point_cols,
#         lwd = config$lwd %||% plot_params$track_line_lwd
#       )
#     },
#     area = {
#       baseline <- config$baseline %||% 0
#       circos.lines(
#         x_pos,
#         col_data,
#         col = point_cols,
#         area = TRUE,
#         baseline = baseline,
#         border = NA
#       )
#     },
#     lollipop = {
#       baseline <- config$baseline %||% 0
#       circos.segments(
#         x0 = x_pos,
#         y0 = baseline,
#         x1 = x_pos,
#         y1 = col_data,
#         col = point_cols,
#         lwd = config$lwd %||% 1
#       )
#       circos.points(
#         x_pos,
#         col_data,
#         pch = config$pch %||% 16,
#         cex = config$cex %||% plot_params$track_point_cex,
#         col = point_cols
#       )
#     },
#     boxplot = {
#       circos.boxplot(
#         col_data,
#         pos = CELL_META$xcenter,
#         col = point_cols[1],
#         border = "black",
#         outline = FALSE
#       )
#     },
#     violin = {
#       circos.violin(
#         col_data,
#         pos = CELL_META$xcenter,
#         col = point_cols[1],
#         border = NA
#       )
#     },
#     warning(sprintf("Unknown track type: %s", config$type))
#   )

#   # 添加轨道标签
#   if (CELL_META$sector.numeric.index == n_groups) {
#     val_range <- CELL_META$ylim
#     mid_y <- mean(val_range)
#     circos.lines(
#       c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
#       c(mid_y, mid_y),
#       col = "black",
#       lwd = 1
#     )
#     circos.text(
#       CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
#       mid_y,
#       config$label,
#       cex = plot_params$track_label_cex,
#       adj = c(0, 0.5),
#       facing = "inside"
#     )
#   }
# }

# #' Add Color Legend
# #' @keywords internal
# .add_heatmap_legend <- function(
#   heatmap_col,
#   value_range,
#   title = "Value",
#   x = 0.85,
#   y = 0.85,
#   width = 0.08,
#   height = 0.5
# ) {
#   if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#     warning("ComplexHeatmap package required for legend")
#     return(invisible(NULL))
#   }

#   # Create legend
#   lgd <- ComplexHeatmap::Legend(
#     col_fun = heatmap_col,
#     title = title,
#     at = seq(value_range[1], value_range[2], length.out = 5),
#     labels_gp = grid::gpar(fontsize = 8),
#     title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
#     grid_height = grid::unit(4, "mm"),
#     grid_width = grid::unit(4, "mm")
#   )

#   # Draw legend
#   ComplexHeatmap::draw(
#     lgd,
#     x = grid::unit(x, "npc"),
#     y = grid::unit(y, "npc"),
#     just = c("left", "top")
#   )
# }

# # Null coalescing operator helper
# `%||%` <- function(a, b) if (is.null(a)) b else a

# # ==============================================================================
# # Main Function
# # ==============================================================================

# #' Advanced Circular Heatmap with Tracks
# #'
# #' @param data Data frame with features (rows) and groups
# #' @param group_col Column name for grouping variable
# #' @param feature_col Column name for feature identifiers (must be unique)
# #' @param value_cols Character vector of numeric columns for heatmap
# #' @param group_palette Named character vector of colors for groups
# #' @param heatmap_col Color mapping function from colorRamp2()
# #' @param track_configs List of track configurations
# #' @param links_data Data frame with 'from', 'to', and optional 'color' columns
# #' @param show_rownames Logical. Show feature names?
# #' @param show_colnames Logical. Show column names?
# #' @param gap_degree Numeric. Gap after last sector (in degrees)
# #' @param gap_between Numeric. Gap between other sectors (in degrees)
# #' @param start_degree Numeric. Starting angle
# #' @param value_range Numeric vector. Range for data normalization
# #' @param plot_params List of plotting parameters
# #' @param show_legend Logical. Show heatmap legend?
# #' @param legend_params List of legend parameters
# #' @param ... Additional arguments passed to circos.heatmap()
# #'
# #' @return Invisible list with plotting information
# #' @export

# # ==============================================================================
# # 修复版：智能处理聚类 + 图例冲突
# # ==============================================================================

# plot_advanced_circos <- function(
#   data,
#   group_col = "Group",
#   feature_col = "Gene",
#   value_cols,
#   group_palette = NULL,
#   heatmap_col = NULL,
#   track_configs = list(),
#   links_data = NULL,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   gap_degree = 90,
#   gap_between = 2,
#   start_degree = 135,
#   value_range = c(0, 0.5),
#   plot_params = list(),
#   show_legend = TRUE,
#   legend_params = list(),
#   auto_adjust_heights = TRUE,
#   ...
# ) {
#   tryCatch(circos.clear(), error = function(e) NULL)

#   tryCatch(
#     {
#       # ========================================================================
#       # 1. Input Validation
#       # ========================================================================
#       .validate_circos_input(data, group_col, feature_col, value_cols)

#       # Merge plot parameters
#       default_plot_params <- list(
#         heatmap_rownames_cex = 0.7,
#         heatmap_cell_border = "white",
#         heatmap_cell_lwd = 0.5,
#         heatmap_bg_border = "black",
#         heatmap_bg_lwd = 1,
#         colnames_cex = 0.8,
#         track_point_cex = 1.0,
#         track_line_lwd = 1.5,
#         track_label_cex = 0.8,
#         group_label_cex = 1.2,
#         group_bg_alpha = 0.5,
#         heatmap_track_height = 0.2
#       )
#       plot_params <- utils::modifyList(default_plot_params, plot_params)

#       default_legend_params <- list(
#         x = 0.85,
#         y = 0.85,
#         width = 0.08,
#         height = 0.5,
#         title = "Value"
#       )
#       legend_params <- utils::modifyList(default_legend_params, legend_params)

#       # ========================================================================
#       # 1.5 检查聚类参数并调整
#       # ========================================================================
#       user_args <- list(...)
#       will_cluster <- isTRUE(user_args$cluster)

#       # 【关键修复】如果要聚类，自动调整参数
#       if (will_cluster) {
#         # 聚类需要额外的 dendrogram 轨道空间
#         dend_height <- user_args$dend.track.height %||% 0.1

#         # 调整热图高度为聚类腾出空间
#         if (plot_params$heatmap_track_height > 0.18) {
#           original_height <- plot_params$heatmap_track_height
#           plot_params$heatmap_track_height <- 0.15
#           message(sprintf(
#             "Clustering enabled: reducing heatmap height %.2f -> %.2f to accommodate dendrogram",
#             original_height,
#             plot_params$heatmap_track_height
#           ))
#         }

#         # 聚类时禁用 ComplexHeatmap 图例（会冲突）
#         if (show_legend) {
#           message(
#             "Note: ComplexHeatmap legend disabled when clustering is enabled (use circlize legend instead)"
#           )
#           show_legend <- FALSE
#         }
#       }

#       # ========================================================================
#       # 1.6 验证和调整轨道高度
#       # ========================================================================
#       if (auto_adjust_heights && length(track_configs) > 0) {
#         # 计算聚类占用的额外空间
#         extra_space <- if (will_cluster) {
#           (user_args$dend.track.height %||% 0.1) + 0.02 # dendrogram + 缓冲
#         } else {
#           0
#         }

#         track_configs <- .validate_track_heights(
#           track_configs,
#           heatmap_height = plot_params$heatmap_track_height + extra_space,
#           show_colnames = show_colnames,
#           show_group_labels = TRUE
#         )
#       }

#       # ========================================================================
#       # 2. Data Preparation
#       # ========================================================================
#       color_result <- .prepare_group_colors(data, group_col, group_palette)
#       data <- color_result$data
#       group_colors <- color_result$colors

#       heatmap_col <- .prepare_heatmap_colors(heatmap_col, value_range)

#       data_normalized <- data %>%
#         dplyr::group_by(!!sym(group_col)) %>%
#         dplyr::mutate(dplyr::across(
#           dplyr::all_of(value_cols),
#           ~ scales::rescale(., to = value_range)
#         )) %>%
#         dplyr::ungroup()

#       data_matrix <- data_normalized %>%
#         dplyr::select(dplyr::all_of(value_cols)) %>%
#         as.matrix()
#       rownames(data_matrix) <- data_normalized[[feature_col]]

#       # ========================================================================
#       # 3. Initialize Circos Canvas
#       # ========================================================================
#       circos.clear()

#       group_levels <- levels(data_normalized[[group_col]])
#       n_groups <- length(group_levels)

#       gaps_vector <- rep(gap_between, n_groups)
#       gaps_vector[n_groups] <- gap_degree

#       # 根据轨道数量和聚类状态调整间距
#       track_margin <- if (length(track_configs) > 5 || will_cluster) {
#         c(0.005, 0.005)
#       } else {
#         c(0.01, 0.01)
#       }

#       circos.par(
#         start.degree = start_degree,
#         gap.after = gaps_vector,
#         track.margin = track_margin,
#         cell.padding = c(0, 0, 0, 0)
#       )

#       # ========================================================================
#       # 4. Draw Main Heatmap
#       # ========================================================================
#       default_heatmap_args <- list(
#         mat = data_matrix,
#         split = data_normalized[[group_col]],
#         col = heatmap_col,
#         track.height = plot_params$heatmap_track_height,
#         cluster = FALSE,
#         bg.border = plot_params$heatmap_bg_border,
#         bg.lwd = plot_params$heatmap_bg_lwd,
#         cell.border = plot_params$heatmap_cell_border,
#         cell.lwd = plot_params$heatmap_cell_lwd,
#         rownames.side = ifelse(show_rownames, "outside", "none"),
#         rownames.cex = plot_params$heatmap_rownames_cex
#       )

#       final_heatmap_args <- utils::modifyList(default_heatmap_args, user_args)

#       do.call(circos.heatmap, final_heatmap_args)

#       clustered <- isTRUE(final_heatmap_args$cluster)

#       # ========================================================================
#       # 5. Extract Row Order
#       # ========================================================================
#       row_order_info <- .extract_row_order(
#         data_normalized,
#         group_col,
#         clustered
#       )

#       # ========================================================================
#       # 6. Add Column Names Track
#       # ========================================================================
#       if (show_colnames) {
#         circos.track(
#           track.index = get.current.track.index(),
#           bg.border = NA,
#           panel.fun = function(x, y) {
#             if (CELL_META$sector.numeric.index == n_groups) {
#               cn <- colnames(data_matrix)
#               n <- length(cn)
#               cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) /
#                 n
#               y_coords <- seq(
#                 CELL_META$cell.ylim[1] + cell_height / 2,
#                 CELL_META$cell.ylim[2] - cell_height / 2,
#                 length.out = n
#               )

#               for (i in seq_len(n)) {
#                 circos.lines(
#                   c(
#                     CELL_META$cell.xlim[2],
#                     CELL_META$cell.xlim[2] + convert_x(1, "mm")
#                   ),
#                   c(y_coords[i], y_coords[i]),
#                   col = "black",
#                   lwd = 1
#                 )
#               }

#               circos.text(
#                 rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
#                 y_coords,
#                 cn,
#                 cex = plot_params$colnames_cex,
#                 adj = c(0, 0.5),
#                 facing = "inside"
#               )
#             }
#           }
#         )
#       }

#       # ========================================================================
#       # 7. Add External Tracks
#       # ========================================================================
#       if (length(track_configs) > 0) {
#         for (i in seq_along(track_configs)) {
#           config <- track_configs[[i]]

#           required_fields <- c("column", "type", "height", "label")
#           missing_fields <- setdiff(required_fields, names(config))
#           if (length(missing_fields) > 0) {
#             warning(sprintf(
#               "Track %d missing fields: %s. Skipping.",
#               i,
#               paste(missing_fields, collapse = ", ")
#             ))
#             next
#           }

#           if (!config$column %in% colnames(data_normalized)) {
#             warning(sprintf(
#               "Track column '%s' not found. Skipping.",
#               config$column
#             ))
#             next
#           }

#           track_data <- .process_track_data(data_normalized, config)
#           data_normalized[[config$column]] <- track_data

#           val_range <- range(track_data, na.rm = TRUE)
#           range_span <- diff(val_range)

#           if (range_span == 0) {
#             val_range <- c(val_range[1] - 0.1, val_range[2] + 0.1)
#           } else {
#             if (config$type %in% c("bars", "lollipop", "area")) {
#               if (val_range[1] > 0) {
#                 val_range[1] <- 0
#               }
#               val_range[2] <- val_range[2] + range_span * 0.15
#             } else {
#               val_range[1] <- val_range[1] - range_span * 0.1
#               val_range[2] <- val_range[2] + range_span * 0.1
#             }
#           }

#           tryCatch(
#             {
#               circos.track(
#                 ylim = val_range,
#                 track.height = config$height,
#                 bg.border = "grey90",
#                 panel.fun = function(x, y) {
#                   sector_data <- data_normalized[
#                     data_normalized[[group_col]] == CELL_META$sector.index,
#                   ]

#                   .draw_track(
#                     sector_data,
#                     config,
#                     group_col,
#                     n_groups,
#                     row_order_info,
#                     plot_params
#                   )
#                 }
#               )
#             },
#             error = function(e) {
#               warning(sprintf(
#                 "Failed to draw track %d (%s): %s",
#                 i,
#                 config$label,
#                 e$message
#               ))
#             }
#           )
#         }
#       }

#       # ========================================================================
#       # 8. Add Group Labels Track
#       # ========================================================================
#       circos.track(
#         ylim = c(0, 1),
#         track.height = 0.08,
#         bg.border = NA,
#         bg.col = adjustcolor(
#           group_colors[group_levels],
#           alpha.f = plot_params$group_bg_alpha
#         ),
#         panel.fun = function(x, y) {
#           circos.text(
#             CELL_META$xcenter,
#             0.5,
#             CELL_META$sector.index,
#             facing = "bending.inside",
#             niceFacing = TRUE,
#             cex = plot_params$group_label_cex,
#             font = 2,
#             col = "white",
#             adj = c(0.5, 0.5)
#           )
#         }
#       )

#       # ========================================================================
#       # 9. Add Links
#       # ========================================================================
#       if (!is.null(links_data) && nrow(links_data) > 0) {
#         required_link_cols <- c("from", "to")
#         if (!all(required_link_cols %in% colnames(links_data))) {
#           warning(
#             "links_data must contain 'from' and 'to' columns. Skipping links."
#           )
#         } else {
#           feature_mapping <- .calculate_feature_positions(
#             data_normalized,
#             group_col,
#             feature_col,
#             row_order_info
#           )

#           links_drawn <- 0
#           links_failed <- character()

#           for (i in seq_len(nrow(links_data))) {
#             from_meta <- feature_mapping[
#               feature_mapping$feature == links_data$from[i],
#             ]
#             to_meta <- feature_mapping[
#               feature_mapping$feature == links_data$to[i],
#             ]

#             if (nrow(from_meta) == 1 && nrow(to_meta) == 1) {
#               link_col <- if ("color" %in% colnames(links_data)) {
#                 links_data$color[i]
#               } else {
#                 "#00000050"
#               }

#               tryCatch(
#                 {
#                   circos.link(
#                     sector.index1 = as.character(from_meta$sector),
#                     point1 = from_meta$x_pos,
#                     sector.index2 = as.character(to_meta$sector),
#                     point2 = to_meta$x_pos,
#                     col = link_col,
#                     lwd = 1.5,
#                     directional = 1,
#                     arr.length = 0.2
#                   )
#                   links_drawn <- links_drawn + 1
#                 },
#                 error = function(e) {
#                   links_failed <<- c(
#                     links_failed,
#                     sprintf("%s->%s", links_data$from[i], links_data$to[i])
#                   )
#                 }
#               )
#             } else {
#               links_failed <- c(
#                 links_failed,
#                 sprintf(
#                   "%s->%s (not found)",
#                   links_data$from[i],
#                   links_data$to[i]
#                 )
#               )
#             }
#           }

#           if (links_drawn > 0) {
#             message(sprintf(
#               "Drew %d of %d requested links",
#               links_drawn,
#               nrow(links_data)
#             ))
#           }
#           if (length(links_failed) > 0 && length(links_failed) <= 5) {
#             warning(sprintf(
#               "Failed to draw links: %s",
#               paste(links_failed, collapse = ", ")
#             ))
#           } else if (length(links_failed) > 5) {
#             warning(sprintf("Failed to draw %d links", length(links_failed)))
#           }
#         }
#       }

#       # ========================================================================
#       # 10. Add Legend (仅在非聚类时)
#       # ========================================================================
#       if (show_legend) {
#         tryCatch(
#           {
#             .add_heatmap_legend(
#               heatmap_col,
#               value_range,
#               title = legend_params$title,
#               x = legend_params$x,
#               y = legend_params$y,
#               width = legend_params$width,
#               height = legend_params$height
#             )
#           },
#           error = function(e) {
#             warning(sprintf("Failed to add legend: %s", e$message))
#           }
#         )
#       }

#       # ========================================================================
#       # 11. Return Information
#       # ========================================================================
#       result <- list(
#         data_matrix = data_matrix,
#         data_normalized = data_normalized,
#         group_colors = group_colors,
#         group_levels = group_levels,
#         n_groups = n_groups,
#         clustered = clustered,
#         row_order_info = row_order_info,
#         feature_mapping = if (!is.null(links_data)) {
#           .calculate_feature_positions(
#             data_normalized,
#             group_col,
#             feature_col,
#             row_order_info
#           )
#         } else {
#           NULL
#         },
#         plot_params = plot_params,
#         track_configs = track_configs
#       )

#       return(invisible(result))
#     },
#     error = function(e) {
#       tryCatch(circos.clear(), error = function(e2) NULL)
#       stop(sprintf("Error in plot_advanced_circos: %s", e$message))
#     },
#     finally = {
#       circos.clear()
#     }
#   )
# }

library(circlize)
library(RColorBrewer)
library(dplyr)
library(scales)
library(ComplexHeatmap)
library(grid)

# ==============================================================================
# 1. Helper Functions (工具模块)
# ==============================================================================

#' Validate Input Data
.validate_circos_input <- function(data, group_col, feature_col, value_cols) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (nrow(data) == 0) {
    stop("'data' is empty")
  }

  missing_cols <- setdiff(c(group_col, feature_col, value_cols), colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # 检查特征列是否有重复 (环状图要求特征唯一)
  if (anyDuplicated(data[[feature_col]])) {
    stop(sprintf(
      "Duplicate values found in '%s'. Please consolidate rows.",
      feature_col
    ))
  }

  invisible(TRUE)
}

#' Prepare Colors
.prepare_colors <- function(data, group_col, group_palette) {
  # 确保分组按照因子顺序排列，防止颜色错乱
  if (!is.null(group_palette)) {
    data[[group_col]] <- factor(
      data[[group_col]],
      levels = names(group_palette)
    )
  } else {
    data[[group_col]] <- as.factor(data[[group_col]])
    group_levels <- levels(data[[group_col]])
    n <- length(group_levels)
    pal <- if (n <= 8) {
      brewer.pal(max(3, n), "Set2")[1:n]
    } else {
      scales::hue_pal()(n)
    }
    group_palette <- setNames(pal, group_levels)
  }
  return(list(data = data, palette = group_palette))
}

#' Generate Legends
.create_legends <- function(heatmap_col, track_configs, group_palette) {
  leg_list <- list()

  # 1. 热图图例
  leg_list[[1]] <- Legend(
    title = "Expression",
    col_fun = heatmap_col,
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )

  # 2. 轨道图例 (自动推断连续或离散)
  for (cfg in track_configs) {
    if (!is.null(cfg$label)) {
      if (is.function(cfg$color)) {
        # 连续型图例 (简化处理，假设是colorRamp2)
        # 实际场景可能需要更复杂的判断
        try(
          {
            # 尝试生成一个简单的连续图例
            at_vals <- seq(0, 1, length.out = 5) # 假定归一化后的范围
            leg_list[[length(leg_list) + 1]] <- Legend(
              title = cfg$label,
              col_fun = cfg$color,
              at = c(0, 0.5, 1),
              labels = c("Low", "Mid", "High")
            )
          },
          silent = TRUE
        )
      }
    }
  }
  return(leg_list)
}

# ==============================================================================
# 2. Main Function (核心函数)
# ==============================================================================

#' Optimized Advanced Circular Heatmap
#'
#' @export
plot_advanced_circos <- function(
  data,
  group_col = "Group",
  feature_col = "Gene",
  value_cols,
  group_palette = NULL,
  heatmap_col = NULL,
  track_configs = list(),
  links_data = NULL,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gap_degree = 90,
  start_degree = 135,
  value_range = c(0, 0.5),
  plot_params = list(),
  draw_legends = TRUE, # 新增：控制是否绘制图例
  ...
) {
  # --- 1. 初始化与校验 ---
  .validate_circos_input(data, group_col, feature_col, value_cols)

  # 合并绘图参数
  pp <- modifyList(
    list(
      heatmap_bg_border = "black",
      heatmap_bg_lwd = 1,
      heatmap_cell_border = "white",
      heatmap_cell_lwd = 0.5,
      track_height = 0.08,
      label_cex = 0.8
    ),
    plot_params
  )

  # 准备分组颜色与因子水平
  prep_res <- .prepare_colors(data, group_col, group_palette)
  data <- prep_res$data
  group_colors <- prep_res$palette

  # --- 2. 数据标准化 ---
  # 仅对热图数据进行标准化
  data_norm <- data %>%
    group_by(!!sym(group_col)) %>%
    mutate(across(
      all_of(value_cols),
      ~ scales::rescale(., to = value_range)
    )) %>%
    ungroup()

  mat <- as.matrix(data_norm[, value_cols])
  rownames(mat) <- data_norm[[feature_col]]

  # 准备热图颜色
  if (is.null(heatmap_col)) {
    heatmap_col <- colorRamp2(
      seq(value_range[1], value_range[2], length.out = 100),
      colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    )
  }

  # --- 3. 初始化画布 ---
  circos.clear()
  n_groups <- length(group_colors)
  # 确保缺口逻辑严密：如果有多个组，前n-1个是小缝，最后一个是大缝
  gaps <- if (n_groups > 1) c(rep(2, n_groups - 1), gap_degree) else gap_degree

  circos.par(
    start.degree = start_degree,
    gap.after = gaps,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  # --- 4. 绘制热图 (Capture Row Order) ---
  # 使用 tryCatch 确保出错也能清理画布
  tryCatch(
    {
      # 获取用户传入的 heatmap 参数，合并默认值
      args <- list(...)
      hm_args <- modifyList(
        list(
          mat = mat,
          split = data_norm[[group_col]],
          col = heatmap_col,
          track.height = 0.2,
          bg.border = pp$heatmap_bg_border,
          bg.lwd = pp$heatmap_bg_lwd,
          cell.border = pp$heatmap_cell_border,
          cell.lwd = pp$heatmap_cell_lwd,
          rownames.side = if (show_rownames) "outside" else "none",
          rownames.cex = 0.7,
          cluster = FALSE # 默认为 FALSE，如果用户传入 TRUE 会被 modifyList 覆盖
        ),
        args
      )

      do.call(circos.heatmap, hm_args)

      # ✨ 核心修复：捕获聚类后的基因顺序 ✨
      # circos.heatmap 执行后，如果 cluster=TRUE，扇区内的基因顺序会改变。
      # 我们必须遍历所有扇区，获取新的 row_order，用于后续的 Tracks 对齐和 Links 坐标。

      clustered <- isTRUE(hm_args$cluster)
      sector_levels <- levels(data_norm[[group_col]])

      # 建立全局映射表：Feature -> (Sector, New_Index_in_Sector)
      # 这对正确绘制连线至关重要
      feature_pos_map <- list()

      for (sec in sector_levels) {
        # 获取该扇区在 data_norm 中的原始子集
        sec_data <- data_norm[data_norm[[group_col]] == sec, ]
        original_features <- sec_data[[feature_col]]

        # 获取 circos 当前状态下的 row_order (这是相对于该扇区数据的索引)
        # get.cell.meta.data("row_order") 返回的是当前扇区绘图时数据的排列索引
        # 注意：需要在对应的 track index 下获取。circos.heatmap 通常占用 Track 1
        current_ro <- get.cell.meta.data(
          "row_order",
          sector.index = sec,
          track.index = 1
        )

        # 记录映射关系
        # 真实的物理顺序 features：original_features[current_ro]
        ordered_feats <- original_features[current_ro]

        # 存入 map
        for (i in seq_along(ordered_feats)) {
          feature_pos_map[[ordered_feats[i]]] <- list(
            sector = sec,
            x_idx = i,
            n_total = length(ordered_feats)
          )
        }
      }

      # --- 5. 绘制列名 ---
      if (show_colnames) {
        circos.track(
          track.index = get.current.track.index(),
          bg.border = NA,
          panel.fun = function(x, y) {
            if (CELL_META$sector.numeric.index == n_groups) {
              cn <- colnames(mat)
              n <- length(cn)
              # 简单的列名绘制逻辑
              y_pos <- seq(
                CELL_META$cell.ylim[1],
                CELL_META$cell.ylim[2],
                length.out = n
              )
              circos.text(
                rep(CELL_META$cell.xlim[2] + convert_x(2, "mm"), n),
                y_pos,
                cn,
                cex = 0.8,
                adj = c(0, 0.5),
                facing = "inside"
              )
            }
          }
        )
      }

      # --- 6. 绘制外部轨道 ---
      if (length(track_configs) > 0) {
        for (cfg in track_configs) {
          # 数据准备
          track_raw <- data_norm[[cfg$column]]
          # 处理 Scale
          if (!is.null(cfg$scale)) {
            if (cfg$scale == "zscore") {
              track_raw <- as.numeric(scale(track_raw))
            }
            if (cfg$scale == "minmax") {
              track_raw <- scales::rescale(track_raw, to = c(0, 1))
            }
          }
          # 处理 Outliers
          if (isTRUE(cfg$cap_outliers)) {
            q <- quantile(track_raw, c(0.01, 0.99), na.rm = TRUE)
            track_raw[track_raw < q[1]] <- q[1]
            track_raw[track_raw > q[2]] <- q[2]
          }
          data_norm[[cfg$column]] <- track_raw # 更新临时数据

          # 确定 Y 轴范围
          yrange <- range(track_raw, na.rm = TRUE)
          if (cfg$type %in% c("bars", "area", "lollipop") && yrange[1] > 0) {
            yrange[1] <- 0
          }

          circos.track(
            ylim = yrange,
            track.height = cfg$height %||% pp$track_height,
            bg.border = "grey90",
            panel.fun = function(x, y) {
              sec_id <- CELL_META$sector.index
              # 获取该扇区的原始数据
              sec_df <- data_norm[data_norm[[group_col]] == sec_id, ]

              # ✨ 关键：使用与热图一致的 row_order 对数据进行重排 ✨
              # 无论是否聚类，get.cell.meta.data("row_order") 都是最保险的
              ro <- get.cell.meta.data(
                "row_order",
                sector.index = sec_id,
                track.index = 1
              )

              # 提取并重排数据
              y_data <- sec_df[[cfg$column]][ro]

              # 计算 X 坐标 (居中)
              n <- length(y_data)
              x_center <- CELL_META$xlim[1] +
                (CELL_META$xlim[2] - CELL_META$xlim[1]) * (1:n - 0.5) / n

              # 颜色映射
              cols <- if (is.function(cfg$color)) {
                cfg$color(y_data)
              } else {
                cfg$color
              }

              # 绘图逻辑 (Switch Case)
              switch(
                cfg$type,
                points = circos.points(
                  x_center,
                  y_data,
                  pch = 16,
                  cex = 0.8,
                  col = cols
                ),
                bars = circos.barplot(
                  y_data,
                  pos = x_center,
                  col = cols,
                  bar_width = 0.6,
                  border = NA
                ),
                lines = circos.lines(x_center, y_data, col = cols, lwd = 1.5),
                area = circos.lines(
                  x_center,
                  y_data,
                  col = cols,
                  area = TRUE,
                  baseline = yrange[1],
                  border = NA
                ),
                boxplot = circos.boxplot(
                  y_data,
                  pos = CELL_META$xcenter,
                  col = cols[1],
                  outline = FALSE
                ),
                violin = circos.violin(
                  y_data,
                  pos = CELL_META$xcenter,
                  col = cols[1],
                  border = NA
                )
              )

              # 添加轨道标签 (仅在最后一个扇区)
              if (CELL_META$sector.numeric.index == n_groups) {
                circos.text(
                  CELL_META$cell.xlim[2] + convert_x(2, "mm"),
                  mean(yrange),
                  cfg$label,
                  cex = pp$label_cex,
                  adj = c(0, 0.5),
                  facing = "inside"
                )
              }
            }
          )
        }
      }

      # --- 7. 分组标签 ---
      circos.track(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = NA,
        panel.fun = function(x, y) {
          circos.text(
            CELL_META$xcenter,
            0.5,
            CELL_META$sector.index,
            facing = "bending.inside",
            niceFacing = TRUE,
            font = 2,
            cex = 1.2
          )
        }
      )

      # --- 8. 绘制连线 (使用修正后的坐标) ---
      if (!is.null(links_data) && nrow(links_data) > 0) {
        for (i in 1:nrow(links_data)) {
          u <- links_data$from[i]
          v <- links_data$to[i]

          # 从 map 中查找真实位置
          pos_u <- feature_pos_map[[u]]
          pos_v <- feature_pos_map[[v]]

          if (!is.null(pos_u) && !is.null(pos_v)) {
            # 计算物理坐标： (index - 0.5)
            pt1 <- pos_u$x_idx - 0.5
            pt2 <- pos_v$x_idx - 0.5

            lcol <- if ("color" %in% names(links_data)) {
              links_data$color[i]
            } else {
              "#00000040"
            }

            circos.link(
              pos_u$sector,
              pt1,
              pos_v$sector,
              pt2,
              col = lcol,
              lwd = 1.5,
              directional = 1
            )
          }
        }
      }
    },
    error = function(e) {
      message("Error occurred in plot_advanced_circos: ", e$message)
    },
    finally = {
      # 这里的 clear 只能在所有绘图结束后调用，但如果我们要画 Legend，
      # circos.clear() 会重置视口，导致 Legend 画不上去。
      # 所以通常不在这里 clear，而是让用户画完手动 clear，或者我们只 clear 参数不 clear 画布
      # 为了安全起见，我们仅在出错时完全重置，正常结束保留状态以便后续添加 grid 元素
    }
  )

  # --- 9. 绘制图例 (使用 ComplexHeatmap Grid 系统) ---
  if (draw_legends) {
    lgd_list <- .create_legends(heatmap_col, track_configs, group_palette)
    # 将图例打包
    pd <- packLegend(list = lgd_list)
    # 绘制到右侧
    draw(pd, x = unit(1, "npc") - unit(2, "mm"), just = "right")
  }

  invisible(list(feature_map = feature_pos_map))
}

# 简化的空值合并算子
`%||%` <- function(a, b) if (is.null(a)) b else a
