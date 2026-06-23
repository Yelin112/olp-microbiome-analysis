# 1. 加载必要的包
library(tidyverse)
library(tidyplots)

# 2. 读取数据 (请将之前提供的整合CSV保存为 'flow_data_combined.csv')
# 或者直接使用以下代码创建数据框：
data <- read.csv("data\\1-29流式.csv")

# 3. 数据预处理：计算 CD8/CD4 比值
data <- data %>%
  mutate(
    CD8_div_CD4 = CD8_pct_CD3 / CD4_pct_CD3,
    # 确保 Group 的顺序是 Control -> F -> P -> R -> S
    Group = factor(Group, levels = c("Control", "F", "P", "R", "S"))
  )
# 4. 定义需要绘制的指标及其对应的 Y 轴标签
# 格式: "列名" = "图表标题/Y轴标签"
metrics_map <- list(
  "CD3_pct" = "CD3+ T Cell Proportion (%)",
  "CD4_pct_CD3" = "CD4+ T Cell (% of CD3)",
  "CD8_pct_CD3" = "CD8+ T Cell (% of CD3)",
  "CD8_div_CD4" = "CD8 / CD4 Ratio",
  "IL17_pct_CD4" = "IL-17+ (% of CD4)",
  "IFNg_pct_CD4" = "IFN-g+ (% of CD4)",
  "IL17_MFI" = "IL-17 MFI (Intensity)",
  "IFNg_MFI" = "IFN-g MFI (Intensity)"
)

# 5. 获取所有组织类型
tissues <- unique(data$Tissue)

# 6. 开始双重循环：遍历组织 -> 遍历指标
for (tis in tissues) {
  # --- 步骤 A: 创建文件夹 ---
  # 去除组织名称中的空格，例如 "Lymph Node" -> "Lymph_Node"
  folder_name <- gsub(" ", "_", tis)

  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
    message(paste("📁 Created folder:", folder_name))
  }

  # --- 步骤 B: 筛选该组织的数据 ---
  df_tissue <- data %>%
    filter(Tissue == tis)

  message(paste("👉 Processing Tissue:", tis))

  # --- 步骤 C: 循环绘制每个指标 ---
  for (metric_col in names(metrics_map)) {
    metric_label <- metrics_map[[metric_col]]

    # 检查数据是否全为空 (避免报错)
    if (all(is.na(df_tissue[[metric_col]]))) {
      message(paste("   ⚠️ Skipping", metric_col, "- No data available"))
      next
    }

    # 过滤掉当前指标为 NA 的行 (统计检验需要)
    df_plot <- df_tissue %>%
      filter(!is.na(!!sym(metric_col)))

    # 开始绘图
    p <- df_plot %>%
      tidyplot(x = Group, y = !!sym(metric_col), color = Group) %>%
      add_mean_bar(alpha = 0.5) %>%
      add_sem_errorbar() %>%
      add_data_points_beeswarm(
        alpha = 0.7,
        white_border = TRUE,
        size = 2 # 稍微调大点点的大小
      ) %>%
      # 尝试添加 P 值 (如果 Control 组有数据且有方差)
      # 使用 tryCatch 防止因数据缺失导致的统计报错中断循环
      {
        tryCatch(
          add_test_pvalue(
            .,
            ref.group = "Control",
            padding_top = 0.1,
            label.size = 3.5
          ),
          error = function(e) {
            message(paste("   Notice: Could not calc p-value for", metric_col))
            . # 返回原图对象，不添加 P 值
          }
        )
      } %>%
      adjust_x_axis_title("Group") %>%
      adjust_y_axis_title(metric_label) %>%
      adjust_title(paste(tis, "-", metric_label)) %>%
      # adjust_colors(colors_discrete_tableau10) %>%
      # 设置单个图的尺寸 (单位 mm)
      adjust_size(width = 80, height = 80)
    # --- 步骤 D: 保存 PDF ---
    file_name <- file.path(folder_name, paste0(metric_col, ".pdf"))
    save_plot(p, file_name)
  }
}

message("✅ All plots generated successfully!")


library(tidyverse)
library(tidyplots)

# ==============================================================================
# 自定义主题与配色设置函数
# ==============================================================================
set_my_theme <- function(
  plot,
  # 1. 默认配色：你可以修改这里的默认 Hex 码
  palette = c(
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
  font_family = "Arial",
  font_face = "bold",
  line_width = 1
) {
  # --- A. 处理颜色逻辑 ---
  # 既然 tidyplots 需要 color_scheme 对象，我们在这里动态生成它
  # 如果传入的是颜色向量 (character)，则转换为 scheme
  if (is.character(palette)) {
    my_scheme <- new_color_scheme(palette, name = "Custom Palette")
  } else {
    # 如果传入的已经是 scheme 对象，直接使用
    my_scheme <- palette
  }

  # --- B. 应用样式 ---
  plot %>%
    # 1. 应用颜色
    adjust_colors(my_scheme) %>%

    # 2. 字体设置 (Arial + Bold)
    adjust_font(family = font_family, face = font_face) %>%

    # # 3. 移除刻度线 (根据你的偏好)
    # remove_x_axis_ticks() %>%
    # remove_y_axis_ticks() %>%

    # 4. 边框设置 (黑色闭合外框)
    adjust_theme_details(
      # panel.border 控制绘图区外框 (黑色，无填充)
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = line_width
      ),
      # axis.line 置空，避免与外框重叠
      axis.line = ggplot2::element_blank(),
      # 确保网格线存在但简洁
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank()
    )
}


# 6. 开始双重循环：遍历组织 -> 遍历指标
for (tis in tissues) {
  # --- 步骤 A: 创建文件夹 ---
  # 去除组织名称中的空格，例如 "Lymph Node" -> "Lymph_Node"
  folder_name <- gsub(" ", "_", tis)

  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
    message(paste("📁 Created folder:", folder_name))
  }

  # --- 步骤 B: 筛选该组织的数据 ---
  df_tissue <- data %>%
    filter(Tissue == tis)

  message(paste("👉 Processing Tissue:", tis))

  # --- 步骤 C: 循环绘制每个指标 ---
  for (metric_col in names(metrics_map)) {
    metric_label <- metrics_map[[metric_col]]

    # 检查数据是否全为空 (避免报错)
    if (all(is.na(df_tissue[[metric_col]]))) {
      message(paste("   ⚠️ Skipping", metric_col, "- No data available"))
      next
    }

    # 过滤掉当前指标为 NA 的行 (统计检验需要)
    df_plot <- df_tissue %>%
      filter(!is.na(!!sym(metric_col)))

    # 开始绘图
    p <- df_plot %>%
      tidyplot(x = Group, y = !!sym(metric_col), color = Group) %>%
      add_mean_bar(alpha = 0.6) %>%
      add_sem_errorbar() %>%
      add_data_points_beeswarm(
        alpha = 0.7,
        white_border = TRUE,
        size = 2 # 稍微调大点点的大小
      ) %>%
      # 尝试添加 P 值 (如果 Control 组有数据且有方差)
      # 使用 tryCatch 防止因数据缺失导致的统计报错中断循环
      {
        tryCatch(
          add_test_pvalue(
            .,
            ref.group = "Control",
            padding_top = 0.1,
            label.size = 3.5
          ),
          error = function(e) {
            message(paste("   Notice: Could not calc p-value for", metric_col))
            . # 返回原图对象，不添加 P 值
          }
        )
      } %>%
      adjust_x_axis_title("Group") %>%
      adjust_y_axis_title(metric_label) %>%
      adjust_title(paste(tis, "-", metric_label)) %>%
      adjust_size(width = 80, height = 80) |>
      set_my_theme()
    # --- 步骤 D: 保存 PDF ---
    file_name <- file.path(folder_name, paste0(metric_col, ".pdf"))
    save_plot(p, file_name)
  }
}


set_my_theme <- function(
  plot,
  palette = c(
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
  font_family = "sans", # 改为 "sans" 以避免 PDF 字体兼容性问题
  font_face = "bold",
  line_width = 1
) {
  # 颜色处理
  if (is.character(palette)) {
    my_scheme <- new_color_scheme(palette, name = "Custom Palette")
  } else {
    my_scheme <- palette
  }

  plot %>%
    adjust_colors(my_scheme) %>%
    adjust_font(family = font_family, face = font_face) %>%

    # 关键：显式定义主题细节，强制显示坐标轴
    adjust_theme_details(
      # 1. 外框：黑色，无填充
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = line_width
      ),

      # 2. 坐标轴线：设为空，避免与外框重叠
      axis.line = ggplot2::element_blank(),

      # 3. 坐标轴刻度文字：强制显示为黑色 (防止被隐藏)
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = 10,
        margin = ggplot2::margin(t = 5)
      ),
      axis.text.y = ggplot2::element_text(
        color = "black",
        size = 10,
        margin = ggplot2::margin(r = 5)
      ),

      # 4. 坐标轴刻度线：强制显示
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),

      # 5. 网格线
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank()
    )
}


# df <- read.csv(
#   "E:\\打工人\\Zeng\\Help\\LM\\Timer2_FTSJD2_STST3.csv",
#   header = TRUE
# )

# # 3. 计算相关性 (Pearson)
# # 计算全局相关性，并在图表标题中展示
# cor_res <- cor.test(df$FTSJD2, df$STAT3)
# title_text <- paste0(
#   "Correlation: R = ",
#   round(cor_res$estimate, 3),
#   ", p = ",
#   signif(cor_res$p.value, 3)
# )

# # 4. 绘图
# p <- df %>%
#   na.omit() |>
#   # 初始化绘图：x轴为FTSJD2，y轴为STAT3，颜色按subtype分组
#   tidyplot(x = FTSJD2, y = STAT3) %>%

#   # 添加散点
#   add_data_points(size = 2, white_border = TRUE) %>%

#   # 添加线性拟合曲线 (Linear Regression) 和置信区间
#   # method = "lm" 表示线性模型
#   # 默认会显示 95% 置信区间 (灰色阴影区域)
#   add_curve_fit(method = "lm", alpha = 0.2) %>%

#   # 调整标题和标签
#   adjust_title("FTSJD2 vs STAT3 Linear Regression") %>%
#   adjust_x_axis_title("FTSJD2 Expression") %>%
#   adjust_y_axis_title("STAT3 Expression") %>%

#   # 添加相关性统计结果到副标题
#   adjust_caption(title_text) %>%

#   # 美化
#   adjust_colors(colors_discrete_friendly) %>%
#   adjust_size(width = 100, height = 100) %>%
#   set_my_theme()

# # 5. 保存或显示
# print(p)
# # save_plot(p, "linear_correlation_plot.pdf")

# # ==============================================================================
# # 3. 绘图函数 (封装通用逻辑)
# # ==============================================================================
# plot_correlation <- function(data, title_suffix, file_name) {
#   # 计算相关性
#   cor_res <- cor.test(data$FTSJD2, data$STAT3)
#   stats_text <- paste0(
#     "Pearson R = ",
#     round(cor_res$estimate, 3),
#     ", p = ",
#     signif(cor_res$p.value, 3)
#   )

#   # 绘图
#   p <- data %>%
#     tidyplot(x = STAT3, y = FTSJD2, color = subtype) %>% # 即使单组也设置color以应用配色

#     # 添加图层
#     add_data_points(size = 1.5, white_border = TRUE, alpha = 0.5) %>%
#     add_curve_fit(method = "lm", alpha = 0.2) %>% # 线性拟合 + 置信区间

#     # 调整标题和标签
#     adjust_title(paste("Correlation Analysis -", title_suffix)) %>%
#     adjust_caption(stats_text) %>% # 在图下显示统计值
#     adjust_x_axis_title("STAT3 Expression") %>%
#     adjust_y_axis_title("CMTR1 Expression") %>%

#     # 调整尺寸
#     adjust_size(width = 100, height = 100) %>%

#     # 应用自定义主题 (Seaside配色)
#     set_my_theme(palette = colors_discrete_friendly)

#   # 保存
#   # 使用 cairo_pdf 以确保字体正确渲染
#   save_plot(p, filename = file_name)

#   print(p) # 在 RStudio 中显示预览
# }

# # ==============================================================================
# # 4. 执行绘图
# # ==============================================================================

# # (1) 整体数据 (不论 Subset)
# message("👉 Plotting Overall Data...")
# df <- df |> na.omit()
# # 为了让 Overall 图也有统一颜色，我们给数据加一个虚拟分组或者直接应用配色
# df_overall <- df |> mutate(subtype = "All Samples")
# plot_correlation(df_overall, "All Samples", "Corr_Overall.svg")

# # (2) 按 Subset 分开
# subtypes <- unique(df$subtype)

# for (st in subtypes) {
#   message(paste("👉 Plotting Subset:", st, "(SVG)..."))

#   df_subset <- df %>% filter(subtype == st)

#   # ✅ 关键修复：显式替换 + 和 - 为 Pos 和 Neg
#   clean_name <- st %>%
#     gsub("\\+", "_Pos", .) %>% # 将 + 替换为 _Pos (注意转义)
#     gsub("-", "_Neg", .) %>% # 将 - 替换为 _Neg
#     gsub("[^a-zA-Z0-9_]", "", .) # 清理其他非字母数字字符

#   # 最终文件名类似于: Corr_HPV_Pos.svg 和 Corr_HPV_Neg.svg
#   file_name <- paste0("Corr_", clean_name, ".svg")

#   plot_correlation(df_subset, st, file_name)
# }

# message("✅ All SVG files generated successfully (Filename conflict resolved)!")

# 1. 确保主题系统代码已加载
source("Scripts\\tidyplot_themes.R")


MY_LAB_THEME <- create_custom_theme(
  name = "Lab Standard Theme",
  palette = c(
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148"
  ), # ✅ 原来的9色配色
  font_family = "sans", # ✅ 与原代码一致
  font_face = "bold", # ✅ 与原代码一致
  font_size = 10, # ✅ 与原代码一致
  border_width = 1, # ✅ 与原代码一致
  border_color = "black",
  grid_major = "grey92", # ✅ 与原代码一致
  grid_minor = FALSE,
  axis_line = FALSE, # ✅ 与原代码一致
  legend_position = "right",
  background = "white"
)

# 2. 运行分析脚本
source("流式\\Control_R.R")

# 3. 查看结果
# 所有图表会自动保存在 Control_vs_R_Analysis/ 文件夹中
