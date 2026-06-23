setwd("E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\help_others")

# 加载必要的包
library(ggplot2)
library(ggthemes)
source('E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\Scripts\\my_themes.R')

# 1. 读取数据
# 请确保工作目录已设置到文件所在的路径，或者提供完整路径
df <- read.csv("redigitized_scatter_data.csv")

# 2. 计算 Pearson 相关系数 (R) 和 P 值
cor_test <- cor.test(df$log2_PSME3_TPM, df$log2_PGP9_5_TPM, method = "pearson")
r_value <- round(cor_test$estimate, 3)
p_value <- format.pval(cor_test$p.value, digits = 3, eps = 0.001) # 格式化P值，极小值显示为 <0.001

# 构建要显示在图上的标注文本
annot_text <- sprintf(
  "R = %s\nP %s",
  r_value,
  ifelse(grepl("<", p_value), p_value, paste0("= ", p_value))
)

# 3. 绘制图形
p <- ggplot(df, aes(x = log2_PSME3_TPM, y = log2_PGP9_5_TPM)) +
  # 绘制散点：shape=21 允许同时设置边框和填充色
  # fill 设置内部颜色，color 设置边框为白色，alpha 设置透明度，stroke 设置边框粗细
  geom_point(
    shape = 21,
    fill = "#3498db",
    color = "white",
    size = 3,
    alpha = 0.7,
    stroke = 0.8
  ) +

  # 添加线性回归拟合线及置信区间 (se = TRUE 默认开启)
  geom_smooth(
    method = "lm",
    color = "#3f5da9ff",
    fill = "gray70",
    linetype = "solid"
  ) +

  # 在图的左上角添加 R 和 P 值标注
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = annot_text,
    hjust = -0.2,
    vjust = 1.5,
    size = 5,
    fontface = "italic"
  ) +

  # 设置坐标轴标签
  labs(
    x = bquote("log"[2] * "(PSME3 TPM)"),
    y = bquote("log"[2] * "(PGP9.5 TPM)")
  ) +

  # 设置整体主题风格
  theme_par() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  fix_panel(
    width = 7,
    height = 5
  )

# 4. 显示图形
print(p)

# 如果需要保存为高清 PDF 或 PNG，可以取消下面代码的注释：
ggsave("Correlation_plot.svg", plot = p)


# 1. 读取数据
# 请确保工作目录已设置到文件所在的路径，或者提供完整路径
df <- read.csv("PA28g_neural_scatter_approx_data.csv")

# 2. 计算 Pearson 相关系数 (R) 和 P 值
cor_test <- cor.test(df$log2_PSME3_TPM, df$log2_TUBB3_TPM, method = "pearson")
r_value <- round(cor_test$estimate, 3)
p_value <- format.pval(cor_test$p.value, digits = 3, eps = 0.001) # 格式化P值，极小值显示为 <0.001

# 构建要显示在图上的标注文本
annot_text <- sprintf(
  "R = %s\nP %s",
  r_value,
  ifelse(grepl("<", p_value), p_value, paste0("= ", p_value))
)

# 3. 绘制图形
p <- ggplot(df, aes(x = log2_PSME3_TPM, y = log2_TUBB3_TPM)) +
  # 绘制散点：shape=21 允许同时设置边框和填充色
  # fill 设置内部颜色，color 设置边框为白色，alpha 设置透明度，stroke 设置边框粗细
  geom_point(
    shape = 21,
    fill = "#3498db",
    color = "white",
    size = 3,
    alpha = 0.7,
    stroke = 0.8
  ) +

  # 添加线性回归拟合线及置信区间 (se = TRUE 默认开启)
  geom_smooth(
    method = "lm",
    color = "#3f5da9ff",
    fill = "gray70",
    linetype = "solid"
  ) +

  # 在图的左上角添加 R 和 P 值标注
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = annot_text,
    hjust = -0.2,
    vjust = 1.5,
    size = 5,
    fontface = "italic"
  ) +

  # 设置坐标轴标签
  labs(
    x = bquote("log"[2] * "(PSME3 TPM)"),
    y = bquote("log"[2] * "(TUBB3_TPM)")
  ) +

  # 设置整体主题风格
  theme_par() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  fix_panel(
    width = 7,
    height = 5
  )

# 4. 显示图形
print(p)

# 如果需要保存为高清 PDF 或 PNG，可以取消下面代码的注释：
ggsave("Correlation_plot_2.svg", plot = p)
