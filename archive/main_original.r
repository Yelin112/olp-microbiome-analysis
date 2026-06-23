# 25-9-21绘制RmSo生长曲线
library(readxl)
library(tidyverse)
library(tidyplots)

# 读取数据
data <- read_excel("data\\RmSoCurve.xlsx", sheet = "Sheet2")

data <- pivot_longer(
  data,
  cols = matches("Rm*|So*"),
  names_to = "Sample",
  values_to = "Value"
) |>
  mutate(Group = ifelse(grepl("Rm", Sample), "Rm", "So"))

pdf("9-21RmSoCurve_2.pdf")
data |>
  tidyplot(x = Time, y = Value, color = Group, dodge_width = 0) |>
  add_curve_fit(orientation = "x") |>
  adjust_x_axis(breaks = seq(0, 15, by = 1)) |>
  adjust_y_axis(title = 'OD600') |>
  apply_tidyplot_theme(MY_LAB_THEME, palette = c('#00A087'))
dev.off()


# 25-9-25 绘制Rm生长曲线  厌氧菌不生长数据用不了
library(readxl)
library(tidyverse)
library(tidyplots)
data <- read_excel("data\\20250923生长曲线4种.xlsx", sheet = "Sheet1")

data <- pivot_longer(
  data,
  cols = matches("Rm*"),
  names_to = "Sample",
  values_to = "Value"
) |>
  mutate(
    Group = ifelse(grepl("Rm", Sample), "Rm", ""),
    Time = factor(Time, levels = sort(as.numeric(as.character(unique(Time))))),
    Value = as.numeric(Value) + 0.15
  )

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


pdf("9-25RmCurve.pdf")
data |>
  tidyplot(x = Time, y = Value, color = Group, dodge_width = 0) |>
  # add_curve_fit(orientation = "x") |>
  add_mean_dot() |>
  # add_line() |>
  add_sem_errorbar(width = 0.6, linewidth = 0.35, ) |>
  adjust_x_axis(breaks = seq(0, 15, by = 1)) |>
  adjust_y_axis(title = 'OD600') |>
  apply_tidyplot_theme(MY_LAB_THEME, palette = c('#3C5488')) |>
  adjust_size(width = 80, height = 60)
dev.off()


# 加载必要的包
library(tidyverse)
library(ggplot2)
library(broom)

# 创建数据框
data <- tibble(
  OD = c(0.102, 0.14, 0.206, 0.246, 0.35),
  CFU = c(3E5, 6E6, 1E7, 1.2E8, 2E8),
  time = c(8, 11, 16, 18, 24)
)

# 第一部分：OD和CFU的关系分析（对CFU进行log10变换）
# 对CFU进行log10变换
data <- data %>%
  mutate(log10_CFU = log10(CFU))

# 执行线性回归（log10(CFU) ~ OD）
lm_model <- lm(log10_CFU ~ OD, data = data)

# 显示回归结果
print("OD和CFU的线性回归结果（log10(CFU) ~ OD）：")
summary(lm_model)

# 获取回归系数
coefficients <- coef(lm_model)
intercept <- coefficients[1]
slope <- coefficients[2]

# 输出回归方程
cat(
  "\n回归方程: log10(CFU) =",
  round(intercept, 4),
  "+",
  round(slope, 4),
  "* OD\n"
)
cat("等价于: CFU = 10^(", round(intercept, 4), "+", round(slope, 4), "* OD)\n")

# 获取回归统计信息
regression_stats <- glance(lm_model)
cat("R-squared:", round(regression_stats$r.squared, 4), "\n")
cat("Adjusted R-squared:", round(regression_stats$adj.r.squared, 4), "\n")

# 可视化1：OD和CFU的关系（对数尺度）
p1 <- ggplot(data, aes(x = OD, y = CFU)) +
  geom_point(size = 3, color = "blue", alpha = 0.7) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    aes(y = 10^stat(y)), # 将预测值转换回原始尺度
    color = "red",
    se = TRUE,
    fill = "lightpink",
    alpha = 0.3
  ) +
  scale_y_continuous(
    trans = 'log10',
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    title = "OD和CFU的关系",
    subtitle = paste0(
      "log10(CFU) = ",
      round(intercept, 4),
      " + ",
      round(slope, 4),
      " * OD"
    ),
    x = "OD值",
    y = "CFU (对数尺度)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# 可视化2：OD和log10(CFU)的线性关系
p2 <- ggplot(data, aes(x = OD, y = log10_CFU)) +
  geom_point(size = 3, color = "blue", alpha = 0.7) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "red",
    se = TRUE,
    fill = "lightpink",
    alpha = 0.3
  ) + # span控制平滑度
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf, # -Inf 和 Inf 代表绘图区的极左和极上
    hjust = -0.1,
    vjust = 1.5, # 微调偏移量，防止文本紧贴着边框
    label = "log[10](CFU) == 4.89 + 10.72 %*% OD", # R的数学表达式语法
    parse = TRUE, # 开启公式解析渲染
    size = 5.5, # 字体大小
    color = "black",
    family = "sans" # 保持与全局字体一致
  ) +
  labs(title = "", x = expression(OD[600]), y = expression(log[10](CFU))) +
  theme_my_stat(base_size = 14, grid_type = "horizontal")
ggsave("RmCurve_CFU.pdf", p2, width = 5, height = 4)


# 第二部分：OD随时间变化的平滑曲线
# 可视化3：OD随时间变化的散点图和平滑曲线
library(ggplot2)
source('my_themes.R')
p3 <- ggplot(data, aes(x = time, y = OD)) +
  geom_point(size = 3, color = "darkgreen", alpha = 0.7) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    color = "orange",
    se = TRUE,
    fill = "lightyellow",
    alpha = 0.3,
    span = 1.2
  ) + # span控制平滑度
  labs(title = "OD随时间的变化", x = "时间", y = "OD值") +
  theme_my_stat(base_size = 14, grid_type = "horizontal")

# OLP基线表 -----------------------------------------------------------------

library(readxl)
library(dplyr)
data_HC <- read_excel(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\患者管理\\Clinics_Normal.xlsx"
)
data_OLP <- read_excel(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\患者管理\\Clinics_OLP更新.xlsx"
)
colnames(data_HC)
colnames(data_OLP)


data_all <- bind_rows(data_HC, data_OLP) |>
  select(all_of(intersect(colnames(data_HC), colnames(data_OLP)))) |>
  select(-c(1:4, 6))
library(tableone)


# 创建统计表
vars <- c(
  "Type",
  "Age",
  "Sex",
  "Smoke",
  "Alcohol",
  "Betel chewing",
  "Periodontal disease",
  "History of infection"
) # 要比较的变量
categorical_vars <- c(
  "Type",
  "Sex",
  "Smoke",
  "Alcohol",
  "Betel chewing",
  "Periodontal disease",
  "History of infection"
) # 分类变量

# 生成TableOne对象
table1 <- CreateTableOne(
  vars = vars,
  strata = "Group", # 分组变量
  data = data_all, # 您的数据框
  factorVars = categorical_vars,
  addOverall = TRUE
)

# 打印结果
print(
  table1,
  showAllLevels = TRUE,
  missing = TRUE, # 显示缺失值情况
  test = TRUE, # 进行组间比较检验
  formatOptions = list(big.mark = ","),
  catDigits = 1, # 分类变量小数位数
  contDigits = 1
) # 连续变量小数位数

write.csv(
  print(
    table1,
    showAllLevels = TRUE,
    missing = TRUE, # 显示缺失值情况
    test = TRUE, # 进行组间比较检验
    formatOptions = list(big.mark = ","),
    catDigits = 1, # 分类变量小数位数
    contDigits = 1
  ), # 连续变量小数位,
  "tableone.csv"
)

# 如果需要非参数检验（如果年龄不服从正态分布）
print(
  table1,
  nonnormal = "年龄", # 指定年龄使用非参数检验
  showAllLevels = TRUE,
  formatOptions = list(big.mark = ",")
)


colnames(data_OLP) |> dput()

vars <- c(
  "Age",
  "Sex",
  "Type",
  "Scores",
  "Sampling site",
  "Course(M)",
  "Probing hemorrhage",
  "Gingival index",
  "Smoke",
  "Alcohol",
  "Betel chewing",
  "Periodontal disease",
  "History of infection"
) # 要比较的变量
categorical_vars <- c(
  "Sex",
  "Type",
  "Sampling site",
  "Probing hemorrhage",
  "Smoke",
  "Alcohol",
  "Betel chewing",
  "Periodontal disease",
  "History of infection"
) # 分类变量

# 生成TableOne对象
table1 <- CreateTableOne(
  vars = vars,
  strata = "Type", # 分组变量
  data = data_OLP, # 您的数据框
  factorVars = categorical_vars,
  addOverall = TRUE
)

# 打印结果
print(
  table1,
  showAllLevels = TRUE,
  missing = TRUE, # 显示缺失值情况
  test = TRUE, # 进行组间比较检验
  formatOptions = list(big.mark = ","),
  catDigits = 1, # 分类变量小数位数
  contDigits = 1
) # 连续变量小数位数


write.csv(
  print(
    table1,
    showAllLevels = TRUE,
    missing = TRUE, # 显示缺失值情况
    test = TRUE, # 进行组间比较检验
    formatOptions = list(big.mark = ","),
    catDigits = 1, # 分类变量小数位数
    contDigits = 1
  ), # 连续变量小数位,
  "tableone_OLP.csv"
)


# qPCR数据美化 ---------------------------------------------------------------

data <- readxl::read_xlsx(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\实验数据\\qPCR\\qPCR7_8\\PCR.xlsx",
  sheet = 2
)

library(readxl)
library(ggplot2)
library(patchwork) # 用于组合多个图

source('Scripts\\my_themes.R')
# source('Scripts\\compare.R')
source('Scripts\\compare_plot_optimized.R')


palette_two <- list(Group = c(OLP = "#E88471FF", HC = "#39B185FF"))
genes <- unique(data$Gene)

# 3. 循环绘制每个基因
plot_list <- list()

for (gene in genes) {
  cat("正在绘制:", gene, "\n")

  # 筛选当前基因数据
  gene_data <- data[data$Gene == gene, ]

  # 绘图
  p <- compare_plot(
    data = gene_data,
    value.var = "RelativeExp",
    group.by = "Group",
    fill.by = "Group",
    plot_type = "bar",
    add_box = TRUE,
    add_point = TRUE,
    add_stat = "t.test",
    comparisons = list(c("Rm", "Control")),
    palette = c("#E88471FF", "#39B185FF"),
    xlab = "Group",
    ylab = "Relative Expression",
    title = gene,
    theme_use = theme_my_stat
  )
  p <- p +
    fix_panel(width = 4, height = 4) +
    theme(
      legend.position = "none",
      text = element_text(face = "bold") # 全局加粗
    )
  # 保存到列表
  plot_list[[gene]] <- p
}
plot_list[[1]]

# 创建输出文件夹
if (!dir.exists("qPCR")) {
  dir.create("qPCR")
}

# 批量保存为PDF
for (gene in names(plot_list)) {
  cat("保存:", gene, "\n")

  ggsave(
    filename = file.path("qPCR", paste0(gene, ".pdf")),
    plot = plot_list[[gene]],
    width = 4,
    height = 4
  )
}


# 小鼠模型qPCR定植 -------------------------------------------------------------
data <- readxl::read_xlsx(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\实验数据\\qPCR\\qPCR_Rm定植\\Rm定植.xlsx",
  sheet = 1
)
p <- compare_plot(
  data = data,
  value.var = "Value",
  group.by = "Group",
  fill.by = "Group",
  plot_type = "bar",
  add_point = TRUE,
  add_stat = "t.test",
  comparisons = list(c("Rm", "Control")),
  palette = rm_colors_two,
  xlab = "Group",
  ylab = "Relative Abundance of R.m",
  title = '',
  theme_use = theme_my_stat
) +
  fix_panel(width = 4, height = 5) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold") # 全局加粗
  )
ggsave("qPCR//Rm模型定植.pdf", plot = p, width = 6, height = 6)


# FISH荧光统计 ---------------------------------------------------------------

olp_colors_two <- c(OLP = "#F4A58295", HC = "#4393C395")
data <- readxl::read_xlsx(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\FISH\\FISH_servbio.xlsx"
)
p <- compare_plot(
  data = data,
  value.var = "Area_Ratio",
  group.by = "Group",
  fill.by = "Group",
  plot_type = "bar",
  add_point = TRUE,
  add_stat = "t.test",
  comparisons = list(c("OLP", "HC")),
  palette = olp_colors_two,
  xlab = "Group",
  ylab = "Relative Bacterial Load",
  title = 'FISH_N=6',
  theme_use = theme_my_stat
) +
  fix_panel(width = 4, height = 5) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold") # 全局加粗
  )
ggsave("Fish//FISH_Rm共定位.pdf", plot = p, width = 6, height = 6)


p <- compare_plot(
  data = data,
  value.var = "Intensity_Norm",
  group.by = "Group",
  fill.by = "Group",
  plot_type = "bar",
  add_point = TRUE,
  add_stat = "t.test",
  comparisons = list(c("OLP", "HC")),
  palette = olp_colors_two,
  xlab = "Group",
  ylab = "Normalized Fluorescence Intensity",
  title = 'FISH_N=6',
  theme_use = theme_my_stat
) +
  fix_panel(width = 4, height = 5) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold") # 全局加粗
  )
ggsave("Fish//FISH_Rm共定位_荧光强度.pdf", plot = p, width = 6, height = 6)


rm_colors_two <- c(Rm = "#E88471FF", Control = "#39B185FF")

data <- readxl::read_xlsx(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\FISH\\Model\\Mouse_model.xlsx"
)
p <- compare_plot(
  data = data,
  value.var = "adhesion_index",
  group.by = "Group",
  fill.by = "Group",
  plot_type = "bar",
  add_point = TRUE,
  add_stat = "t.test",
  comparisons = list(c("Rm", "Control")),
  palette = rm_colors_two,
  xlab = "Group",
  ylab = "Relative Bacterial Adhesion",
  title = 'Quantitative analysis of Rm adhesion',
  theme_use = theme_my_stat
) +
  fix_panel(width = 4, height = 5) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold") # 全局加粗
  )
ggsave("Fish//Model_FISH_Rm共定位.pdf", plot = p, width = 6, height = 6)


# 免疫组化统计 -----------------------------------------------------------------

data <- readxl::read_xlsx(
  "IHC\\CD3_IHC.xlsx"
)
rm_colors_two <- c(Rm = "#E88471FF", Control = "#39B185FF")


p <- compare_plot(
  data = data,
  value.var = "Positive_Ratio",
  group.by = "Group",
  fill.by = "Group",
  plot_type = "bar",
  add_point = TRUE,
  add_stat = "t.test",
  comparisons = list(c("Rm", "Control")),
  palette = rm_colors_two,
  xlab = "Group",
  ylab = "CD3 Positive Ratio",
  title = 'IHC_N=5',
  theme_use = theme_my_stat
) +
  fix_panel(width = 4, height = 5) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold") # 全局加粗
  )
ggsave("IHC//IHC_CD3_pos.pdf", plot = p, width = 6, height = 6)


# 宏基因组网络分析 ---------------------------------------------------------------

source('Scripts\\filter_abundance.R')
library(tidyverse)
metadata <- read.csv(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\metagenome\\metadata.csv",
  row.names = 1
)

library(readxl)
bac_species_profile <- read.delim(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\metagenome\\data\\bac\\tax\\species.profile.xls",
  header = TRUE,
  row.names = 1
)
bac_spe <- bac_species_profile[, intersect(
  colnames(bac_species_profile),
  metadata$ID_PBS
)]
bac_spe <- filter_for_network(bac_spe) |> pluck('filtered')

tax_tab <- read.delim(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\metagenome\\data\\bac\\tax\\species.tax.txt",
  header = TRUE
)

bac_spe_OLP <- bac_spe[,
  colnames(bac_spe) %in% metadata[which(metadata$Group == "OLP"), ][['ID_PBS']]
]
bac_spe_HC <- bac_spe[,
  colnames(bac_spe) %in% metadata[which(metadata$Group == "HC"), ][['ID_PBS']]
]


library(ggNetView)
graph_obj <- build_graph_from_mat(
  mat = bac_spe,
  node_annotation = tax_tab,
  transfrom.method = "none",
  method = "WGCNA",
  proc = "BH",
  r.threshold = 0.7,
  p.threshold = 0.05,
  module.method = "Fast_greedy",
  top_modules = 15,
  seed = 1115
)

graph_obj_olp <- build_graph_from_mat(
  mat = bac_spe_OLP,
  node_annotation = tax_tab,
  transfrom.method = "none",
  method = "WGCNA",
  proc = "BH",
  r.threshold = 0.7,
  p.threshold = 0.05,
  module.method = "Fast_greedy",
  top_modules = 15,
  seed = 1115
)


graph_obj_HC <- build_graph_from_mat(
  mat = bac_spe_HC,
  node_annotation = tax_tab,
  transfrom.method = "none",
  method = "WGCNA",
  proc = "BH",
  r.threshold = 0.7,
  p.threshold = 0.05,
  module.method = "Fast_greedy",
  top_modules = 15,
  seed = 1115
)


library(igraph)

# 1. 提取所有节点的名称向量
all_node_names <- V(graph_obj)$name

# 2. 寻找精确名称的行号
target_idx <- which(all_node_names == "Rothia_mucilaginosa")
print(target_idx)

# # 为特定细菌创建高亮标记
# V(graph_obj)$highlight <- ifelse(
#   V(graph_obj)$name == "Rothia_mucilaginosa",
#   "Rothia_mucilaginosa",
#   "Other"
# )

# # 假设节点名称存储在 "name" 属性中
# V(graph_obj)$label_show <- ifelse(
#   V(graph_obj)$name == "Rothia_mucilaginosa",
#   "Rothia_mucilaginosa",
#   ""
# )

# # 然后使用 pointlabel 参数
# p <- ggNetView(graph_obj,
#                layout = "gephi",
#                pointlabel = "label_show",  # 使用创建的标签列
#                pointlabelsize = 5)

# # 使用 color.by 参数高亮
# p <- ggNetView(graph_obj,
#                layout = "gephi",
#                color.by = "highlight",
#                color = c("Rothia_mucilaginosa" = "red", "Other" = "black"),
#                pointsize = c(3, 15))  # 可以设置更大的点

# 1. 基础网络图
p1 <- ggNetView(
  graph_obj = graph_obj_HC,
  layout = "petal",
  center = FALSE,
  shrink = 0.65,
  layout.module = "adjacent",
  group.by = "Modularity",
  k_nn = 70,
  # center = TRUE,
  # pointlabel = "top2",
  idx = 19,
  # add_outer = TRUE,
  mapping_line = TRUE, # 开启连线映射 (自动区分正负相关)
  curve = TRUE, # 使用弧线
  linealpha = 0.4,
  fill.by = "Modularity"
)

# 2. 基础网络图 (添加外圈和标签)
p1_2 <- ggNetView(
  graph_obj = graph_obj,
  layout = "gephi",
  center = FALSE,
  shrink = 0.85,
  layout.module = "adjacent",
  group.by = "Modularity",
  fill.by = "Modularity",
  k_nn = 50,
  add_outer = TRUE,
  label = TRUE
)

# --- 定义统一的自定义颜色面板 ---
my_colors <- c(
  '#006d2c',
  '#810f7c',
  '#0868ac',
  '#d7301f',
  '#02818a',
  '#e7298a',
  '#f768a1',
  '#41ab5d',
  '#1d91c0',
  '#6a3d9a',
  '#ffff99',
  '#b15928',
  '#8dd3c7',
  '#ffffb3',
  '#bebada',
  '#bdbdbd'
)

# 3. 自定义颜色网络图
p2 <- ggNetView(
  graph_obj = graph_obj,
  layout = "gephi",
  center = FALSE,
  shrink = 0.85,
  layout.module = "adjacent",
  group.by = "Modularity",
  fill.by = "Modularity",
  color = my_colors,
  fill = my_colors
)

# 4. 自定义颜色网络图 (添加外圈和标签)
p2_2 <- ggNetView(
  graph_obj = graph_obj,
  layout = "gephi",
  center = FALSE,
  shrink = 0.85,
  layout.module = "adjacent",
  group.by = "Modularity",
  fill.by = "Modularity",
  color = my_colors,
  fill = my_colors,
  add_outer = TRUE,
  label = TRUE
)


source('Scripts\\NetworkAnalyzer.R')
bac_spe_OLP <- bac_spe[,
  colnames(bac_spe) %in% metadata[which(metadata$Group == "OLP"), ][['ID_PBS']]
]
bac_spe_HC <- bac_spe[,
  colnames(bac_spe) %in% metadata[which(metadata$Group == "HC"), ][['ID_PBS']]
]

net_olp <- NetworkAnalyzer$new(
  data = bac_spe_OLP,
  data_type = "abundance",
  filter_threshold = 0.001 # 过滤低丰度特征
)

net_olp$build_network(
  method = "correlation",
  cor_method = "spearman",
  cor_threshold = 0.6,
  use_abs = TRUE,
  keep_negative = TRUE,
  p_threshold = 0.05
)

# 检测模块
net_olp$detect_modules(method = "louvain", use_abs_weight = TRUE)

# 计算网络属性
net_olp$calculate_properties(calculate_roles = TRUE)

# 查看结果
print(net_olp)


net_HC <- NetworkAnalyzer$new(
  data = bac_spe_HC,
  data_type = "abundance",
  filter_threshold = 0.001 # 过滤低丰度特征
)

net_HC$build_network(
  method = "correlation",
  cor_method = "spearman",
  cor_threshold = 0.6,
  p_threshold = 0.05
)

# 检测模块
net_HC$detect_modules(method = "louvain")

# 计算网络属性
net_HC$calculate_properties(calculate_roles = TRUE)

# 查看结果
print(net_HC)

comparison <- net_olp$compare_with(net_HC, permutations = 1000)

# 查看比较结果
print(comparison)

# 可视化比较结果
plot(comparison, type = "properties")
plot(comparison, type = "roles")
plot(comparison, type = "degree")


library(ggraph)
library(ggplot2)

p <- net_olp$plot_network(
  method = "ggraph",
  layout = "fr",
  color_by = "module",
  size_by = "degree"
)

# 进一步自定义
p +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "Microbial Network",
    subtitle = paste(
      vcount(net$network),
      "nodes,",
      ecount(net$network),
      "edges"
    )
  ) +
  theme(legend.position = "bottom")


# 细菌lefse重绘 --------------------------------------------------------------

# 加载必要的包
library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# 1. 录入所有数据
# ==========================================
aligned <- data.frame(
  Taxon = c(
    "R. mucilaginosa",
    "T. forsythia",
    "P. intermedia",
    "P. jejuni",
    "F. pseudoperiodonticum",
    "V. atypica",
    "Nanosynbacter sp.",
    "L. parvula",
    "S. oralis",
    "S. mitis"
  ),
  Traditional = c(4.3, 3.2, 3.52, 3.4, 3.48, 3.1, 2.8, 3.0, -4.0, -4.3),
  MAGs = c(4.4, 3.8, 3.6, 3.3, 3.1, 2.6, 3.0, 2.3, -3.3, -3.6),
  Category = c("Target", rep("Shared", 9))
)

unaligned_trad <- data.frame(
  Taxon = c(
    "S. odontolytica",
    "N. mucosa",
    "P. gingivalis",
    "P. endodontalis",
    "P. multiformis",
    "P. denticola",
    "P. veroralis",
    "S. meyeri",
    "M. diversum",
    "F. alocis",
    "P. micra",
    "Prevotella_sp_299",
    "C. concisus",
    "Leptotrichia_sp_221",
    "L. umeaense",
    "O. uli",
    "S. cardiffensis",
    "M. paludicola",
    "C. turicensis",
    "D. oralis",
    "C. aerofaciens",
    "S. loihica",
    "Erwinia_sp.",
    "P. scopos",
    "Nanosynbacter_sp.",
    "Streptococcus_sp_192",
    "F. periodonticum",
    "B. punctulatus",
    "P. rubescens",
    "Pantoea_sp.",
    "Streptococcus_sp_1643",
    "K. oralis",
    "P. histicola",
    "Streptococcus_sp_116",
    "S. pseudopneumoniae",
    "S. pneumoniae"
  ),
  Traditional = c(
    4.1,
    3.8,
    3.5,
    3.4,
    3.3,
    3.3,
    3.0,
    2.8,
    2.8,
    2.8,
    2.8,
    2.7,
    2.7,
    2.5,
    2.5,
    2.3,
    2.3,
    2.2,
    2.2,
    2.2,
    2.1,
    2.1,
    2.1,
    2.1,
    2.1,
    2.1,
    2.1,
    -2.0,
    -2.1,
    -2.3,
    -2.3,
    -2.7,
    -3.0,
    -3.1,
    -3.3,
    -4.0
  ),
  MAGs = NA,
  Category = "Unaligned"
)

unaligned_mags <- data.frame(
  Taxon = c(
    "Pauljensenia_sp_725",
    "Pauljensenia_sp_545",
    "Neisseria_sp_165",
    "A. graevenitzii",
    "S. salivarius",
    "S. parasanguinis_F",
    "Pauljensenia_sp_dur",
    "S. rubneri",
    "S. longum_A",
    "Nanogingivalis_sp_795",
    "Lachnoanaerobaculum_sp_675",
    "S. parasanguinis",
    "P. salivae",
    "Pauljensenia_sp_hk",
    "C. bilenii",
    "C. ochracea",
    "G. adiacens",
    "A. massiliensis",
    "A. gerencseriae",
    "V. parvula_A",
    "H. haemolyticus",
    "N. subflava_A",
    "Unclassified"
  ),
  Traditional = NA,
  MAGs = c(
    4.1,
    3.75,
    3.7,
    3.65,
    3.4,
    3.1,
    3.0,
    2.8,
    2.8,
    2.4,
    2.1,
    -2.6,
    -2.7,
    -2.9,
    -3.1,
    -3.1,
    -3.2,
    -3.2,
    -3.3,
    -3.3,
    -3.6,
    -3.8,
    -4.3
  ),
  Category = "Unaligned"
)

all_data <- bind_rows(aligned, unaligned_trad, unaligned_mags)

# ==========================================
# 2. 为所有元素分配等距坐标 (保持防重叠)
# ==========================================
all_data$Y_left <- NA
all_data$Y_right <- NA

# 左侧间距 (保留跨度以分散所有背景噪点)
idx_t_pos <- which(all_data$Traditional > 0)
all_data$Y_left[idx_t_pos[order(all_data$Traditional[idx_t_pos])]] <- seq(
  0.5,
  12,
  length.out = length(idx_t_pos)
)
idx_t_neg <- which(all_data$Traditional < 0)
all_data$Y_left[idx_t_neg[order(all_data$Traditional[idx_t_neg])]] <- seq(
  -12,
  -0.5,
  length.out = length(idx_t_neg)
)

# 右侧间距
idx_m_pos <- which(all_data$MAGs > 0)
all_data$Y_right[idx_m_pos[order(all_data$MAGs[idx_m_pos])]] <- seq(
  0.5,
  12,
  length.out = length(idx_m_pos)
)
idx_m_neg <- which(all_data$MAGs < 0)
all_data$Y_right[idx_m_neg[order(all_data$MAGs[idx_m_neg])]] <- seq(
  -12,
  -0.5,
  length.out = length(idx_m_neg)
)

# ==========================================
# 3. 提取绘图数据与文字排版参数
# ==========================================
left_data <- all_data %>%
  filter(!is.na(Traditional)) %>%
  mutate(
    hjust_val = ifelse(Traditional > 0, 0, 1),
    nudge_x = ifelse(Traditional > 0, 0.2, -0.2),
    t_size = case_when(
      Category == "Target" ~ 5,
      Category == "Shared" ~ 4,
      TRUE ~ 3
    ),
    t_face = case_when(
      Category == "Target" ~ "bold.italic",
      Category == "Shared" ~ "bold.italic",
      TRUE ~ "italic"
    )
  )

right_data <- all_data %>%
  filter(!is.na(MAGs)) %>%
  mutate(
    hjust_val = ifelse(MAGs > 0, 0, 1),
    nudge_x = ifelse(MAGs > 0, 0.2, -0.2),
    t_size = case_when(
      Category == "Target" ~ 5,
      Category == "Shared" ~ 4,
      TRUE ~ 3
    ),
    t_face = case_when(
      Category == "Target" ~ "bold.italic",
      Category == "Shared" ~ "bold.italic",
      TRUE ~ "italic"
    )
  )

# 配色方案 (背景色调浅至 grey80，使其安静地做背景)
my_colors <- c(
  "Target" = "#D32F2F",
  "Shared" = "grey20",
  "Unaligned" = "grey80"
)
my_theme <- theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none"
  )

# ==========================================
# 4. 绘制左侧：Traditional 棒棒图
# ==========================================
p_left <- ggplot(left_data, aes(y = Y_left)) +
  geom_segment(aes(
    x = 0,
    xend = Traditional,
    yend = Y_left,
    color = Category,
    size = ifelse(
      Category == "Target",
      1.8,
      ifelse(Category == "Shared", 1, 0.5)
    ),
    alpha = ifelse(
      Category == "Target",
      1,
      ifelse(Category == "Shared", 0.9, 0.5)
    )
  )) +
  geom_point(aes(
    x = Traditional,
    color = Category,
    size = ifelse(
      Category == "Target",
      4,
      ifelse(Category == "Shared", 2.5, 1.2)
    ),
    alpha = ifelse(
      Category == "Target",
      1,
      ifelse(Category == "Shared", 0.9, 0.5)
    )
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  # 【关键修改】：仅对 Target 和 Shared 渲染文字
  geom_text(
    data = filter(left_data, Category %in% c("Target", "Shared")),
    aes(
      x = Traditional + nudge_x,
      label = Taxon,
      color = Category,
      size = t_size,
      fontface = t_face,
      hjust = hjust_val
    )
  ) +
  scale_color_manual(values = my_colors) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-7.5, 7.5)) +
  scale_y_continuous(limits = c(-12.5, 12.5), expand = c(0, 0)) +
  my_theme +
  labs(title = "Traditional (16S/Meta)") +
  fix_panel(8, 12)

ggsave(
  "lefse//left_Plot.pdf",
  plot = p_left,
  width = 14,
  height = 14,
  dpi = 600
)
# ==========================================
# 5. 绘制中间：排位交互 Slopegraph
# ==========================================
p_mid <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(
    data = filter(all_data, Category == "Unaligned" & !is.na(Y_left)),
    aes(x = 1, y = Y_left),
    color = my_colors["Unaligned"],
    size = 1,
    alpha = 0.4
  ) +
  geom_point(
    data = filter(all_data, Category == "Unaligned" & !is.na(Y_right)),
    aes(x = 2, y = Y_right),
    color = my_colors["Unaligned"],
    size = 1,
    alpha = 0.4
  ) +
  geom_segment(
    data = filter(all_data, Category %in% c("Target", "Shared")),
    aes(
      x = 1,
      xend = 2,
      y = Y_left,
      yend = Y_right,
      color = Category,
      size = ifelse(Category == "Target", 2, 0.8)
    )
  ) +
  geom_point(
    data = filter(all_data, Category %in% c("Target", "Shared")),
    aes(
      x = 1,
      y = Y_left,
      color = Category,
      size = ifelse(Category == "Target", 4, 2.5)
    )
  ) +
  geom_point(
    data = filter(all_data, Category %in% c("Target", "Shared")),
    aes(
      x = 2,
      y = Y_right,
      color = Category,
      size = ifelse(Category == "Target", 4, 2.5)
    )
  ) +
  scale_color_manual(values = my_colors) +
  scale_size_identity() +
  scale_x_continuous(
    limits = c(0.8, 2.2),
    breaks = c(1, 2),
    labels = c("Trad.", "MAGs")
  ) +
  scale_y_continuous(limits = c(-12.5, 12.5), expand = c(0, 0)) +
  my_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(color = "white")
  ) +
  labs(title = "X") +
  fix_panel(3, 12)
ggsave("lefse//mid_Plot.pdf", plot = p_mid, width = 14, height = 14, dpi = 600)
# ==========================================
# 6. 绘制右侧：MAGs 棒棒图
# ==========================================
p_right <- ggplot(right_data, aes(y = Y_right)) +
  geom_segment(aes(
    x = 0,
    xend = MAGs,
    yend = Y_right,
    color = Category,
    size = ifelse(
      Category == "Target",
      1.8,
      ifelse(Category == "Shared", 1, 0.5)
    ),
    alpha = ifelse(
      Category == "Target",
      1,
      ifelse(Category == "Shared", 0.9, 0.5)
    )
  )) +
  geom_point(aes(
    x = MAGs,
    color = Category,
    size = ifelse(
      Category == "Target",
      4,
      ifelse(Category == "Shared", 2.5, 1.2)
    ),
    alpha = ifelse(
      Category == "Target",
      1,
      ifelse(Category == "Shared", 0.9, 0.5)
    )
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  # 【关键修改】：仅对 Target 和 Shared 渲染文字
  geom_text(
    data = filter(right_data, Category %in% c("Target", "Shared")),
    aes(
      x = MAGs + nudge_x,
      label = Taxon,
      color = Category,
      size = t_size,
      fontface = t_face,
      hjust = hjust_val
    )
  ) +
  scale_color_manual(values = my_colors) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-7.5, 7.5)) +
  scale_y_continuous(limits = c(-12.5, 12.5), expand = c(0, 0)) +
  my_theme +
  labs(title = "MAGs Reconstruction") +
  fix_panel(8, 12)
ggsave(
  "lefse//right_Plot.pdf",
  plot = p_right,
  width = 14,
  height = 14,
  dpi = 600
)
# ==========================================
# 7. 拼接并输出
# ==========================================
final_plot <- p_left + p_mid + p_right + plot_layout(widths = c(3, 1, 3))

final_plot <- final_plot +
  plot_annotation(
    title = "Comprehensive Cross-Validation of Microbial Biomarkers",
    subtitle = "Background lines indicate unaligned taxa. R. mucilaginosa remains the robust consensus target.",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(
        size = 13,
        hjust = 0.5,
        color = "grey30",
        margin = margin(b = 10)
      )
    )
  )

print(final_plot)

# 保存代码示例（推荐尺寸，保证留白与美观）
ggsave(
  "lefse//Composite_Plot.pdf",
  plot = final_plot,
  width = 14,
  height = 14,
  dpi = 600
)


# 韦恩图 --------------------------------------------------------------------

# 加载必要的包
library(ggplot2)
library(ggvenn)

# ==========================================
# 1. 使用您的统计数字生成模拟基因列表
# ==========================================
# 这样可以欺骗算法，让它直接按照您的数字画出完美的对称圆
genes_olp_unique <- paste0("OLP_", 1:5455)
genes_hc_unique <- paste0("HC_", 1:5473)
genes_shared <- paste0("Shared_", 1:11289)

# 构建两组的完整集合
data_list <- list(
  OLP = c(genes_olp_unique, genes_shared),
  HC = c(genes_hc_unique, genes_shared)
)
olp_colors_two <- c(OLP = "#F4A58295", HC = "#4393C395")
# ==========================================
# 2. 开始绘制对称美观的韦恩图
# ==========================================
# 继续使用我们之前设定的 OLP (砖红/亮红) 和 HC (森林绿) 配色体系
p <- ggvenn(
  data_list,
  columns = c("OLP", "HC"),

  # 颜色与透明度设置
  fill_color = olp_colors_two,
  fill_alpha = 0.65,

  # 边框美化：使用白色细边框代替默认的黑色粗边框，质感瞬间提升
  stroke_color = "white",
  stroke_size = 1.5,
  stroke_alpha = 1,

  # 文字排版
  set_name_size = 7, # 组名 (OLP, HC) 的字号
  text_size = 6, # 圆圈内部数字的字号
  text_color = "black", # 数字颜色

  # 隐藏无用的百分比显示（可选，如果想显示改为TRUE即可）
  show_percentage = T
)

# 打印图像
print(p)

# ==========================================
# 3. 导出为高清 PDF 供 AI 后期排版
# ==========================================
# ggsave("Symmetrical_Venn.pdf", plot = p, width = 6, height = 6, dpi = 300)

# 宏转录组LEfse可视化 -----------------------------------------------------------

# 加载必要的包
library(ggplot2)
library(dplyr)

# 1. 读取并清洗数据
file_path <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/lefse/宏转录组/Species_tpm_out.xls"
lefse_data <- read.table(
  file_path,
  sep = "\t",
  header = T,
  stringsAsFactors = FALSE,
  quote = ""
)
colnames(lefse_data) <- c("Taxon", "LogMaxMean", "Group", "LDA", "Pvalue")

clean_data <- lefse_data %>%
  filter(!is.na(LDA) & Group != "")

# 提取各组 Top 15 (为避免环形图过于拥挤，Top 15 是非常理想的数量)
top_data <- clean_data %>%
  group_by(Group) %>%
  top_n(15, wt = LDA) %>%
  ungroup() %>%
  arrange(Group, LDA)

# 2. 为环形图构建极坐标参数 (为了在两组之间留出缺口)
empty_bar <- 3 # 组与组之间留出 3 个空位的间距
to_add <- data.frame(matrix(
  NA,
  empty_bar * nlevels(as.factor(top_data$Group)),
  ncol(top_data)
))
colnames(to_add) <- colnames(top_data)
to_add$Group <- rep(levels(as.factor(top_data$Group)), each = empty_bar)

# 合并数据并重置行号
plot_data_polar <- rbind(top_data, to_add)
plot_data_polar <- plot_data_polar %>% arrange(Group)
plot_data_polar$id <- seq(1, nrow(plot_data_polar))

# 3. 计算文字标签的精确旋转角度
# 根据条形图所在的极坐标角度，智能翻转文字，确保所有文字均可正向阅读
angle <- 90 - 360 * (plot_data_polar$id - 0.5) / nrow(plot_data_polar)
plot_data_polar$hjust <- ifelse(angle < -90, 1, 0)
plot_data_polar$angle <- ifelse(angle < -90, angle + 180, angle)

# 4. 配色方案
color_hc <- "#228B22" # 森林绿
color_olp <- "#D32F2F" # 砖红

# 5. 绘制高级环形图
p_circular <- ggplot(
  plot_data_polar,
  aes(x = as.factor(id), y = LDA, fill = Group)
) +
  # 柱状图基础层
  geom_bar(stat = "identity", alpha = 0.85, color = "white", linewidth = 0.5) +

  # 设置 Y 轴范围：为了中心留空，将下限设为负数 (如 -3)
  scale_y_continuous(
    limits = c(-3, max(plot_data_polar$LDA, na.rm = TRUE) + 2)
  ) +
  scale_fill_manual(
    values = c("HC" = color_hc, "OLP" = color_olp),
    na.translate = FALSE
  ) +

  # 极坐标转换
  coord_polar(start = 0) +

  # 添加精确角度的文字标签
  geom_text(
    data = plot_data_polar,
    aes(x = id, y = LDA + 0.2, label = Taxon, hjust = hjust),
    color = "black",
    fontface = "italic",
    alpha = 0.9,
    size = 3.5,
    angle = plot_data_polar$angle,
    inherit.aes = FALSE
  ) +

  # 极简主题消除背景
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p_circular)

# 导出建议 (保存为正方形最佳)
# ggsave("Circular_Barplot_Top15.pdf", plot = p_circular, width = 8, height = 8, dpi = 300)

# 接上面的 clean_data 与 top_n 提取步骤
# 为了实现双向对称，将 HC 组的 LDA 设为负数
top_data_lollipop <- clean_data %>%
  group_by(Group) %>%
  top_n(15, wt = LDA) %>%
  ungroup() %>%
  mutate(Plot_LDA = ifelse(Group == "HC", -LDA, LDA)) %>%
  arrange(Plot_LDA) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))

# 绘制极简双向棒棒糖图
p_lollipop <- ggplot(
  top_data_lollipop,
  aes(x = Taxon, y = Plot_LDA, color = Group)
) +

  # 画背景虚线基准面
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.8
  ) +

  # 画棒棒糖的“棍子”
  geom_segment(
    aes(x = Taxon, xend = Taxon, y = 0, yend = Plot_LDA),
    linewidth = 1.2,
    alpha = 0.8
  ) +

  # 画棒棒糖的“糖” (数据点)
  geom_point(size = 4) +

  # 翻转坐标系
  coord_flip() +

  # 设定颜色体系
  scale_color_manual(values = c("HC" = color_hc, "OLP" = color_olp)) +
  scale_y_continuous(limits = c(-5.5, 5.5), breaks = seq(-4, 4, by = 2)) +

  labs(
    title = "Top Active Microbial Drivers in OLP vs HC",
    x = NULL,
    y = "LDA Score (log10)"
  ) +

  # 定制主题，营造干净利落的版面风格
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(), # 去除横向网格
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey80"),

    axis.text.y = element_text(size = 11, face = "italic", color = "black"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(
      size = 13,
      face = "bold",
      margin = margin(t = 15)
    ),

    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 20)
    ),

    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold")
  )

print(p_lollipop)

# 导出建议
ggsave(
  "Lollipop_Chart_Top15.pdf",
  plot = p_lollipop,
  width = 9,
  height = 8,
  dpi = 300
)


# 溶菌酶可视化 -----------------------------------------------------------------

# 1. 读取数据 (请确保文件名与您的工作目录一致)
# 使用 read_csv 可以更好地保留列名中的空格
# 如果尚未安装这些包，请取消注释并运行以下代码：
# install.packages(c("tidyverse", "ggplot2", "igraph", "ggraph", "ggsci"))

library(tidyverse)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggsci) # 用于科研配色
df <- read_csv(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\公司标准分析报告\\重组\\Rm溶菌酶\\宿主信息.csv"
)

# 2. 数据清洗与分组
df_clean <- df %>%
  # 使用正则表达式提取种水平注释 (s__后面的内容)
  mutate(
    Species = str_extract(`Host taxonomy`, "s__[^;]+"),
    Species = str_remove(Species, "s__"),
    # 如果没有种水平注释(NA)，则标记为仅属水平
    Species = ifelse(is.na(Species), "Rothia (Genus level only)", Species),

    # 划分核心高亮组：目标细菌 vs 其他Rothia
    Target_Group = case_when(
      str_detect(Species, "Rothia mucilaginosa") ~ "R. mucilaginosa",
      TRUE ~ "Other Rothia species"
    )
  ) %>%
  # 按照噬菌体和Rank排序，每个vOTU只保留得分最高(Rank 1)的一个预测结果
  group_by(Virus) %>%
  arrange(Rank) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  # 将核心靶点放在因子图例的前面
  mutate(
    Target_Group = factor(
      Target_Group,
      levels = c("R. mucilaginosa", "Other Rothia species")
    )
  )


# 绘制散点图
p_dot <- ggplot(
  df_clean,
  aes(x = reorder(Virus, as.numeric(Target_Group)), y = Species)
) +
  # 绘制气泡：颜色代表是否为靶标细菌，形状代表预测Method
  geom_point(
    aes(color = Target_Group, shape = Method),
    size = 4,
    alpha = 0.85
  ) +
  # 自定义配色：R. mucilaginosa 用显眼的珊瑚红/橙色，其他用柔和蓝色
  scale_color_manual(
    values = c(
      "R. mucilaginosa" = "#FF7043",
      "Other Rothia species" = "#64B5F6"
    )
  ) +
  # 美化主题
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10, face = "italic"), # 细菌名使用斜体
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "grey90"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "Phage vOTUs (n = 21)",
    y = "Predicted Host Species",
    title = "Phage-Host Predictions in Genus Rothia",
    color = "Host Category",
    shape = "Prediction Method"
  )

# 显示并保存图片
print(p_dot)
ggsave("Phage_Host_DotPlot.svg", p_dot, width = 8, height = 4)


# 1. 构建网络所需的边列表 (Edge list)
edges <- df_clean %>%
  select(Virus, Species) %>%
  rename(from = Virus, to = Species)

# 2. 构建节点属性 (Node attributes)
# 提取所有独特的噬菌体和宿主节点
nodes <- data.frame(name = unique(c(edges$from, edges$to))) %>%
  mutate(
    # 定义节点类型
    Type = case_when(
      str_detect(name, "vOTU") ~ "Phage vOTU",
      str_detect(name, "Rothia mucilaginosa") ~ "Target Host (R. mucilaginosa)",
      TRUE ~ "Other Host (Rothia spp.)"
    )
  )

# 3. 创建 igraph 对象
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

# 4. 使用 ggraph 绘制网络图
p_net <- ggraph(net, layout = 'fr') + # fr 布局(Fruchterman-Reingold)会自动聚类
  # 画边，添加连线样式
  geom_edge_link(
    aes(edge_alpha = 0.5),
    edge_width = 0.8,
    color = "grey50",
    arrow = arrow(length = unit(2, 'mm')),
    end_cap = circle(3, 'mm')
  ) +
  # 画节点，根据类型赋予颜色和大小
  geom_node_point(aes(color = Type, size = Type)) +
  # 添加节点标签（过滤掉vOTU名字，只显示宿主名字以免太拥挤；如果您想全显示可去掉ifelse）
  geom_node_text(
    aes(label = ifelse(Type != "Phage vOTU", name, name)),
    repel = TRUE,
    size = 3.5,
    fontface = "italic"
  ) +
  # 设置节点颜色
  scale_color_manual(
    values = c(
      "Target Host (R. mucilaginosa)" = "#E53935", # 靶标红色
      "Other Host (Rothia spp.)" = "#1E88E5", # 其他宿主蓝色
      "Phage vOTU" = "#BDBDBD" # 噬菌体灰色
    )
  ) +
  # 设置节点大小差异
  scale_size_manual(
    values = c(
      "Target Host (R. mucilaginosa)" = 8,
      "Other Host (Rothia spp.)" = 6,
      "Phage vOTU" = 3
    )
  ) +
  theme_void() + # 去除背景坐标轴
  theme(legend.position = "bottom") +
  labs(
    title = "Virus-Host Interaction Network",
    subtitle = "Highlighting specific vOTU targeting Rothia mucilaginosa",
    color = "Node Type",
    size = "Node Type"
  )

# 显示并保存网络图
print(p_net)
ggsave("Phage_Host_Network.svg", p_net, width = 8, height = 6)


# 美吉生物16SrRNA ------------------------------------------------------------

# 1. 构建数据框
rm_data <- data.frame(
  Relative_Abundance = c(
    0.05284,
    0.04542,
    0.07091,
    0.1066,
    0.1733,
    0.01311,
    0.005525,
    0.07501,
    0.008717,
    0.01765
  ),
  Group = factor(rep(c("OLP", "HC"), each = 5), levels = c("OLP", "HC"))
)

# 2. 设定分组颜色
olp_colors_two <- c(OLP = "#F4A58295", HC = "#4393C395")

# 3. 使用自定义函数进行可视化与统计检验
p_rm <- compare_plot(
  data = rm_data,
  value.var = "Relative_Abundance",
  group.by = "Group",
  fill.by = "Group",
  plot_type = "bar", # 对于 n=5 的微生物组数据，也可以考虑改成 "box" (箱线图) 看数据分布
  add_point = TRUE,
  add_stat = "wilcox.test", # 考虑到每组 n=5，如果担心不满足正态分布，可以替换为非参数检验 "wilcox.test"
  stat_label = c("p.format"),
  comparisons = list(c("OLP", "HC")),
  palette = olp_colors_two,
  xlab = "Group",
  ylab = "Relative Abundance of R. mucilaginosa",
  title = "16S rRNA Sequencing (n=5)",
  theme_use = theme_my_stat
) +
  fix_panel(width = 4, height = 5) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold") # 全局加粗
  )

# 4. 自动检查目录并保存 PDF
if (!dir.exists("16SrRNA")) {
  dir.create("16SrRNA")
}

ggsave("16SrRNA/Rm_Relative_Abundance.pdf", plot = p_rm, width = 6, height = 6)
