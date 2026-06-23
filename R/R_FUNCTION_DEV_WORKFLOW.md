# R函数优化项目 - 基于 r-function-dev 工作流

## 🎯 项目目标

优化 `compare_plot_optimized` 函数，实现基于每组样本量的自动可视化策略切换。

---

## 📂 项目结构（遵循 r-function-dev 规范）

### 建议在现有Scripts目录下创建

```
E:\打工人\Zeng\2023-3 OLP Microbiome\分析\Scripts\
│
├── [现有的其他文件和文件夹...]
│
└── r_functions/                    # 新建：R函数管理根目录
    ├── dev/                        # 开发区：正在开发的函数
    │   └── compare_plot_adaptive/  # 本次优化项目
    │       ├── function.R          # 自适应绘图主函数
    │       ├── test.R              # 测试脚本
    │       ├── example.R           # 使用示例
    │       ├── notes.md            # 开发笔记
    │       └── strategies/         # 策略函数（子模块）
    │           ├── detect_strategy.R
    │           ├── pure_scatter.R
    │           ├── bar_points.R
    │           ├── boxplot_points.R
    │           └── pure_boxplot.R
    │
    ├── lib/                        # 函数库：成熟函数
    │   ├── visualization.R         # 归档后放这里
    │   ├── data_processing.R
    │   └── statistics.R
    │
    └── archive/                    # 归档：旧版本备份
        └── compare_plot_optimized_v1_20250405.R
```

---

## 🔄 开发流程（遵循 r-function-dev）

### Phase 1: 项目初始化

#### 1.1 创建开发目录
```r
# 在Scripts目录下创建结构
dir.create("r_functions/dev/compare_plot_adaptive/strategies", recursive = TRUE)
dir.create("r_functions/lib", recursive = TRUE)
dir.create("r_functions/archive", recursive = TRUE)
```

#### 1.2 备份现有函数
```r
# 将原函数复制到archive
file.copy(
  "compare_plot_optimized.R",
  "r_functions/archive/compare_plot_optimized_v1_20250405.R"
)
```

#### 1.3 创建开发笔记
在 `r_functions/dev/compare_plot_adaptive/notes.md` 中记录：
- 开发目标
- 参考资料（对话中的参考图）
- 当前进度
- 已知问题

---

### Phase 2: 模块化开发（在 dev/ 中）

#### 2.1 策略检测模块
**文件**: `strategies/detect_strategy.R`

```r
#' 检测最佳可视化策略
#'
#' 根据每组样本量自动选择最适合的可视化策略
#'
#' @param data 数据框
#' @param group_col 分组列名（字符串）
#' @param threshold_small 小样本阈值（默认5）
#' @param threshold_medium 中等样本阈值（默认10）
#' @param threshold_large 大样本阈值（默认20）
#'
#' @return 包含策略名称和样本量统计的列表
#'
#' @importFrom dplyr group_by summarise pull
#' @importFrom rlang sym
#'
#' @examples
#' data <- data.frame(group = rep(c("A", "B"), each = 5), value = rnorm(10))
#' detect_strategy(data, "group")
#'
#' @export
detect_strategy <- function(data, 
                             group_col,
                             threshold_small = 5,
                             threshold_medium = 10,
                             threshold_large = 20) {
  
  # 参数验证
  if (!group_col %in% names(data)) {
    stop("Column '", group_col, "' not found in data")
  }
  
  # 计算每组样本量
  n_per_group <- data %>%
    group_by(!!sym(group_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    pull(n)
  
  # 统计量
  max_n <- max(n_per_group)
  median_n <- median(n_per_group)
  min_n <- min(n_per_group)
  
  # 策略决策树
  if (max_n < threshold_small) {
    strategy <- "pure_scatter"
    message("Strategy: Pure scatter (max n=", max_n, " < ", threshold_small, ")")
  } else if (max_n < threshold_medium) {
    strategy <- "bar_points"
    message("Strategy: Bar + points (max n=", max_n, " < ", threshold_medium, ")")
  } else if (max_n < threshold_large) {
    strategy <- "boxplot_points"
    message("Strategy: Boxplot + points (max n=", max_n, " < ", threshold_large, ")")
  } else {
    strategy <- "pure_boxplot"
    message("Strategy: Pure boxplot (max n=", max_n, " >= ", threshold_large, ")")
  }
  
  # 返回结果
  return(list(
    strategy = strategy,
    n_per_group = n_per_group,
    max_n = max_n,
    median_n = median_n,
    min_n = min_n
  ))
}
```

#### 2.2 策略实现模块
**文件**: `strategies/pure_scatter.R`, `bar_points.R`, `boxplot_points.R`, `pure_boxplot.R`

每个策略一个独立文件，包含完整的Roxygen2注释。

---

#### 2.3 主函数
**文件**: `function.R`

```r
#' 自适应组间比较可视化
#'
#' 根据每组样本量自动选择最佳可视化策略，生成发表级图形
#'
#' @param data 数据框，包含分组变量和数值变量
#' @param group_col 分组列名（字符串）
#' @param value_col 数值列名（字符串）
#' @param strategy 策略选择。"auto"（默认）自动检测，或手动指定：
#'   "pure_scatter", "bar_points", "boxplot_points", "pure_boxplot"
#' @param threshold_small 样本量阈值：n < 此值使用纯散点（默认5）
#' @param threshold_medium 样本量阈值：n < 此值使用柱状图+散点（默认10）
#' @param threshold_large 样本量阈值：n < 此值使用箱线图+散点（默认20）
#' @param color_palette 配色方案。"default", "colorblind", "grayscale" 或自定义向量
#' @param show_points 是否显示散点（在纯箱线图策略中）
#' @param verbose 是否打印策略信息
#' @param ... 传递给ggplot2的其他参数
#'
#' @return 包含以下元素的列表：
#'   \item{plot}{ggplot对象}
#'   \item{strategy}{使用的策略名称}
#'   \item{n_per_group}{各组样本量}
#'   \item{metadata}{其他元数据}
#'
#' @importFrom ggplot2 ggplot aes geom_jitter geom_bar geom_boxplot geom_violin
#' @importFrom ggplot2 stat_summary scale_fill_manual scale_color_manual theme_classic
#' @importFrom dplyr group_by summarise pull
#' @importFrom rlang sym
#'
#' @examples
#' # 示例1：小样本（自动使用纯散点）
#' data_small <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   value = rnorm(9, mean = rep(c(30, 35, 28), each = 3), sd = 2)
#' )
#' result <- compare_plot_adaptive(data_small, "group", "value")
#' print(result$plot)
#'
#' # 示例2：大样本（自动使用箱线图）
#' data_large <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 25),
#'   value = rnorm(75, mean = rep(c(30, 35, 28), each = 25), sd = 5)
#' )
#' result <- compare_plot_adaptive(data_large, "group", "value")
#' print(result$plot)
#'
#' # 示例3：手动指定策略
#' result <- compare_plot_adaptive(data_small, "group", "value", strategy = "bar_points")
#'
#' @export
compare_plot_adaptive <- function(data,
                                   group_col,
                                   value_col,
                                   strategy = "auto",
                                   threshold_small = 5,
                                   threshold_medium = 10,
                                   threshold_large = 20,
                                   color_palette = "default",
                                   show_points = TRUE,
                                   verbose = TRUE,
                                   ...) {
  
  # 1. 参数验证
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  if (!group_col %in% names(data)) {
    stop("Column '", group_col, "' not found in data")
  }
  if (!value_col %in% names(data)) {
    stop("Column '", value_col, "' not found in data")
  }
  
  # 2. 加载策略函数
  source("strategies/detect_strategy.R")
  source("strategies/pure_scatter.R")
  source("strategies/bar_points.R")
  source("strategies/boxplot_points.R")
  source("strategies/pure_boxplot.R")
  
  # 3. 策略检测
  if (strategy == "auto") {
    detection <- detect_strategy(data, group_col, 
                                  threshold_small, threshold_medium, threshold_large)
    strategy <- detection$strategy
    n_info <- detection
  } else {
    # 手动指定策略
    if (!strategy %in% c("pure_scatter", "bar_points", "boxplot_points", "pure_boxplot")) {
      stop("Invalid strategy. Choose from: auto, pure_scatter, bar_points, boxplot_points, pure_boxplot")
    }
    if (verbose) message("Manual strategy override: ", strategy)
    n_info <- detect_strategy(data, group_col)  # 仍需获取样本量信息
  }
  
  # 4. 配色方案
  colors <- get_color_palette(color_palette, n_groups = length(unique(data[[group_col]])))
  
  # 5. 根据策略绘图
  p <- switch(strategy,
    pure_scatter = plot_pure_scatter(data, group_col, value_col, colors),
    bar_points = plot_bar_points(data, group_col, value_col, colors),
    boxplot_points = plot_boxplot_points(data, group_col, value_col, colors, use_violin = FALSE),
    pure_boxplot = plot_pure_boxplot(data, group_col, value_col, colors, show_points = show_points)
  )
  
  # 6. 返回结果
  return(list(
    plot = p,
    strategy = strategy,
    n_per_group = n_info$n_per_group,
    metadata = list(
      max_n = n_info$max_n,
      median_n = n_info$median_n,
      min_n = n_info$min_n,
      thresholds = c(small = threshold_small, 
                     medium = threshold_medium, 
                     large = threshold_large)
    )
  ))
}
```

---

### Phase 3: 测试文件（test.R）

```r
# 测试文件：compare_plot_adaptive
# 创建日期：2025-04-05

# 加载函数
source("function.R")

# 加载必要的包
library(dplyr)
library(ggplot2)

# === 测试1: 极小样本（n < 5）===
cat("测试1: 极小样本 (n=3)\n")
data_small <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 3),
  value = c(
    rnorm(3, 35, 3),
    rnorm(3, 38, 4),
    rnorm(3, 25, 2),
    rnorm(3, 30, 3)
  )
)

result1 <- compare_plot_adaptive(data_small, "group", "value")
stopifnot(result1$strategy == "pure_scatter")
cat("✓ 策略正确: pure_scatter\n")
cat("✓ 样本量:", result1$n_per_group, "\n\n")

# 保存测试图
ggsave("test_output/test1_pure_scatter.png", result1$plot, width = 6, height = 4)

# === 测试2: 小样本（5 ≤ n < 10）===
cat("测试2: 小样本 (n=7)\n")
data_medium_small <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 7),
  value = c(
    rnorm(7, 35, 5),
    rnorm(7, 38, 6),
    rnorm(7, 25, 4),
    rnorm(7, 30, 5)
  )
)

result2 <- compare_plot_adaptive(data_medium_small, "group", "value")
stopifnot(result2$strategy == "bar_points")
cat("✓ 策略正确: bar_points\n")
cat("✓ 样本量:", result2$n_per_group, "\n\n")

ggsave("test_output/test2_bar_points.png", result2$plot, width = 6, height = 4)

# === 测试3: 中等样本（10 ≤ n < 20）===
cat("测试3: 中等样本 (n=15)\n")
data_medium_large <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 15),
  value = c(
    rnorm(15, 35, 5),
    rnorm(15, 38, 6),
    rnorm(15, 25, 4),
    rnorm(15, 30, 5)
  )
)

result3 <- compare_plot_adaptive(data_medium_large, "group", "value")
stopifnot(result3$strategy == "boxplot_points")
cat("✓ 策略正确: boxplot_points\n")
cat("✓ 样本量:", result3$n_per_group, "\n\n")

ggsave("test_output/test3_boxplot_points.png", result3$plot, width = 6, height = 4)

# === 测试4: 大样本（n ≥ 20）===
cat("测试4: 大样本 (n=30)\n")
data_large <- data.frame(
  group = rep(c("A", "B", "C", "D"), each = 30),
  value = c(
    rnorm(30, 35, 5),
    rnorm(30, 38, 6),
    rnorm(30, 25, 4),
    rnorm(30, 30, 5)
  )
)

result4 <- compare_plot_adaptive(data_large, "group", "value")
stopifnot(result4$strategy == "pure_boxplot")
cat("✓ 策略正确: pure_boxplot\n")
cat("✓ 样本量:", result4$n_per_group, "\n\n")

ggsave("test_output/test4_pure_boxplot.png", result4$plot, width = 6, height = 4)

# === 测试5: 混合样本量 ===
cat("测试5: 混合样本量 (3, 8, 15, 30)\n")
data_mixed <- data.frame(
  group = c(rep("A", 3), rep("B", 8), rep("C", 15), rep("D", 30)),
  value = c(
    rnorm(3, 35, 3),
    rnorm(8, 38, 4),
    rnorm(15, 25, 3),
    rnorm(30, 30, 5)
  )
)

result5 <- compare_plot_adaptive(data_mixed, "group", "value")
stopifnot(result5$strategy == "pure_boxplot")  # 按max_n=30决定
cat("✓ 策略正确: pure_boxplot (按max_n决定)\n")
cat("✓ 样本量:", result5$n_per_group, "\n\n")

ggsave("test_output/test5_mixed.png", result5$plot, width = 6, height = 4)

# === 测试6: 边界情况（n恰好等于阈值）===
cat("测试6: 边界情况 (n=5, 10, 20)\n")

# n=5 应该用 bar_points
data_boundary1 <- data.frame(group = rep("A", 5), value = rnorm(5))
result6a <- compare_plot_adaptive(data_boundary1, "group", "value")
stopifnot(result6a$strategy == "bar_points")

# n=10 应该用 boxplot_points
data_boundary2 <- data.frame(group = rep("A", 10), value = rnorm(10))
result6b <- compare_plot_adaptive(data_boundary2, "group", "value")
stopifnot(result6b$strategy == "boxplot_points")

# n=20 应该用 pure_boxplot
data_boundary3 <- data.frame(group = rep("A", 20), value = rnorm(20))
result6c <- compare_plot_adaptive(data_boundary3, "group", "value")
stopifnot(result6c$strategy == "pure_boxplot")

cat("✓ 所有边界情况正确\n\n")

# === 测试7: 手动覆盖策略 ===
cat("测试7: 手动指定策略\n")
result7 <- compare_plot_adaptive(data_large, "group", "value", strategy = "bar_points")
stopifnot(result7$strategy == "bar_points")
cat("✓ 手动覆盖成功\n\n")

# === 测试8: 缺失值处理 ===
cat("测试8: 缺失值处理\n")
data_na <- data.frame(
  group = rep(c("A", "B"), each = 5),
  value = c(rnorm(4), NA, rnorm(4), NA)
)

# 应该在函数内部处理NA（删除或提示）
tryCatch({
  result8 <- compare_plot_adaptive(data_na, "group", "value")
  cat("✓ 缺失值处理正常\n\n")
}, error = function(e) {
  cat("✓ 缺失值正确报错:", e$message, "\n\n")
})

cat("\n=== 所有测试通过！ ===\n")
```

---

### Phase 4: 使用示例（example.R）

```r
# 使用示例：compare_plot_adaptive

# 加载函数
source("function.R")
library(dplyr)
library(ggplot2)

# === 示例1: 基本使用（自动策略）===
# 场景：小型实验，每组3个生物学重复
cat("示例1: 小样本实验数据\n")

data_exp1 <- data.frame(
  group = rep(c("Control", "Treatment_A", "Treatment_B"), each = 3),
  expression = c(
    c(5.2, 5.5, 5.1),     # Control
    c(7.3, 7.8, 7.5),     # Treatment_A
    c(4.1, 4.3, 4.2)      # Treatment_B
  )
)

result1 <- compare_plot_adaptive(
  data = data_exp1,
  group_col = "group",
  value_col = "expression"
)

print(result1$plot)
# 自动使用 pure_scatter 策略

# === 示例2: 大型实验数据 ===
# 场景：流式细胞术数据，每组30个样本
cat("\n示例2: 大样本流式数据\n")

data_flow <- data.frame(
  cell_type = rep(c("CD4+", "CD8+", "NK", "B"), each = 30),
  percentage = c(
    rnorm(30, 40, 8),
    rnorm(30, 25, 6),
    rnorm(30, 15, 4),
    rnorm(30, 20, 5)
  )
)

result2 <- compare_plot_adaptive(
  data = data_flow,
  group_col = "cell_type",
  value_col = "percentage"
)

print(result2$plot +
  labs(title = "Cell Type Distribution",
       y = "Percentage (%)",
       x = "Cell Type"))
# 自动使用 pure_boxplot 策略

# === 示例3: 自定义配色 ===
cat("\n示例3: 使用色盲友好配色\n")

result3 <- compare_plot_adaptive(
  data = data_exp1,
  group_col = "group",
  value_col = "expression",
  color_palette = "colorblind"
)

print(result3$plot)

# === 示例4: 手动指定策略 ===
cat("\n示例4: 强制使用柱状图\n")

result4 <- compare_plot_adaptive(
  data = data_flow,  # 虽然是大样本
  group_col = "cell_type",
  value_col = "percentage",
  strategy = "bar_points",  # 但强制用柱状图
  verbose = TRUE
)

print(result4$plot)

# === 示例5: 自定义阈值 ===
cat("\n示例5: 调整策略切换点\n")

result5 <- compare_plot_adaptive(
  data = data_exp1,
  group_col = "group",
  value_col = "expression",
  threshold_small = 6,    # 改为n<6用纯散点
  threshold_medium = 12,  # n<12用柱状图
  threshold_large = 25    # n<25用箱线图+散点
)

print(result5$plot)

# === 示例6: 进一步定制ggplot对象 ===
cat("\n示例6: 在结果基础上继续定制\n")

p <- result1$plot +
  labs(
    title = "Gene Expression Comparison",
    subtitle = paste("Strategy:", result1$strategy),
    y = "Expression Level (log2 CPM)",
    x = NULL
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50")
  )

print(p)
ggsave("my_final_plot.png", p, width = 7, height = 5, dpi = 300)
```

---

### Phase 5: 开发笔记（notes.md）

```markdown
# 开发笔记：compare_plot_adaptive

## 项目信息
- **开始日期**: 2025-04-05
- **原函数**: compare_plot_optimized.R
- **目标**: 添加基于样本量的自适应可视化策略

---

## 版本历史
- **v1** (归档于 archive/): 原始版本，手动调整参数
- **v2** (当前开发): 自适应版本，自动策略切换

---

## 开发目标

### 核心需求
根据每组样本量自动选择可视化策略：
- n < 5: 纯散点图
- 5 ≤ n < 10: 柱状图 + 散点
- 10 ≤ n < 20: 箱线图/小提琴图 + 散点  
- n ≥ 20: 仅箱线图/小提琴图

### 设计原则
1. 按 max_n（最大组样本量）决定策略
2. 用户可自定义阈值
3. 支持手动覆盖
4. 向后兼容（如果可能）

---

## 参考资料

### 对话记录
- Claude chat conversation (2025-04-05)
- 需求文档: CLAUDE_CODE_REQUIREMENTS.md
- 快速指南: QUICK_START_GUIDE.md

### 参考图片
1. **Dnase2a图**: 纯散点 + 误差线风格
   - 大散点 (size=3)
   - 黑色边框 (stroke=0.5)
   - 粗误差线 (size=1.2)

2. **彩色4组柱状图**: 柱状图 + 散点风格
   - 柱子宽度 0.7
   - 黑色散点 (fill="black")
   - 配色: #4DBBD5, #E64B35, #00A087, #F39B7F

3. **RBC箱线图**: 箱线图 + 散点风格
   - 箱线图 alpha=0.7
   - 小散点 (size=1.5, alpha=0.6)
   - 不重复显示离群点

---

## 当前进度

### 已完成 ✓
- [x] Phase 1: 项目结构建立
- [x] Phase 2: 策略检测模块
- [ ] Phase 3: 策略1实现 (pure_scatter)
- [ ] Phase 4: 策略2实现 (bar_points)
- [ ] Phase 5: 策略3实现 (boxplot_points)
- [ ] Phase 6: 策略4实现 (pure_boxplot)
- [ ] Phase 7: 主函数整合
- [ ] Phase 8: 完整测试
- [ ] Phase 9: 归档到 lib/

---

## 技术细节

### 依赖包
```r
library(dplyr)      # 数据处理
library(ggplot2)    # 绘图
library(rlang)      # NSE处理
```

### 关键函数
1. `detect_strategy()`: 样本量检测和策略选择
2. `plot_pure_scatter()`: 纯散点绘图
3. `plot_bar_points()`: 柱状图+散点绘图
4. `plot_boxplot_points()`: 箱线图+散点绘图
5. `plot_pure_boxplot()`: 纯箱线图绘图
6. `compare_plot_adaptive()`: 主函数

---

## 已知问题

### 待解决
- [ ] 缺失值处理策略（删除 vs 提示）
- [ ] 极端值对箱线图的影响
- [ ] 统计检验集成（ggpubr vs ggsignif）
- [ ] 是否需要支持分面（facet）

### 已解决 ✓
- [x] NSE列名处理 (使用 !!sym())
- [x] 策略决策逻辑

---

## 性能考虑
- 小数据集 (<1000行): 性能不是瓶颈
- 大数据集 (>10000行): 散点可能需要采样

---

## 下一步计划
1. 实现四种策略的绘图函数
2. 编写完整测试用例
3. 生成示例图片
4. 决定归档策略（替换原函数 vs 新增函数）
5. 更新文档

---

## 笔记
- 用户希望在VSCode中使用Claude Code插件进行开发
- 遵循 r-function-dev workflow
- 强调不要一次性重写，而是逐步测试
```

---

## 📋 Claude Code 开发清单

### 提示词序列（Phase by Phase）

#### Phase 1: 初始化项目
```
请帮我在 E:\打工人\Zeng\2023-3 OLP Microbiome\分析\Scripts\ 
目录下建立以下结构：

r_functions/
├── dev/compare_plot_adaptive/
│   ├── strategies/
│   └── test_output/
├── lib/
└── archive/

然后：
1. 复制 compare_plot_optimized.R 到 archive/ 并重命名为 compare_plot_optimized_v1_20250405.R
2. 创建 dev/compare_plot_adaptive/notes.md（使用我提供的模板内容）

确认完成后告诉我。
```

#### Phase 2: 实现策略检测
```
在 dev/compare_plot_adaptive/strategies/ 目录下创建 detect_strategy.R

要求：
1. 实现detect_strategy()函数（参考我提供的模板）
2. 包含完整的Roxygen2注释
3. 参数验证
4. 返回策略名称和样本量统计

完成后用以下数据测试：
- 测试1: n=3 → 应返回"pure_scatter"
- 测试2: n=7 → 应返回"bar_points"
- 测试3: n=15 → 应返回"boxplot_points"
- 测试4: n=30 → 应返回"pure_boxplot"

给我看测试结果。
```

#### Phase 3-6: 实现四种策略
```
现在依次实现四种绘图策略。

先从策略1开始：在 strategies/ 下创建 pure_scatter.R

实现 plot_pure_scatter() 函数：
- 大散点 (size=3, shape=21, black边框)
- 中心线 (均值, crossbar)
- 误差线 (标准误, size=1.2)
- theme_classic主题

参考需求文档中的"策略1"详细参数。

完成后用测试数据验证，保存图片到test_output/
```

（类似地依次完成策略2、3、4）

#### Phase 7: 整合主函数
```
所有策略都测试通过后，创建 dev/compare_plot_adaptive/function.R

实现主函数 compare_plot_adaptive()：
1. 参数验证
2. 加载所有策略函数（source()）
3. 调用detect_strategy()
4. 根据策略调用相应绘图函数
5. 返回list(plot, strategy, metadata)

参考我提供的function.R模板。
```

#### Phase 8: 完整测试
```
创建 dev/compare_plot_adaptive/test.R

包含8个测试用例：
1. 极小样本 (n=3)
2. 小样本 (n=7)
3. 中等样本 (n=15)
4. 大样本 (n=30)
5. 混合样本量
6. 边界情况 (n=5,10,20)
7. 手动覆盖策略
8. 缺失值处理

运行测试并保存所有测试图到test_output/

告诉我哪些测试通过，哪些失败。
```

#### Phase 9: 创建示例
```
创建 dev/compare_plot_adaptive/example.R

包含6个实际使用示例：
1. 基本使用（自动策略）
2. 大型实验数据
3. 自定义配色
4. 手动指定策略
5. 自定义阈值
6. 进一步定制ggplot对象

参考我提供的example.R模板。
```

#### Phase 10: 归档决策
```
现在所有功能都已完成并测试通过。

请帮我分析：
1. 是否应该替换原来的 compare_plot_optimized.R？
2. 还是作为新函数 compare_plot_adaptive.R 共存？
3. 如何保持向后兼容？

给出建议和实施方案。
```

---

## ✅ 验收标准

### 功能性
- [ ] 各样本量梯度自动识别策略
- [ ] 手动指定策略能覆盖自动检测
- [ ] 所有测试用例通过
- [ ] 图形符合参考风格

### 代码质量
- [ ] 所有函数有Roxygen2注释
- [ ] 参数验证完整
- [ ] 错误处理友好
- [ ] 代码风格一致

### 文档质量
- [ ] notes.md记录开发过程
- [ ] example.R有实际使用示例
- [ ] test.R覆盖所有场景

---

## 🎯 关键优势

遵循 **r-function-dev** 规范的好处：
1. ✅ **结构清晰**: dev/lib/archive分离
2. ✅ **易于管理**: 每个函数独立文件夹
3. ✅ **版本控制**: archive自动备份
4. ✅ **测试完善**: test.R强制测试
5. ✅ **文档完整**: notes.md记录全过程

---

## 📌 下一步行动

1. **在VSCode中打开项目目录**
2. **启动Claude Code**
3. **按Phase顺序执行**（使用上面的提示词）
4. **每完成一个Phase就测试验证**
5. **在notes.md中更新进度**

**准备好开始了吗？** 🚀
