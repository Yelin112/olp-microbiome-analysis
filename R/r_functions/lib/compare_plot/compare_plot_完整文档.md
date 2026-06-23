# compare_plot() 函数说明文档

> 文件：`Scripts/compare_plot_optimized.R`  
> 更新日期：2025-04-05

---

## 概述

`compare_plot()` 是基于 ggplot2 的组间比较可视化封装函数，核心特性：

- **自适应策略**：根据每组样本量自动选择最合适的图形类型
- **灵活配置**：通过 `*_args` 参数精细控制各图层样式
- **科研配色**：内置 NPG、AAAS、NEJM 等期刊配色方案
- **统计检验**：集成 wilcox、t.test、anova 等检验及显著性标注

---

## 依赖包

```r
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
```

---

## 函数签名

```r
compare_plot(
  data, value.var, group.by,
  fill.by = NULL, split.by = NULL,

  # 自适应策略
  strategy = "auto",
  threshold_small = 5, threshold_medium = 10, threshold_large = 20,

  # 主图类型（strategy="none" 时手动指定）
  plot_type = "violin",

  # 散点控制
  add_point = FALSE,
  point_shape = 21, point_size = 2, point_alpha = 1,
  point_fill = "same", point_color = "black",
  point_stroke = 0.7, point_jitter = 0.2,

  # 图层细节参数
  box_args = list(), violin_args = list(), bar_args = list(),
  errorbar_args = list(), crossbar_args = list(),

  # 叠加元素
  add_box = FALSE, add_trend = FALSE, add_bg = FALSE,
  bg_color = "#0000000D",

  # 统计检验
  add_stat = "none", stat_label = "p.signif",
  comparisons = NULL, hide_ns = FALSE,
  step_increase = 0.12, y_expand = 0.15,

  # 外观
  palette = "NPG", theme_use = theme_classic,
  title = NULL, xlab = NULL, ylab = NULL, ...
)
```

---

## 参数详解

### 数据参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `data` | data.frame | — | 输入数据框 |
| `value.var` | 字符串 | — | Y 轴数值列名 |
| `group.by` | 字符串 | — | X 轴分组列名 |
| `fill.by` | 字符串 | `NULL` | 填充颜色列名。`NULL` 时与 `group.by` 相同；指定不同列时开启**组内簇状比较**模式 |
| `split.by` | 字符串 | `NULL` | 分面列名，使用 `facet_wrap` |

---

### 自适应策略参数

`strategy = "auto"` 时，函数自动根据 **最大组样本量（max n）** 选择图形类型：

| max n | 自动策略 | 图形效果 |
|-------|---------|---------|
| n < 5 | `pure_scatter` | 所有散点 + 均值线（crossbar）+ 标准误差线 |
| 5 ≤ n < 10 | `bar_points` | 均值柱状图 + 误差线 + 所有散点 |
| 10 ≤ n < 20 | `boxplot_points` | 箱线图 + 所有散点 |
| n ≥ 20 | `pure_boxplot` | 箱线图 + 均值点 |

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `strategy` | 字符串 | `"auto"` | `"auto"` 自动检测；`"none"` 禁用（使用 `plot_type`）；或直接指定四种策略之一 |
| `threshold_small` | 数值 | `5` | pure_scatter 与 bar_points 的分界（n < 此值用散点） |
| `threshold_medium` | 数值 | `10` | bar_points 与 boxplot_points 的分界 |
| `threshold_large` | 数值 | `20` | boxplot_points 与 pure_boxplot 的分界 |

> **注意**：`strategy` 会覆盖 `plot_type` 和 `add_point` 的设置。如需完全手动控制，设置 `strategy = "none"`。

---

### 主图类型（`strategy = "none"` 时生效）

| `plot_type` | 图形效果 |
|-------------|---------|
| `"violin"` | 小提琴图 |
| `"box"` | 箱线图 + 均值白点 |
| `"bar"` | 均值柱状图 + 标准误差线 |
| `"dot"` | 均值点 + 标准误差线 |
| `"scatter"` | 均值线（crossbar）+ 标准误差线（需配合 `add_point = TRUE`） |

---

### 散点参数

仅在 `add_point = TRUE`（或策略自动开启）时生效。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `point_shape` | `21` | 散点形状。21=填充圆，16=实心圆，22=填充方，23=填充菱形 |
| `point_size` | `2` | 点大小 |
| `point_alpha` | `1` | 透明度（1=不透明） |
| `point_fill` | `"same"` | 填充色。`"same"` 继承分组颜色；`"white"`、`"black"` 或任意色值 |
| `point_color` | `"black"` | 边框颜色（shape 21-25 的圈线颜色） |
| `point_stroke` | `0.7` | 边框粗细（shape 21-25 的圈线宽度） |
| `point_jitter` | `0.2` | 水平抖动幅度。`0` = 无抖动 |

---

### 图层细节参数（`*_args`）

通过传入命名列表覆盖对应图层的默认样式，**只需写想改的字段**，其余保持默认。

#### `box_args` — 箱线图样式

作用于 `plot_type = "box"` 的 `geom_boxplot()`。

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.7`（单组）/ `0.6`（组内） | 箱体宽度 |
| `alpha` | `0.6` | 填充透明度 |
| `color` | `"black"` | 边框颜色 |
| `linewidth` | `0.8` | 边框/须线宽度 |
| `outlier.shape` | `NA` | 离群点形状（`NA` = 隐藏，与散点层统一展示） |
| `fill` | — | 设为 `NA` 可使箱体透明（仅保留边框） |
| `notch` | `FALSE` | 是否显示置信缺口 |

#### `violin_args` — 小提琴图样式

作用于 `plot_type = "violin"` 的 `geom_violin()`。

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `alpha` | `0.6` | 填充透明度 |
| `color` | `"black"` | 轮廓颜色 |
| `linewidth` | `0.5` | 轮廓线宽 |
| `trim` | `FALSE` | 是否裁剪至数据范围 |
| `scale` | `"width"` | 宽度标准化方式（`"width"`, `"area"`, `"count"`） |

#### `bar_args` — 柱状图样式

作用于 `plot_type = "bar"` 的均值柱（`stat_summary(geom="bar")`）。

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.7` | 柱宽 |
| `alpha` | `0.6` | 填充透明度 |
| `color` | `"black"` | 边框颜色 |
| `linewidth` | `0.5` | 边框线宽 |
| `fill` | — | 设为 `NA` 可使柱体透明（空心柱） |

#### `errorbar_args` — 误差线样式

作用于 `bar`、`dot`、`scatter` 的误差线（`stat_summary(geom="errorbar")`）。

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.25` | 误差线帽宽 |
| `linewidth` | `0.6` | 误差线线宽 |

#### `crossbar_args` — 均值线样式

作用于 `scatter` 策略的均值 crossbar（`stat_summary(geom="crossbar")`）。

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.3` | crossbar 宽度 |
| `linewidth` | `1.2` | 线宽 |
| `fatten` | `1` | 中心线相对粗细 |

---

### 统计检验参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `add_stat` | `"none"` | 统计方法：`"mean"`, `"median"`, `"t.test"`, `"wilcox.test"`, `"anova"`, `"kruskal.test"` |
| `stat_label` | `"p.signif"` | 标注格式：`"p.signif"`（`*`/`**`/`***`）或 `"p.format"`（`p=0.023`） |
| `comparisons` | `NULL` | 指定比较对，如 `list(c("A","B"), c("A","C"))`。`NULL` 时自动生成全部两两比较 |
| `hide_ns` | `FALSE` | 是否隐藏 ns（不显著）标注 |
| `step_increase` | `0.12` | 显著性括号的垂直间距，多组时调大避免重叠 |
| `y_expand` | `0.15` | Y 轴顶部扩展比例，为标注留空间 |

---

### 外观参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `palette` | `"NPG"` | 配色方案。内置：`"NPG"`, `"AAAS"`, `"NEJM"`, `"Lancet"`, `"JCO"`, `"JAMA"`, `"D3"`；或传入颜色向量 |
| `theme_use` | `theme_classic` | ggplot2 主题函数 |
| `title` | `NULL` | 图标题 |
| `xlab` | `group.by` | X 轴标签（默认使用列名） |
| `ylab` | `value.var` | Y 轴标签（默认使用列名） |
| `add_bg` | `FALSE` | 斑马纹背景（交替浅色条带） |
| `bg_color` | `"#0000000D"` | 背景条带颜色（默认 5% 透明黑） |
| `add_box` | `FALSE` | 在小提琴图内叠加细箱线图 |
| `add_trend` | `FALSE` | 叠加均值连线（适合时间序列数据） |
| `...` | — | 传递给 `theme_use()` 的额外参数，如 `base_size = 14` |

---

## 使用示例

### 基础用法（自动策略）

```r
source("Scripts/compare_plot_optimized.R")

# 自动根据样本量选择图形
compare_plot(data, value.var = "abundance", group.by = "group")
```

### 禁用自动策略，手动指定图形类型

```r
compare_plot(data, value.var = "abundance", group.by = "group",
             strategy  = "none",
             plot_type = "violin",
             add_box   = TRUE,
             add_point = TRUE)
```

### 调整自动策略的阈值

```r
# n=12 时仍用柱状图+散点（threshold_medium 提高到 15）
compare_plot(data, value.var = "value", group.by = "group",
             threshold_medium = 15)
```

### 精细控制各图层样式

```r
# 箱线图：细箱体 + 加粗边框
compare_plot(data, value.var = "value", group.by = "group",
             strategy = "none", plot_type = "box",
             box_args = list(width = 0.5, linewidth = 1.2))

# 空心柱状图 + 黑色散点
compare_plot(data, value.var = "value", group.by = "group",
             strategy   = "none", plot_type = "bar",
             add_point  = TRUE,
             bar_args   = list(fill = NA, linewidth = 0.8),
             point_fill = "black", point_size = 1.5)

# 空心箱线图 + 有色散点
compare_plot(data, value.var = "value", group.by = "group",
             add_point = TRUE,
             box_args  = list(fill = NA))

# 误差线加粗、帽宽调窄
compare_plot(data, value.var = "value", group.by = "group",
             strategy      = "none", plot_type = "bar",
             errorbar_args = list(width = 0.15, linewidth = 1.0))
```

### 添加统计检验

```r
# 自动两两比较（Wilcoxon），显示星号
compare_plot(data, value.var = "value", group.by = "group",
             add_stat = "wilcox.test")

# 指定比较对，显示 p 值，隐藏 ns
compare_plot(data, value.var = "value", group.by = "group",
             add_stat    = "wilcox.test",
             stat_label  = "p.format",
             comparisons = list(c("Control", "Treatment")),
             hide_ns     = TRUE)
```

### 组内簇状比较

```r
# fill.by 不同于 group.by 时，自动切换为簇状模式
compare_plot(data,
             value.var = "score",
             group.by  = "timepoint",
             fill.by   = "treatment",
             strategy  = "none",
             plot_type = "box",
             add_stat  = "t.test")
```

### 自定义配色 + 主题

```r
compare_plot(data, value.var = "value", group.by = "group",
             palette   = c("#2196F3", "#FF5722", "#4CAF50"),
             theme_use = theme_bw,
             base_size = 14,
             title     = "My Plot",
             ylab      = "Expression (log2 CPM)")
```

### 分面

```r
compare_plot(data, value.var = "value", group.by = "group",
             split.by = "batch",
             add_bg   = TRUE)
```

### 在返回对象上继续定制

```r
p <- compare_plot(data, value.var = "value", group.by = "group")

p + theme(legend.position = "none") +
    scale_y_log10() +
    geom_hline(yintercept = 1, linetype = "dashed")
```

---

## 各策略的散点默认参数

策略自动激活时会覆盖以下散点参数，用户仍可手动传参覆盖：

| 策略 | `point_size` | `point_stroke` | `point_fill` | `add_point` |
|------|-------------|----------------|--------------|-------------|
| `pure_scatter` | 3 | 0.5 | `"same"` | `TRUE` |
| `bar_points` | 3 | 0.7（默认） | `"white"` | `TRUE` |
| `boxplot_points` | 3 | 0.7（默认） | `"same"` | `TRUE` |
| `pure_boxplot` | 2（默认） | 0.7（默认） | `"same"` | `FALSE` |

---

## 返回值

返回 **ggplot 对象**，可继续用 `+` 叠加图层或主题。

---

## 注意事项

1. **缺失值**：`value.var` 或 `group.by` 中的 `NA` 会自动移除并在 console 提示
2. **因子顺序**：分组变量会按数据中的出现顺序自动转为因子，如需自定义顺序，请在传入前手动设置 `factor(data$group, levels = c(...))`
3. **`fill = NA` 与散点颜色**：对 `box_args`/`bar_args` 等设置 `fill = NA` 可使主图透明，散点颜色正常保留（通过 `scale_fill_manual(na.value = "transparent")` 处理）
4. **strategy vs plot_type**：`strategy != "none"` 时会覆盖 `plot_type` 和 `add_point`；需要手动控制图形类型时请设置 `strategy = "none"`
