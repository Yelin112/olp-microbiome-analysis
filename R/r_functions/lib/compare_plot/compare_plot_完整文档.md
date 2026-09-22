# compare_plot() 函数说明文档

> 文件：`R/r_functions/lib/compare_plot/compare_plot_optimized.R`  
> 更新日期：2026-06-27

---

## 概述

`compare_plot()` 是基于 ggplot2 的组间比较可视化封装函数，核心特性：

- **自适应策略**：根据每组样本量自动选择最合适的图形类型
- **统计检验**：集成 wilcox、t.test、anova、kruskal 检验及显著性标注
- **灵活配置**：通过 `*_args` 参数精细控制各图层样式
- **统一配色**：通过 `R/utils/palette_system.R` 获取颜色——支持内置期刊配色、RColorBrewer、paletteer、自定义注册
- **统一主题**：通过 `R/utils/theme_system.R` 的 `theme_pub_base()` / `theme_nature()` 等一键切换
- **固定面板**：配合 `R/utils/panel_fix.R` 的 `fix_panel_size()` 统一多图 panel 尺寸

---

## 依赖与加载

```r
# 需要先加载工具模块（颜色 + 主题系统）
library(ggplot2)
library(ggpubr)
library(dplyr)

# 方式 1：绝对路径（推荐，最可靠）
PROJ <- "e:/打工人/Zeng/2023-3 OLP Microbiome/分析"
source(file.path(PROJ, "R/utils/helpers.R"))
source(file.path(PROJ, "R/utils/palette_system.R"))
source(file.path(PROJ, "R/r_functions/lib/compare_plot/compare_plot_optimized.R"))

# 方式 2：函数自带守卫加载
# compare_plot_optimized.R 内部会自动加载 utils，但如果调用方已经 source 过则跳过
source("R/r_functions/lib/compare_plot/compare_plot_optimized.R")

# 可选：主题系统
source(file.path(PROJ, "R/utils/theme_system.R"))

# 可选：固定面板尺寸
source(file.path(PROJ, "R/utils/panel_fix.R"))
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
  palette = NULL, theme_use = theme_classic,
  title = NULL, xlab = NULL, ylab = NULL, ...
)
```

---

## 配色系统（palette 参数）

### 来源优先级

`palette` 参数会传递给 `palette_system.R` 的 `get_colors()`，按以下优先级解析：

1. **自定义注册** — `set_palette("my_colors", c(...))` 注册的色板
2. **内置库** — NPG, AAAS, NEJM, Lancet, JCO, JAMA, D3, ColorBlind, RWB, Viridis
3. **RColorBrewer** — 所有 Brewer 色板名（如 "Set1", "Paired", "Spectral"…）
4. **paletteer** — 50+ 包的色板（如 `"ggsci::nrc_npg"`, `"viridis::viridis"`）
5. **回退** — 全局默认色板（默认 NPG）

### 使用方式

```r
# NULL — 使用全局默认色板（可在项目启动时 set_default_palette() 设定）
compare_plot(data, value.var = "value", group.by = "group")  # palette = NULL

# 内置色板名
compare_plot(data, ..., palette = "JAMA")
compare_plot(data, ..., palette = "NPG")

# 自定义颜色向量
compare_plot(data, ..., palette = c("#E41A1C", "#377EB8", "#4DAF4A"))

# RColorBrewer 色板
compare_plot(data, ..., palette = "Set1")

# 运行时注册的色板
set_palette("proj_colors", c("#E64B35", "#4DBBD5", "#00A087"))
compare_plot(data, ..., palette = "proj_colors")
```

---

## 主题系统（theme_use 参数）

```r
# pub_ 主题系统
source("R/utils/theme_system.R")

# 基石主题
compare_plot(..., theme_use = theme_pub_base)           # 白底黑框
compare_plot(..., theme_use = theme_pub_base, base_size = 14)  # 调字号

# 期刊预设
compare_plot(..., theme_use = theme_nature)             # Nature 风格
compare_plot(..., theme_use = theme_cell)               # Cell 风格
compare_plot(..., theme_use = theme_jama)               # JAMA 风格
compare_plot(..., theme_use = theme_science)            # Science 风格

# 场景主题
compare_plot(..., theme_use = theme_pub_stat)           # 统计图：虚线网格
compare_plot(..., theme_use = theme_pub_present)        # 演示：大字体

# 原生 ggplot2 主题也兼容
compare_plot(..., theme_use = theme_bw)
compare_plot(..., theme_use = theme_minimal)
```

---

## 固定面板尺寸

```r
source("R/utils/panel_fix.R")

# 生成图
p <- compare_plot(data, value.var = "value", group.by = "group")

# 固定 panel 为 12×8 cm，自动计算总尺寸后保存
g <- fix_panel_size(p, width = 12, height = 8)
dw <- convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
dh <- convertHeight(sum(g$heights), "cm", valueOnly = TRUE)
ggsave("output.pdf", g, width = dw, height = dh)
```

---

## 自适应策略参数

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
| `threshold_small` | 数值 | `5` | pure_scatter 与 bar_points 的分界 |
| `threshold_medium` | 数值 | `10` | bar_points 与 boxplot_points 的分界 |
| `threshold_large` | 数值 | `20` | boxplot_points 与 pure_boxplot 的分界 |

> **注意**：`strategy` 会覆盖 `plot_type` 和 `add_point`。如需完全手动控制，设置 `strategy = "none"`。

---

## 主图类型（`strategy = "none"` 时生效）

| `plot_type` | 图形效果 |
|-------------|---------|
| `"violin"` | 小提琴图 |
| `"box"` | 箱线图 + 均值白点 |
| `"bar"` | 均值柱状图 + 标准误差线 |
| `"dot"` | 均值点 + 标准误差线 |
| `"scatter"` | 均值线（crossbar）+ 标准误差线（需配合 `add_point = TRUE`） |

---

## 散点参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `point_shape` | `21` | 21=填充圆，16=实心圆，22=填充方，23=填充菱形 |
| `point_size` | `2` | 点大小 |
| `point_alpha` | `1` | 透明度 |
| `point_fill` | `"same"` | `"same"` 继承分组颜色；`"white"`、`"black"` 或任意色值 |
| `point_color` | `"black"` | 边框颜色 |
| `point_stroke` | `0.7` | 边框粗细 |
| `point_jitter` | `0.2` | 水平抖动幅度。`0` = 无抖动 |

---

## 图层细节参数（`*_args`）

通过传入命名列表覆盖对应图层的默认样式，**只需写想改的字段**。

### `box_args` — 箱线图

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.7` / 组内 `0.6` | 箱体宽度 |
| `alpha` | `0.6` | 填充透明度 |
| `color` | `"black"` | 边框颜色 |
| `linewidth` | `0.8` | 边框线宽 |
| `outlier.shape` | `NA` | 离群点形状（`NA` = 隐藏） |
| `fill` | — | 设为 `NA` 可使箱体透明 |

### `violin_args` — 小提琴图

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `alpha` | `0.6` | 填充透明度 |
| `color` | `"black"` | 轮廓颜色 |
| `linewidth` | `0.5` | 轮廓线宽 |
| `trim` | `FALSE` | 裁剪至数据范围 |
| `scale` | `"width"` | 宽度标准化 |

### `bar_args` — 柱状图

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.7` | 柱宽 |
| `alpha` | `0.6` | 填充透明度 |
| `color` | `"black"` | 边框颜色 |
| `linewidth` | `0.5` | 边框线宽 |
| `fill` | — | 设为 `NA` 使柱体透明 |

### `errorbar_args` — 误差线

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.25` | 帽宽 |
| `linewidth` | `0.6` | 线宽 |

### `crossbar_args` — 均值线（scatter 策略）

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `width` | `0.3` | 线宽 |
| `linewidth` | `1.2` | 线宽 |
| `fatten` | `1` | 中心线粗细 |

---

## 统计检验参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `add_stat` | `"none"` | `"t.test"`, `"wilcox.test"`, `"anova"`, `"kruskal.test"` |
| `stat_label` | `"p.signif"` | `"p.signif"`（`*`/`**`/`***`）或 `"p.format"`（`p=0.023`） |
| `comparisons` | `NULL` | 指定比较对，如 `list(c("A","B"), c("A","C"))`。NULL 时自动生成 |
| `hide_ns` | `FALSE` | 隐藏不显著标注 |
| `step_increase` | `0.12` | 显著性括号垂直间距 |
| `y_expand` | `0.15` | Y 轴顶部扩展比例 |

---

## 外观参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `palette` | `NULL` | NULL = 全局默认；色板名；颜色向量。通过 palette_system.R 解析 |
| `theme_use` | `theme_classic` | ggplot2 主题函数。推荐 `theme_pub_base()`、`theme_nature()` 等 |
| `title` | `NULL` | 图标题 |
| `xlab` | `group.by` | X 轴标签 |
| `ylab` | `value.var` | Y 轴标签 |
| `add_bg` | `FALSE` | 斑马纹背景 |
| `bg_color` | `"#0000000D"` | 背景条带颜色 |
| `add_box` | `FALSE` | 小提琴图内叠加细箱线图 |
| `add_trend` | `FALSE` | 叠加均值连线 |
| `...` | — | 传递给 `theme_use()` 的额外参数 |

---

## 使用示例

### 基础用法

```r
# 完全默认（自适应策略 + 全局默认配色）
compare_plot(iris, value.var = "Sepal.Length", group.by = "Species")

# 指定配色 + 主题
compare_plot(iris, value.var = "Sepal.Length", group.by = "Species",
             palette = "JAMA", theme_use = theme_nature())
```

### 添加统计检验

```r
# 整体 ANOVA
compare_plot(data, value.var = "thickness", group.by = "group",
             add_stat = "anova")

# 成对比较（Control vs 各处理组）
compare_plot(data, value.var = "thickness", group.by = "group",
             add_stat    = "wilcox.test",
             comparisons = list(c("Control","OXA"), c("Control","OXA_MN")),
             stat_label  = "p.signif", hide_ns = TRUE,
             y_expand    = 0.35)
```

### 组内簇状比较

```r
compare_plot(data, value.var = "score", group.by = "timepoint",
             fill.by = "treatment", strategy = "none",
             plot_type = "box", add_stat = "t.test")
```

### 精细控制图层 + 固定面板

```r
p <- compare_plot(data, value.var = "value", group.by = "group",
                  strategy = "none", plot_type = "box",
                  add_point = TRUE,
                  box_args  = list(width = 0.5, linewidth = 1.2),
                  palette   = "NPG", theme_use = theme_pub_base)

g <- fix_panel_size(p, width = 10, height = 7)
ggsave("output.pdf", g,
  width  = convertWidth(sum(g$widths), "cm", valueOnly = TRUE),
  height = convertHeight(sum(g$heights), "cm", valueOnly = TRUE))
```

### 分面

```r
compare_plot(data, value.var = "value", group.by = "group",
             split.by = "batch", add_bg = TRUE)
```

---

## 返回值

返回 **ggplot 对象**，可继续用 `+` 叠加图层或主题。

对比 `compare_plot_adaptive()`（dev 版）返回 `list(plot, strategy, n_per_group, metadata)`，本函数仅返回 ggplot 对象。如需元数据，使用 dev 版。

---

## 相关文件

| 文件 | 说明 |
|------|------|
| `R/utils/palette_system.R` | 颜色唯一真相来源 |
| `R/utils/theme_system.R` | 主题工厂 + 期刊预设 |
| `R/utils/panel_fix.R` | 固定面板尺寸 |
| `R/utils/README.md` | 主题与配色系统完整手册 |
| `R/utils/export_pptx.R` | 图表导出为可编辑 PPTX |
| `R/PLOTTING_CONVENTIONS.md` | 分析图表产出规范：导出格式、图内文字用英文、组间比较统计设计（vs Control、Kruskal-Wallis 放 subtitle、百分比轴上限等） |
| `R/r_functions/dev/compare_plot_adaptive/` | 管线式 v3 实验版 |

---

## 注意事项

1. **缺失值**：`value.var` 或 `group.by` 中的 `NA` 自动移除
2. **因子顺序**：按数据出现顺序自动转因子，自定义顺序请在传入前设置
3. **strategy vs plot_type**：`strategy != "none"` 时覆盖 `plot_type` 和 `add_point`
4. **palette = NULL**：使用 palette_system.R 的全局默认色板，由 `set_default_palette()` 控制
5. **theme_use**：推荐使用 `theme_pub_base()` 替代 `theme_classic()`，获得统一白底黑框风格
