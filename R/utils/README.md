# 主题与配色系统 — 使用手册

> **适用场景**: 本项目的所有 ggplot2 制图，包括日常 ad-hoc 分析和上层绘图函数开发。  
> **最后更新**: 2026-06-27
>
> 完整分析图表产出规范（导出格式、图内文字用英文、组间比较统计设计）见 [../PLOTTING_CONVENTIONS.md](../PLOTTING_CONVENTIONS.md)。

---

## 1. 系统概览

```
R/utils/
├── helpers.R              # 零依赖基础工具 (%||%, %ni%)
├── palette_system.R       # 颜色唯一真相来源 (调色板库 + 解析引擎 + scale 函数)
├── theme_system.R         # 主题工厂 (基石 → 场景 → 期刊预设)
├── panel_fix.R            # 固定面板尺寸 (独立工具)
├── export_pptx.R          # ggplot2 → 可编辑 PPTX 导出 (officer + rvg::dml())
├── FUNCTION_TEMPLATE.R    # 新绘图函数开发模板
├── archive/               # 已废弃的旧主题/配色系统备份
└── README.md              # 本文件
```

**依赖关系**: `helpers.R` → `palette_system.R` ← `theme_system.R`；`export_pptx.R` 依赖 `panel_fix.R`。

> 领域特定的分析函数库（微生物丰度过滤/标准化/差异检验、网络分析等）不属于"通用工具"，
> 已迁移至 `R/r_functions/lib/`（如 `abundance_processing/`、`network_analysis/`），
> 遵循 `R_FUNCTION_DEV_WORKFLOW.md` 的 dev → lib → archive 流程管理。
> 本目录只保留与绘图/主题/配色/导出相关、跨项目通用的工具。

---

## 2. 日常 ggplot2 使用 (最重要)

两个文件就够了。任何 `.R` 脚本开头加两行：

```r
source("R/utils/helpers.R")
source("R/utils/palette_system.R")   # 颜色
source("R/utils/theme_system.R")     # 主题 (可选)
```

### 2.1 快速出图 (最常用)

```r
# 箱线图 — Nature 风格 + NPG 配色
ggplot(iris, aes(Species, Sepal.Length, fill = Species)) +
  geom_boxplot() +
  scale_fill_pub_d("NPG") +
  theme_pub_base()

# 柱状图 — JAMA 风格
ggplot(mtcars, aes(factor(cyl), mpg, fill = factor(cyl))) +
  stat_summary(fun = mean, geom = "bar") +
  scale_fill_pub_d("JAMA") +
  theme_jama()

# 散点图 — Cell 风格 (color 而非 fill)
ggplot(iris, aes(Sepal.Length, Petal.Length, color = Species)) +
  geom_point(size = 2) +
  scale_color_pub_d("NPG") +
  theme_cell()

# 热图 — 连续色标
ggplot(faithfuld, aes(waiting, eruptions, fill = density)) +
  geom_tile() +
  scale_fill_pub_c("RWB") +
  theme_pub_base(grid = FALSE)
```

### 2.2 极简启动 (三个字母记住全部 API)

| 操作 | 离散 fill | 离散 color | 连续 fill | 连续 color |
|---|---|---|---|---|
| 配色 | `scale_fill_pub_d("NPG")` | `scale_color_pub_d("NPG")` | `scale_fill_pub_c("RWB")` | `scale_color_pub_c("RWB")` |

| 场景 | 主题函数 | 特点 |
|---|---|---|
| 基础 | `theme_pub_base()` | 白底黑字，可通过参数控制边框/网格 |
| 统计图 | `theme_pub_stat()` | 虚线网格 + 右侧图例 |
| 演示 | `theme_pub_present()` | 大字体 + 无网格 |
| 海报 | `theme_pub_poster()` | 超大字体 |
| Nature | `theme_nature()` | 黑细边框 + 右侧图例 |
| Cell | `theme_cell()` | 极简 + 浅灰网格 |
| JAMA | `theme_jama()` | 无边框 + 粗轴线 |
| Science | `theme_science()` | 紧凑 + 薄边框 |

### 2.3 进阶：直接用颜色向量

```r
# 不用色板名，直接传颜色
ggplot(df, aes(x, y, fill = group)) +
  geom_boxplot() +
  scale_fill_pub_d(c("#E41A1C", "#377EB8", "#4DAF4A"))
```

### 2.4 主题参数微调

```r
# theme_pub_base 的三个开关 + ... 透传
theme_pub_base(base_size = 14, border = FALSE, grid = TRUE)
theme_pub_base(legend.position = "none")   # ... 传给 theme()
```

---

## 3. 配色管理

### 3.1 查看所有可用色板

```r
list_palettes("all")       # 内置 + 自定义
list_palettes("builtin")   # 只看内置
```

### 3.2 预览色板效果

```r
preview_palette("NPG")           # 预览某一套
preview_palette("JCO")           # 换一套看
```

### 3.3 注册自定义配色 (运行时)

```r
# 单个注册
set_palette("my_project", c("#FF6B6B", "#4ECDC4", "#45B7D1"))

# 批量注册
set_palette(list(
  "microbiome"  = c("#E64B35", "#4DBBD5", "#00A087", "#F39B7F"),
  "tcell"       = c("#BC3C29", "#0072B5", "#E18727", "#20854E"),
  "cytokine"    = c("#440154", "#3B528B", "#21908C", "#5DC863")
))

# 设为自己最常用的默认
set_default_palette("microbiome")

# 之后不传 palette 参数就用这个默认值
scale_fill_pub_d()   # 自动用 "microbiome"
```

### 3.4 永久添加新配色

在 `palette_system.R` 的 `.builtin_palettes` 列表中添加条目：

```r
.builtin_palettes <- list(
  # ... 已有内容 ...
  
  # 你新增的 — 仿照已有格式
  "MyJournal" = c("#AABB11", "#BBCC22", "#CCDD33", "#DDEE44")
)
```

### 3.5 颜色解析优先级

```
用户直接传颜色向量 → 自定义注册 set_palette → 内置库 → RColorBrewer → paletteer → 回退 NPG
```

这意味着你可以在任何层级覆盖：比如注册一个和内置同名的色板会优先使用你的版本。

---

## 4. 主题管理

### 4.1 永久添加新期刊主题

在 `theme_system.R` 的 Layer 3 区域添加：

```r
#' 你的新期刊风格
#'
#' @description 一段描述
#' @inheritParams theme_pub_base
#' @export
theme_my_journal <- function(base_size = 12, ...) {
  theme_pub_base(base_size, border = TRUE, grid = TRUE) +     # 选择基石参数
    theme(
      # 你的特殊设置
      legend.position = "top",
      panel.grid.major = element_line(color = "grey95"),
      # ...
    ) +
    theme(...)
}
```

**不需要写任何 ggplot2 之外的代码。** 这就是一个纯 theme() 的包装函数。

### 4.2 核心的三个参数

`theme_pub_base()` 提供了三个切换开关：

| 参数 | 默认值 | 效果 |
|---|---|---|
| `base_size` | 12 | 所有文字同步缩放 |
| `border` | TRUE | TRUE = 全黑框 / FALSE = 坐标轴线 |
| `grid` | TRUE | TRUE = 虚线网格 / FALSE = 无线条 |

新主题只需要选择这三个参数的值，然后在 `theme(...)` 中覆盖需要微调的部分。

---

## 5. 在上层绘图函数中使用

### 5.1 标准模式 (开发新函数时)

```r
# 在函数文件开头
source("../../utils/helpers.R")
source("../../utils/palette_system.R")

# 在函数内部
my_plot <- function(data, ..., palette = NULL, theme_use = NULL) {
  p <- ggplot(data, aes(...)) +
    geom_xxx() +
    scale_fill_pub_d(palette) +    # ← 一行对接配色
    labs(...)

  if (is.function(theme_use)) {
    p <- p + theme_use(...)        # ← 一行对接主题
  }
  return(p)
}
```

详细模板见 `FUNCTION_TEMPLATE.R`。

### 5.2 在 compare_plot_adaptive 中的用法

```r
compare_plot_adaptive(data, "group", "value",
  palette   = "JAMA",            # 或自定义颜色向量
  theme_use = theme_nature,       # 或 theme_pub_base()
  add_stat  = "wilcox.test"
)
```

### 5.3 在 compare_plot (稳定版) 中的用法

```r
compare_plot(data, value.var = "value", group.by = "group",
  palette   = "NPG",
  theme_use = theme_pub_base
)
```

---

## 6. 工作流速查

### 日常分析：最小启动

```r
source("R/utils/helpers.R")
source("R/utils/palette_system.R")

ggplot(df, aes(x, y, fill = group)) +
  geom_xxx() +
  scale_fill_pub_d("NPG")          # 改这里换配色
```

### 投不同期刊：改两行

```r
# 投 Nature
p + scale_fill_pub_d("NPG") + theme_nature()

# 改投 Cell
p + scale_fill_pub_d("JCO") + theme_cell()

# 改投 JAMA
p + scale_fill_pub_d("JAMA") + theme_jama()
```

### 发现好看的配色：存起来

```r
set_palette("from_paper_x", c("#...", "#...", "#..."))
preview_palette("from_paper_x")      # 看一眼
```

### 新分析项目启动：批量注册

```r
source("R/utils/palette_system.R")

set_palette(list(
  "proj_phyla"  = c("#E64B35", "#4DBBD5", "#00A087", "#F39B7F", "#8491B4"),
  "proj_groups" = c("#BC3C29", "#0072B5", "#E18727")
))
set_default_palette("proj_groups")
```

---

## 7. 文件修改指南 (给 Claude 用)

当需要改动这个系统时，按以下顺序操作：

| 需求 | 改动文件 | 改什么 |
|---|---|---|
| 加一个内置色板 | `palette_system.R` | 在 `.builtin_palettes` 列表末尾加一个命名向量 |
| 加一个新主题 | `theme_system.R` | 在 Layer 3 区域仿照 `theme_nature()` 写一个新函数 |
| 改默认字号/边框风格 | `theme_system.R` | 修改 `theme_pub_base()` 中的 theme() 参数 |
| 加新的 scale 类型 | `palette_system.R` | 仿照 `scale_fill_pub_d()` 写新 scale 函数 |
| 只在这个脚本用一套颜色 | 不改系统文件 | 用 `set_palette()` 运行时注册 |
| 修改颜色解析优先级 | `palette_system.R` | 修改 `get_colors()` 中的 if-else 链 |
| 新绘图函数开发 | 复制 `FUNCTION_TEMPLATE.R` | 填入业务逻辑，保留 palette/theme_use 参数 |

**原则**: 不重复定义调色板列表。如果发现自己在复制颜色向量到新文件，应该改用 `set_palette()` 或更新 `palette_system.R`。

---

## 8. 依赖项

| 文件 | 需要的 R 包 |
|---|---|
| `helpers.R` | 无 |
| `palette_system.R` | ggplot2, scales | 
| `theme_system.R` | ggplot2 |
| `panel_fix.R` | ggplot2, grid |
| `export_pptx.R` | officer, rvg, ggplot2, grid（依赖 `panel_fix.R`） |

可选增强（安装了就能用更多色板）：
- `RColorBrewer` — 自动识别其所有色板名
- `paletteer` — 自动识别 ggsci/viridis/wesanderson 等 50+ 包的色板

---

## 9. 速查卡片

```
# 加载
source("R/utils/helpers.R")
source("R/utils/palette_system.R")
source("R/utils/theme_system.R")

# 配色
scale_fill_pub_d("NPG")       # 离散 fill
scale_color_pub_d("NPG")      # 离散 color
scale_fill_pub_c("RWB")       # 连续 fill

# 主题
theme_pub_base()              # 基础
theme_pub_stat()              # 统计图
theme_nature()                # Nature
theme_cell()                  # Cell
theme_jama()                  # JAMA
theme_science()               # Science

# 管理
list_palettes()               # 列出色板
preview_palette("NPG")        # 预览
set_palette("name", colors)   # 注册
set_default_palette("name")   # 设默认
get_colors("NPG", n = 5)      # 直接取颜色向量
```
