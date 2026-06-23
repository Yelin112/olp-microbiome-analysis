# enrich_plot 开发者笔记

**版本**：v1.1
**创建日期**：2026-04-24
**最后更新**：2026-05-25
**开发环境**：R + ggplot2 生态，面向微生物组/转录组富集分析可视化

---

## 一、项目目标与背景

为微生物组 OLP 项目开发的富集分析（KEGG/GO）可视化工具集。核心需求：

- 主要输入来源是 clusterProfiler 的 `enrichResult`/`gseaResult` 对象
- 同时支持微生物组分析工具的非标准输出（KO 编号、自定义数据库、OTU 标签）
- **ID 转换不在范围内**：基因/特征 ID 原样展示，不做 symbol↔ENTREZID 等转换
- 提供多种发表级图形样式，参数化程度高，方便一行调用

---

## 二、目录结构

```
enrich_plot/
├── load_enrich_plot.R        # 一键加载入口（见"陷阱"部分）
├── README.md                 # 用户使用说明
├── DEVELOPER_NOTES.md        # 本文件
│
├── R/
│   ├── std_format.R          # 核心数据层（S3 转换 + 数据准备函数）
│   ├── registry.R            # 注册表 + plot_enrich() 主函数
│   └── styles/               # 各可视化样式（每个文件独立，可单独 source）
│       ├── bubble.R
│       ├── bar.R
│       ├── lollipop.R
│       ├── combined.R
│       ├── point_bar.R
│       ├── scatter.R
│       └── radial.R
│
├── tests/                    # 各样式测试脚本（在 RStudio 中点 Source 运行）
│   ├── test_std_format.R
│   ├── test_bubble.R
│   ├── test_bar_lollipop.R
│   ├── test_combined.R
│   ├── test_point_bar.R
│   ├── test_scatter.R
│   └── test_radial.R
│
└── _legacy/                  # 归档的旧代码（plot_enrich_trinity 等）
```

---

## 三、架构设计

### 3.1 数据流

```
用户输入
(enrichResult / gseaResult / data.frame / named list)
        ↓
  as_enrich_std()          ← S3 泛型，统一转为 EnrichStd 格式
        ↓
  prep_enrich_data()       ← 单组：筛选/排序/top_n
  prep_enrich_dual()       ← 双组：合并+方向标记（bar/lollipop/scatter 用）
  prep_enrich_combined()   ← 多分类：合并+category factor（combined/point_bar/radial 用）
        ↓
  plot_enrich_xxx()        ← 各样式绘图函数
        ↓
  ggplot / patchwork 对象
```

### 3.2 标准格式 EnrichStd

继承自 `data.frame` 的 S3 类，包含以下列：

| 列名 | 类型 | 必须 | 说明 |
|---|---|---|---|
| `term_id` | chr | 是 | 通路唯一标识（任意格式） |
| `term_name` | chr | 是 | 通路显示名称 |
| `gene_ids` | chr | 是 | "/" 分隔的基因列表 |
| `p_value` | dbl | 是 | 原始 p 值 |
| `p_adjust` | dbl | 是 | 校正 p 值 |
| `gene_count` | int | 是 | 富集基因数 |
| `bg_count` | int | 否 | 背景基因总数 |
| `gene_ratio` | dbl | 否 | gene_count / 输入基因总数 |
| `rich_factor` | dbl | 否 | gene_count / bg_count |

原始数据的其余列会原样保留在末尾（不会丢失）。

### 3.3 注册表机制

```r
# 每个样式文件末尾自注册
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot("bubble", plot_enrich_bubble, "描述", overwrite = TRUE)
}

# 主函数通过注册表调度
plot_enrich(data, type = "bubble", ...)
# → 查注册表 → 调用 plot_enrich_bubble(data, ...)
```

这意味着：加载顺序必须是 `registry.R` 先于各样式文件，否则自注册会静默跳过。

### 3.4 x 轴变量统一计算（`.compute_x_col()`）

所有样式通过同一个内部函数将 `x_var` 字符串转换为实际列值：

```r
.compute_x_col <- function(df, x_var = "p_adjust") {
  if (x_var == "gene_ratio") { ... return(df$gene_ratio) }
  if (x_var == "gene_count") { return(as.numeric(df$gene_count)) }
  -log10(df$p_adjust)   # 默认
}
```

定义在 `std_format.R` 中，各样式文件直接调用，无需各自实现转换逻辑。

---

## 四、关键设计决策（及背后原因）

### 4.1 为什么用 S3 而不是 R5/R6

R 脚本环境，不打包成 package，S3 足够轻量，且与 data.frame 继承兼容，可直接传给 ggplot。

### 4.2 列名自动识别策略（Method B）

`.ENRICH_COL_CANDIDATES` 定义每个标准列的候选列名列表，优先匹配，找不到时要求用户提供 `col_map`。额外列（非标准列）自动附加保留。

**GeneRatio 特殊处理**：clusterProfiler 的 GeneRatio 是 "k/n" 字符串，在 `as_enrich_std.data.frame()` 中有专门的解析逻辑，转为 numeric 后存入 `gene_ratio` 列。

### 4.3 `prep_enrich_combined()` 的 `category_order` 默认逆序

默认 `rev(all_cats)` 是为了让 ggplot y 轴"从下到上"的排列与输入顺序一致（ggplot factor level 1 在底部）。**这是一个反直觉的设计**，后续如果遇到顺序混乱问题首先排查这里。

### 4.4 `point_bar` 的左右面板分离

最初尝试共用 x 轴（左负右正），但 GeneRatio（0-20%）和 -log10(p.adjust)（1-3）量级差异大，视觉严重失衡。最终改为 patchwork 拼图，每侧独立坐标轴，解决了这个问题。

### 4.5 `geom_line(orientation = "y")` 的顺序依赖

ggplot2 的 `geom_line` 按数据行顺序连接点，不是按 factor 层级顺序。`point_bar.R` 中在创建 factor 之前必须先对数据按 `within_id` 排序，否则折线会乱连。这个 bug 曾经存在过，已在 v1.0 修复（`d <- d[order(d$within_id), ]`）。

### 4.6 `combined.R` 左侧面板位置必须相对 `xaxis_max` 缩放

左侧气泡和分类色块的 x 坐标（`left_bubble`、`left_rect_l/r`）必须按 `xaxis_max` 等比缩放，不能使用硬编码绝对值。

**原因**：x 轴可能是 -log10(p.adjust)（典型范围 1-5）或 GeneRatio（典型范围 0.02-0.2），量级相差 10-100 倍。若使用绝对值（如 `left_bubble = -0.8`），切换 x_var 后左侧面板要么完全消失（GeneRatio 模式），要么严重压缩条形图空间（p_adjust 模式）。

当前实现：
```r
xaxis_max   <- max(plot_data$.x_val, na.rm = TRUE) * 1.1
left_bubble <- -bubble_width * xaxis_max
left_rect_r <- -(bubble_width + rect_width * 0.4) * xaxis_max
left_rect_l <- -(bubble_width + rect_width * 1.4) * xaxis_max
```

`bubble_width` 和 `rect_width` 控制左侧面板占总宽的比例，减小这两个值可给条形图留更多空间。

### 4.7 `load_enrich_plot.R` 在嵌套 source 中路径检测失效

`load_enrich_plot.R` 内部用 `sys.frame(1)$ofile` 获取自身路径，在被另一个脚本 `source()` 时，该表达式返回调用方脚本的路径，导致所有子文件路径错误。

**规避方法**：在可视化脚本中直接以绝对路径手动 source，不经过 `load_enrich_plot.R`：
```r
enrich_plot_dir <- "绝对路径/enrich_plot"
source(file.path(enrich_plot_dir, "R/std_format.R"))
source(file.path(enrich_plot_dir, "R/registry.R"))
for (.f in list.files(file.path(enrich_plot_dir, "R/styles"),
                      pattern = "\\.R$", full.names = TRUE)) source(.f)
```

---

## 五、各样式技术要点

### bubble.R
- `geom_point(shape = 21)` 让填充色和边框色分离，实现描边效果
- x_var 自动降级：`gene_ratio → rich_factor → gene_count`，有 warning
- **v1.1 变更**：渐变色参数从 `color_low`/`color_high` 两个参数合并为单一 `palette` 向量，通过 `colorRampPalette(palette)(50)` 支持任意数量的颜色节点，兼容 RColorBrewer 调色板
- **v1.1 变更**：y 轴显示顺序改为按 x_var 降序（最大值在顶部），原实现按 p_adjust 排序与 x 轴变量无关

### bar.R / lollipop.R
- `.is_dual_input(data)` 检测是否双组输入（定义在 std_format.R）
- 双向模式：y 轴用 `geom_text` 替代轴标签，避免 y 轴文字被中间遮挡
- `scale_x_continuous(labels = abs)` 让负数坐标显示为正数
- **v1.1 变更（bar）**：新增 `color_var` 参数，单组模式下可让颜色渐变映射到独立于 x 轴的变量（如 x = GeneRatio，颜色 = p.adjust，实现双变量同时展示）

### scatter.R
- 唯一纯双向图（不支持单组），强制要求 length-2 named list
- 自定义 y 轴刻度线：用 `geom_segment` 在 x=0 位置绘制，避免标准 y 轴被 `geom_vline` 遮挡
- **v1.1 变更**：大幅扩展参数 — `x_var`、`y_var`、`size_var`（均可配置映射列）；`point_shape`、`stroke_color`、`stroke_width`（点外观）
- **shape dispatch 逻辑**：`is_filled <- point_shape %in% 21:25`。实心点（0-20）只有 `colour` 通道，空心点（21-25）有独立 `fill` + `colour` 通道，两者对应不同的 `aes()` 和 `scale_*_manual()`，分支处理：
  ```r
  if (is_filled) {
    # aes(fill = group) + scale_fill_manual + color = stroke_color（固定值）
  } else {
    # aes(colour = group) + scale_colour_manual
  }
  ```

### combined.R
- 依赖 `gground` 包（`geom_round_col`, `geom_round_rect`），未安装则 `use_round_geom = FALSE` 退化
- 依赖 `ggprism`（`theme_prism()`），未安装则 `use_prism_theme = FALSE` 退化
- gene label 位置（bar 下方）在某些数据下可能与 x 轴重叠，已知问题，用户用 Illustrator 后处理
- **v1.1 变更**：新增 `x_var` 参数；左侧面板所有位置改为相对 `xaxis_max` 缩放（见 4.6 节）；`bubble_width`/`rect_width` 默认值缩小（0.32/0.20 → 0.20/0.13），为条形图留更多空间

### point_bar.R
- 依赖 `patchwork`（必须）
- `x_left_scale`：gene_ratio 默认×100 转百分比，其余×1
- 每个分类生成 `p_left + p_right`，`wrap_plots()` 拼网格

### radial.R
- 依赖 `ggraph`, `tidygraph`, `patchwork`（必须），`ggnewscale`, `stringr`（可选）
- 图结构：root("pathway") → category(level=1) → term(level=2)
- 内外圈通过 `patchwork::inset_element()` 叠合，非真正的坐标叠加
- `inner_radius` 控制极坐标内圆大小，默认 `-max(Count)*20`，需根据数据调整

---

## 六、已知问题与局限

| 问题 | 严重性 | 说明 |
|---|---|---|
| combined gene label 遮挡 | 低 | bar 下方基因名可能与 x 轴重叠，暂由用户 Illustrator 修正 |
| 测试无 assertion | 低 | 测试脚本只验证"不报错"，不验证输出内容正确性 |
| gseaResult 未测试 | 低 | `as_enrich_std.gseaResult()` 有实现但无 mock 测试 |
| GeneRatio 仅精确匹配 | 低 | 只识别 `"GeneRatio"` 列名，其他大小写变体不识别 |
| radial inner_radius 需手调 | 中 | 自动计算结果在真实数据上可能仍需调整 |
| load_enrich_plot.R 嵌套失效 | 中 | 在其他脚本中 source 时路径检测失败，需手动 source（见 4.7 节） |

---

## 七、如何新增一种图形样式

1. 在 `R/styles/` 新建 `your_style.R`
2. 文件顶部加 `library(ggplot2)`
3. 实现 `plot_enrich_your_style(data, col_map = NULL, ...)` 函数，内部调用 `as_enrich_std()` 或 `prep_enrich_*()` 处理数据
4. 文件末尾注册：
```r
if (exists("register_enrich_plot", mode = "function")) {
  register_enrich_plot("your_style", plot_enrich_your_style, "描述", overwrite = TRUE)
}
```
5. 在 `tests/` 新建 `test_your_style.R`，参考其他测试文件的 mock 数据结构

### 数据准备参考

```r
# 单组
std <- as_enrich_std(data, col_map = col_map)
std <- prep_enrich_data(std, top_n = top_n, filter_padj = filter_padj)

# 双组（bar/lollipop 模式）
plot_data <- prep_enrich_dual(data, top_n = top_n, col_map = col_map, x_var = x_var)
# 新增列：group(factor), direction(+1/-1), x(按 x_var 计算×direction)

# 多分类（combined/point_bar 模式）
plot_data <- prep_enrich_combined(data, top_n = top_n, col_map = col_map)
# 新增列：category(factor), index(行序号)
```

---

## 八、未来开发方向

- **热图**（`heatmap`）：term × gene 矩阵，`expand_gene_ids()` 已提供长格式展开工具
- **网络图**（`cnet`）：基因-通路关联网络，需要 `igraph` / `ggraph`
- **UpSet 图**：多组通路集合交集，需要 `UpSetR` / `ComplexUpset`
- **testthat 测试框架**：当前测试无断言，迁移到 testthat 可验证输出正确性
- **gseaResult 专用样式**：GSEA running score 图，目前数据层支持但无对应样式
- **修复 load_enrich_plot.R 路径检测**：改用 `tryCatch(normalizePath(sys.frames()[[...]]$ofile))` 或要求调用方传入显式路径参数

---

## 九、v1.0 → v1.1 变更记录

### 破坏性变更

- `bubble.R`：`color_low`/`color_high` 两个参数移除，替换为单一 `palette` 向量。迁移：`color_low = A, color_high = B` → `palette = c(A, B)`

### 新增参数

| 文件 | 参数 | 说明 |
|---|---|---|
| bubble.R | `palette` | N 色渐变向量，替代原 color_low/color_high |
| bar.R | `x_var`, `x_label`, `color_var`, `color_label` | x 轴变量及独立颜色映射变量 |
| lollipop.R | `x_var`, `x_label` | x 轴变量 |
| scatter.R | `x_var`, `y_var`, `size_var`, `point_shape`, `point_size`, `stroke_color`, `stroke_width` | 全面可配置的轴/点映射 |
| combined.R | `x_var`, `base_size` | x 轴变量；基础字号 |

### 内部修复

| 文件 | 修复内容 |
|---|---|
| std_format.R | 新增 `.compute_x_col()` 统一 x 轴计算逻辑 |
| std_format.R | `prep_enrich_combined()` 修复 GO/KEGG 列数不同导致 rbind 报错 |
| combined.R | 左侧面板位置从硬编码改为相对 xaxis_max 缩放 |
| bubble.R | y 轴显示顺序改为按 x_var 降序（最大值在顶） |

---

## 十、与 AI 协作开发说明

本系统由 Claude（Anthropic）辅助开发完成。后续使用 AI 修改时的建议：

**提供上下文的方式**：
- 直接粘贴目标文件内容（AI 不会自动读取文件历史）
- 说明"当前实现是 X，希望改为 Y，原因是 Z"
- 指出报错的具体行号和错误信息

**容易踩坑的地方**（告诉 AI 时重点说明）：
- `prep_enrich_combined()` 的 `category_order` 默认逆序
- `geom_line` 连接顺序依赖数据行顺序，不是 factor 顺序
- `as_enrich_std.data.frame()` 有 GeneRatio 字符串解析逻辑，enrichResult 方法走的是不同分支
- 加载顺序：`std_format.R` → `registry.R` → 各样式文件
- `combined.R` 左侧面板位置必须相对 `xaxis_max` 缩放，不能用绝对值
- `scatter.R` 的 filled shape（21-25）和 solid shape 使用不同 aes 通道，修改 point 相关代码时注意两个分支都要改

**测试验证节奏**：
每个样式独立开发 → 在 RStudio 点 Source 运行对应测试文件 → 确认 PNG 输出正常 → 再集成到 registry
