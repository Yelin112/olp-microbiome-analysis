# enrich_plot：富集分析可视化系统

一套支持多种图形样式的富集分析可视化工具，兼容 clusterProfiler 输出及自定义数据格式。

---

## 快速开始

```r
# 一键加载（推荐，但见"注意事项"中的嵌套 source 陷阱）
source("path/to/enrich_plot/load_enrich_plot.R")

# 查看可用图形类型
list_enrich_plots()

# 直接传入 enrichResult 对象
plot_enrich(kegg_result, type = "bubble")
plot_enrich(kegg_result, type = "bar", top_n = 15)
```

---

## 安装依赖

```r
# 必须
install.packages(c("ggplot2", "patchwork", "ggrepel"))

# radial 图型额外需要
install.packages(c("ggraph", "tidygraph"))

# 可选（增强功能）
install.packages(c("ggnewscale", "stringr"))
```

---

## 支持的输入格式

| 输入类型 | 说明 |
|---|---|
| `enrichResult` | clusterProfiler::enrichKEGG / enrichGO 直接输出 |
| `gseaResult` | clusterProfiler::gseKEGG / gseGO 直接输出 |
| `data.frame` | 标准列名自动识别，非标准列名通过 `col_map` 映射 |
| `named list` | 多分类图型（combined、point_bar、radial）使用 |

**非标准 data.frame 示例：**
```r
plot_enrich(my_df,
  col_map = list(
    term_name = "pathway_name",
    gene_ids  = "gene_list",
    p_adjust  = "FDR"
  )
)
```

**自动识别的列名**（无需 col_map）：
- 通路名：`Description`, `description`, `pathway`, `Name`, `Term` 等
- 基因列：`geneID`, `gene_ids`, `genes`, `core_enrichment` 等
- p 值：`pvalue`, `p.value`, `pval` 等
- 校正 p：`p.adjust`, `padj`, `FDR`, `qvalue` 等
- 基因数：`Count`, `count`, `gene_count`, `overlap` 等

---

## 可视化样式

### 1. bubble — 气泡图（单组）

x 轴为 GeneRatio / Count，颜色映射 -log10(p.adjust)，气泡大小映射基因数。
y 轴通路名称按 x_var 降序排列（最大值在顶部）。

```r
# 基本用法
plot_enrich(kegg_result, type = "bubble", top_n = 15)

# 指定 x 轴变量
plot_enrich(kegg_result, type = "bubble",
  x_var = "gene_ratio"   # 或 "rich_factor", "gene_count"
)

# 自定义两色渐变（颜色顺序：低值→高值）
plot_enrich(kegg_result, type = "bubble",
  palette = c("#d6eaf8", "#1a5276")
)

# 多色渐变（RColorBrewer）
library(RColorBrewer)
plot_enrich(kegg_result, type = "bubble",
  palette = brewer.pal(9, "YlOrRd")
)

# 截取 Spectral 中的几个颜色
plot_enrich(kegg_result, type = "bubble",
  palette = brewer.pal(11, "Spectral")[c(2, 6, 10)]
)
```

### 2. bar — 条形图（单组 / 双向对比）

```r
# 单组，x 轴为 GeneRatio，颜色渐变映射 -log10(p.adjust)
plot_enrich(kegg_result, type = "bar",
  top_n = 15,
  x_var = "gene_ratio"
)

# 单组，x 轴为 GeneRatio，颜色独立映射 p.adjust（双变量展示）
plot_enrich(kegg_result, type = "bar",
  x_var     = "gene_ratio",
  color_var = "p_adjust",
  group_colors = rev(brewer.pal(11, "Spectral"))
)

# 双向对比
plot_enrich(
  list(Up = kegg_up, Down = kegg_down),
  type         = "bar",
  x_var        = "gene_ratio",
  group_colors = c("#009dd3", "#f29c98")
)
```

### 3. lollipop — 棒棒糖图（单组 / 双向对比）

```r
plot_enrich(
  list(Up = kegg_up, Down = kegg_down),
  type         = "lollipop",
  x_var        = "gene_ratio",
  size_var     = "gene_count",    # 点大小映射，NULL 为等大
  group_colors = c("#009dd3", "#f29c98")
)
```

### 4. scatter — 双向散点图

x/y/size 映射变量均可配置，点的填充色与边框色可独立设置。

```r
# 基本用法
plot_enrich(
  list(Up = kegg_up, Down = kegg_down),
  type         = "scatter",
  top_n        = 15,
  group_colors = c(Up = "#3369e7", Down = "#ff6c5f")
)

# 自定义 x/y/size 映射
plot_enrich(
  list(High = kegg_high, Low = kegg_low),
  type         = "scatter",
  x_var        = "gene_count",    # 或 "gene_ratio", "rich_factor"
  y_var        = "p_adjust",      # 或 "p_value"
  size_var     = "gene_ratio",    # 或 "gene_count", NULL（固定大小）
  group_colors = c("#3369e7", "#ff6c5f")
)

# 空心点描边（shape 21-25 支持 fill + 独立边框色）
plot_enrich(
  list(High = kegg_high, Low = kegg_low),
  type         = "scatter",
  point_shape  = 21,
  stroke_color = "black",   # 描边色（shape 21-25 有效）
  stroke_width = 0.5,
  bg_color     = "#f7f7f7",
  group_colors = c("#3369e7", "#ff6c5f")
)
```

### 5. combined — 多分类综合图（GO + KEGG）

```r
plot_enrich(
  list(BP = go_bp, CC = go_cc, MF = go_mf, KEGG = kegg),
  type            = "combined",
  top_n           = 5,
  x_var           = "gene_ratio",   # 或 "p_adjust"（默认）
  category_colors = c(BP = "#437f64", CC = "#d88c51",
                      MF = "#466277", KEGG = "#be5960")
)
```

### 6. point_bar — 分面折线点 + 条形图

```r
plot_enrich(
  list(BP = go_bp, CC = go_cc, MF = go_mf, KEGG = kegg),
  type      = "point_bar",
  top_n     = 5,
  facet_ncol = 2
)
```

### 7. radial — 环形树状图

```r
plot_enrich(
  list(BP = go_bp, CC = go_cc, MF = go_mf, KEGG = kegg),
  type         = "radial",
  top_n        = 8,
  inner_radius = -150    # 负值越大，内圆越大；真实数据常需手动调整
)
```

---

## 常用参数

| 参数 | 适用范围 | 说明 |
|---|---|---|
| `top_n` | 所有 | 每组显示前 N 个 term（按 p_adjust 升序筛选） |
| `filter_padj` | 大部分 | p.adjust 过滤阈值 |
| `x_var` | bubble / bar / lollipop / combined / scatter | x 轴变量：`"p_adjust"`（默认）、`"gene_ratio"`、`"gene_count"` |
| `title` | 所有 | 图标题 |
| `base_size` | 所有 | 基础字号 |
| `group_colors` | 双向 / 单组渐变 | 双向：长度 2 的颜色向量（每组一色）；单组：渐变两端色 |
| `palette` | bubble | N 色渐变向量，支持 RColorBrewer 调色板 |
| `color_var` | bar（单组） | 独立于 x 轴的颜色映射变量，如 `"p_adjust"` |
| `category_colors` | 多分类图 | 分类配色，named vector |
| `category_order` | 多分类图 | 分类排列顺序 |

---

## 保存图形

```r
p <- plot_enrich(kegg_result, type = "bubble")
ggsave("output.pdf", p, width = 8, height = 6)
ggsave("output.png", p, width = 8, height = 6, dpi = 300)
```

---

## 注意事项

- `GeneRatio`（"k/n" 字符串格式）会被自动解析为数值
- `gene_ratio` 不可用时自动降级为 `rich_factor` → `gene_count`，并输出 warning（属正常行为）
- `combined` 样式需要 `gground` 包（圆角图形），不安装则退化为普通矩形；需要 `ggprism` 包（主题），不安装则退化为 `theme_bw()`
- `radial` 图的内圆大小通过 `inner_radius` 参数控制（负值越大，圆孔越大），自动计算值在真实数据上常需手动调整
- **`load_enrich_plot.R` 嵌套 source 陷阱**：若在另一个脚本中 `source("load_enrich_plot.R")`，内部路径检测（`sys.frame(1)$ofile`）会拿到调用方脚本的路径，导致文件找不到。推荐在可视化脚本中直接手动 source 各文件：
  ```r
  enrich_plot_dir <- "绝对路径/enrich_plot"
  source(file.path(enrich_plot_dir, "R/std_format.R"))
  source(file.path(enrich_plot_dir, "R/registry.R"))
  for (.f in list.files(file.path(enrich_plot_dir, "R/styles"),
                        pattern = "\\.R$", full.names = TRUE)) source(.f)
  ```
