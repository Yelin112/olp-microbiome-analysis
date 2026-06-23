# pcoa_plot() 使用说明



## 加载函数

```r
source("e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/lib/pcoa_plot.R")
```

依赖包：
```r
library(ggplot2)           # 必需
install.packages("ggExtra") # density 边际图
install.packages("aplot")   # boxplot 边际图
install.packages("RColorBrewer") # 非 NPG 色板
```

---

## 函数签名

```r
pcoa_plot(
  data,                        # data.frame，含坐标列和分组列
  x,                           # x 轴列名（字符串）
  y,                           # y 轴列名（字符串）
  group.by,                    # 分组列名（字符串）
  marginal  = "density",       # 边际图类型："density" | "boxplot" | "none"
  add_ellipse   = TRUE,        # 是否绘制置信椭圆
  ellipse_level = 0.95,        # 椭圆置信水平
  ellipse_type  = "norm",      # 椭圆统计类型："norm" | "t" | "euclid"
  ellipse_alpha = 0.15,        # 椭圆填充透明度
  stat_text  = NULL,           # 统计注释文字（字符向量，每元素一行）
  stat_x     = Inf,            # 注释 x 坐标（Inf = 右边缘）
  stat_y     = Inf,            # 注释 y 坐标（Inf = 上边缘）
  stat_hjust = 1,              # 注释水平对齐（1 = 右对齐）
  stat_vjust = 1.5,            # 注释垂直偏移
  stat_size  = 3.5,            # 注释字体大小
  stat_italic = TRUE,          # 注释是否斜体
  sig_letters = NULL,          # 显著性字母，仅 boxplot 有效（命名向量）
  point_size  = 2.5,           # 散点大小
  point_alpha = 0.85,          # 散点透明度
  point_shape = 16,            # 散点形状（16=实心圆，21=描边圆）
  palette     = "NPG",         # 配色方案
  marginal_size = 0.25,        # 边际图占主图的比例
  xlab  = NULL,                # x 轴标签（NULL 则用列名）
  ylab  = NULL,                # y 轴标签（NULL 则用列名）
  title = NULL,                # 图标题
  theme_use = theme_classic,   # ggplot2 主题函数
  ...                          # 传给 geom_point() 的其他参数
)
```

---

## 返回值

| `marginal` | 返回类型 | 说明 |
|---|---|---|
| `"none"` | ggplot | 可继续用 `+` 叠加图层 |
| `"density"` | ggExtraPlot | 支持 `ggsave()` |
| `"boxplot"` | aplot | 支持 `ggsave()` |

---

## 使用示例

### 示例1：最简单用法（density 边际图）

```r
source("r_functions/lib/pcoa_plot.R")

# 假设已有 PCoA 结果的 data.frame
# 列：PC1, PC2, Group（以及其他元数据列）
pcoa_plot(pcoa_df, x = "PC1", y = "PC2", group.by = "Group")
```

---

### 示例2：与 beta_calc() 配合使用（推荐工作流）

```r
source("r_functions/lib/pcoa_plot.R")
source("r_functions/dev/beta_calc/function.R")

# 计算距离 + 降维 + 统计
res <- beta_calc(
  data     = my_matrix,   # 样本×特征矩阵
  meta     = my_meta,     # 样本元数据
  group    = "Group",
  method   = "PCoA",
  dist     = "bray"
)

# 直接绘图，轴标签和统计注释自动填入
pcoa_plot(
  res$coords,
  x         = "PC1",
  y         = "PC2",
  group.by  = "Group",
  xlab      = res$axis_labs[1],   # "PCoA1 (39.9%)"
  ylab      = res$axis_labs[2],   # "PCoA2 (12.7%)"
  stat_text = res$stat_text       # "Adonis: R² = 0.33, P < 0.001"
)
```

---

### 示例3：density 边际图（复现图C风格）

```r
pcoa_plot(
  pcoa_df,
  x         = "PC1",
  y         = "PC2",
  group.by  = "Group",
  marginal  = "density",
  palette   = c("#8B7CB3", "#E5A02A"),     # 自定义颜色
  stat_text = "adonis R² = 0.33; P = 0.001",
  xlab      = "PCoA 1 (60.3%)",
  ylab      = "PCoA 2 (12.7%)"
)
```

---

### 示例4：boxplot 边际图 + 显著性字母

```r
pcoa_plot(
  pcoa_df,
  x           = "PCoA1",
  y           = "PCoA2",
  group.by    = "Group",
  marginal    = "boxplot",
  stat_text   = c("Anosim: R = 0.477, P = 0.001",
                  "Adonis: R² = 0.329, P = 0.001"),
  sig_letters = c(SCNC_A = "b", SCNC_N = "a", SCNC_Y = "b"),
  xlab        = "(PCoA1: 39.99%)",
  ylab        = "(PCoA2: 17.69%)"
)
```

---

### 示例5：不要边际图，返回 ggplot 对象后继续修改

```r
p <- pcoa_plot(
  pcoa_df,
  x        = "PC1",
  y        = "PC2",
  group.by = "Group",
  marginal = "none"      # 返回纯 ggplot
)

# 继续叠加图层
p + 
  labs(title = "Beta Diversity - Bray Curtis") +
  theme(legend.position = "bottom")
```

---

### 示例6：保存图片

```r
# density / boxplot 模式
p <- pcoa_plot(pcoa_df, x="PC1", y="PC2", group.by="Group",
               marginal = "density")
ggsave("pcoa_density.pdf", p, width = 6, height = 5)

# none 模式（纯 ggplot）
p <- pcoa_plot(pcoa_df, x="PC1", y="PC2", group.by="Group",
               marginal = "none")
ggsave("pcoa.pdf", p, width = 5, height = 4)
```

---

## 配色设置

### 使用预设色板

```r
palette = "NPG"     # 默认，Nature 系配色（10色）
palette = "Set1"    # RColorBrewer Set1（9色）
palette = "Set2"    # RColorBrewer Set2（8色）
palette = "Paired"  # RColorBrewer Paired（12色，适合多分组）
palette = "Dark2"   # RColorBrewer Dark2（8色，深色调）
# 所有 RColorBrewer 色板均支持，查看全部：
# RColorBrewer::display.brewer.all()
```

### 自定义颜色向量（按顺序）

```r
# 颜色数量需 >= 分组数，按分组水平顺序对应
palette = c("#4DBBD5", "#E64B35", "#00A087")
```

### 自定义颜色向量（按组名指定，推荐）

```r
# names 与分组水平完全对应时，顺序无关紧要
palette = c(Control   = "#4DBBD5",
            Treatment = "#E64B35")
```

---

## 参数速查

### 统计注释位置调整

| 位置 | `stat_x` | `stat_y` | `stat_hjust` |
|---|---|---|---|
| 右上角（默认）| `Inf` | `Inf` | `1` |
| 左上角 | `-Inf` | `Inf` | `0` |
| 右下角 | `Inf` | `-Inf` | `1` |
| 指定坐标 | 数据值 | 数据值 | 自定 |

### 椭圆类型

| `ellipse_type` | 说明 |
|---|---|
| `"norm"` | 正态分布假设（默认）|
| `"t"` | t 分布（小样本更稳健）|
| `"euclid"` | 固定半径圆形 |

### 散点形状

```r
point_shape = 16   # 实心圆（默认，无描边）
point_shape = 21   # 描边圆（fill=组色，可设 stroke 宽度）
point_shape = 17   # 实心三角
point_shape = 15   # 实心方块
```

使用 shape 21-25（描边形状）时，可通过 `...` 传入 `stroke` 参数：
```r
pcoa_plot(..., point_shape = 21, stroke = 0.5)
```
