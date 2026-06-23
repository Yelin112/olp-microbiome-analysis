# beta_calc() 使用说明

## 加载函数

```r
source("e:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/lib/beta_calc/beta_calc.R")
```

依赖包：
```r
library(vegan)                       # 必需
# 可选（微生物组专属距离）
BiocManager::install("phyloseq")
```

---

## 函数定位

通用多元距离计算与排序函数，适用于任何"样本 × 特征"矩阵：

| 数据类型 | 推荐距离 | 推荐方法 |
|---|---|---|
| 微生物组 OTU/ASV 表 | `"bray"` | `PCoA` |
| 转录组（log-CPM/VST）| `"euclidean"` | `PCA` |
| 代谢组、蛋白质组 | `"euclidean"` / `"bray"` | `PCoA` / `PCA` |
| 生态群落多度矩阵 | `"bray"` / `"horn"` | `PCoA` / `NMDS` |
| 二元（有/无）矩阵 | `"jaccard"` | `PCoA` |

输出可直接对接 `pcoa_plot()`。

---

## 函数签名

```r
beta_calc(
  data     = NULL,      # 样本×特征矩阵（通用输入，样本为行）
  otu      = NULL,      # 特征×样本 OTU 表（微生物组，自动转置）
  tree     = NULL,      # 进化树（phylo 对象，UniFrac 距离需要）
  ps       = NULL,      # phyloseq 对象（优先级最高）
  meta     = NULL,      # 样本元数据 data.frame，行名为样本名
  group    = NULL,      # 分组列名（字符串，必填）
  method   = "PCoA",    # 排序方法："PCoA" | "NMDS" | "PCA"
  dist     = "bray",    # 距离类型（PCA 时忽略，自动用欧氏距离）
  ndim     = 2,         # 保留维度数，默认 2
  stat     = "adonis2", # 统计检验："adonis2" | "anosim" | "mrpp" | "none"
                        # 可传向量同时使用多种方法
  pairwise = FALSE,     # 是否进行两两组间比较
  scale    = FALSE,     # PCA 时是否标准化特征（scale.=TRUE）
  seed     = 42,        # 随机种子（NMDS 迭代和置换检验）
  verbose  = TRUE       # 是否打印过程信息
)
```

---

## 三种输入方式

### 方式1：通用矩阵（样本×特征）

```r
# 转录组、代谢组等，样本为行，特征（基因/代谢物）为列
res <- beta_calc(
  data  = expr_matrix,   # nrow=样本数，ncol=特征数
  meta  = sample_meta,   # 行名与 expr_matrix 行名对应
  group = "Group"
)
```

### 方式2：OTU 表（特征×样本）

```r
# 微生物组标准格式：特征（OTU/ASV）为行，样本为列
# 函数内部自动转置
res <- beta_calc(
  otu   = otu_table,     # nrow=OTU数，ncol=样本数
  meta  = sample_meta,
  group = "Group",
  dist  = "bray"
)
```

### 方式3：phyloseq 对象

```r
# 同时包含 OTU 表、metadata、进化树
# 支持 unifrac、wunifrac 等需要进化树的距离
res <- beta_calc(
  ps    = ps_object,
  group = "Group",
  dist  = "wunifrac"
)
```

> 优先级：`ps` > `otu` > `data`（三者同时提供时，ps 优先）

---

## 返回值

```r
res$coords      # data.frame：降维坐标 + 分组列，直接传给 pcoa_plot()
res$dist_mat    # dist 对象：距离矩阵，可供下游复用
res$stat_table  # data.frame：整体统计检验结果
res$pair_table  # data.frame：两两比较结果（pairwise=FALSE 时为 NULL）
res$stat_text   # 字符向量：格式化注释，直接传给 pcoa_plot(stat_text=)
res$axis_labs   # 字符向量：格式化轴标签，如 c("PCoA1 (39.9%)", "PCoA2 (12.7%)")
res$method      # 字符串：使用的排序方法
res$dist_method # 字符串：使用的距离类型
res$stress      # 数值：NMDS stress（非 NMDS 时为 NULL）
```

`stat_table` 列结构：

| 列名 | 说明 |
|---|---|
| `method` | 检验方法名（adonis2/anosim/mrpp）|
| `statistic_name` | 统计量名称（R2/R/delta）|
| `statistic` | 统计量值 |
| `p_value` | P 值 |

`pair_table` 在 `stat_table` 基础上额外增加 `group1`、`group2` 列。

---

## 使用示例

### 示例1：微生物组 PCoA（最常见场景）

```r
source("r_functions/lib/beta_calc/beta_calc.R")

res <- beta_calc(
  data  = otu_relabund,   # 相对丰度表，样本×OTU
  meta  = sample_meta,
  group = "Treatment",
  dist  = "bray",
  method = "PCoA"
)

res$stat_text
# [1] "Adonis: R² = 0.1823, P = 0.001"
```

### 示例2：转录组 PCA

```r
res <- beta_calc(
  data   = log_cpm,       # log-CPM 矩阵，样本×基因
  meta   = sample_meta,
  group  = "Condition",
  method = "PCA",
  scale  = TRUE           # 基因间量纲不同时建议标准化
)
# PCA 模式自动使用欧氏距离做统计检验
```

### 示例3：NMDS + 多种统计方法

```r
res <- beta_calc(
  data   = species_abund,
  meta   = site_meta,
  group  = "Habitat",
  method = "NMDS",
  dist   = "bray",
  stat   = c("adonis2", "anosim")   # 同时跑两种检验
)

res$stat_text
# [1] "Stress = 0.1823"
# [2] "Adonis: R² = 0.2341, P = 0.001"
# [3] "ANOSIM: R = 0.3412, P = 0.001"
```

### 示例4：开启两两比较

```r
res <- beta_calc(
  data     = otu_relabund,
  meta     = sample_meta,
  group    = "Group",
  pairwise = TRUE           # 默认 FALSE，组数多时谨慎开启
)

res$pair_table
#   group1  group2  method  statistic_name  statistic  p_value
#   A       B       adonis2 R2              0.1823     0.001
#   A       C       adonis2 R2              0.2145     0.002
#   B       C       adonis2 R2              0.0934     0.089
```

### 示例5：只计算距离矩阵，跳过统计

```r
res <- beta_calc(
  data  = otu_relabund,
  meta  = sample_meta,
  group = "Group",
  stat  = "none"          # 跳过置换检验（速度快）
)

res$dist_mat    # 距离矩阵，可传给 vegan 其他函数
```

### 示例6：与 pcoa_plot() 完整对接

```r
source("r_functions/lib/beta_calc/beta_calc.R")
source("r_functions/lib/pcoa_plot.R")

# Step 1：计算
res <- beta_calc(
  data     = otu_relabund,
  meta     = sample_meta,
  group    = "Group",
  dist     = "bray",
  stat     = c("adonis2", "anosim"),
  pairwise = TRUE
)

# Step 2：绘图（density 边际）
pcoa_plot(
  res$coords,
  x         = "PC1",
  y         = "PC2",
  group.by  = "Group",
  xlab      = res$axis_labs[1],
  ylab      = res$axis_labs[2],
  stat_text = res$stat_text,
  marginal  = "density"
)

# Step 3：绘图（boxplot 边际 + 显著性字母）
# 从 pair_table 提取显著性字母（需手动整理或用 multcomp）
pcoa_plot(
  res$coords,
  x           = "PC1",
  y           = "PC2",
  group.by    = "Group",
  xlab        = res$axis_labs[1],
  ylab        = res$axis_labs[2],
  stat_text   = res$stat_text,
  marginal    = "boxplot",
  sig_letters = c(A = "a", B = "b", C = "b")
)
```

---

## 支持的距离类型

### 通用（vegan::vegdist）

| 距离 | 适用场景 |
|---|---|
| `"bray"` | 比例/丰度数据，对零值友好（**最常用**）|
| `"euclidean"` | 连续数值，PCA 等价距离 |
| `"jaccard"` | 二元有无数据 |
| `"horn"` | 比例数据，对样本量差异稳健 |
| `"manhattan"` | 绝对差异之和 |
| `"canberra"` | 对低丰度物种敏感 |
| `"gower"` | 混合类型数据 |

更多距离类型见 `?vegan::vegdist`。

### 微生物组专属（需 phyloseq + 进化树）

| 距离 | 说明 |
|---|---|
| `"unifrac"` | 非加权 UniFrac，考虑进化关系 |
| `"wunifrac"` | 加权 UniFrac，考虑相对丰度 |

---

## 注意事项

1. **输入方向**：`data` 参数要求**样本为行**；`otu` 参数接受**特征为行**（OTU 表的常规方向），函数内自动转置
2. **样本对齐**：`data`/`otu` 与 `meta` 通过行名自动匹配；无行名时按行顺序对齐（给出提示）
3. **PCA 的 dist 参数**：传入 `method="PCA"` 时 `dist` 参数被忽略，统计检验自动使用欧氏距离
4. **两两比较的样本量**：当某对分组样本量 < 4 时，adonis2 置换检验可能不稳定，结果仅供参考
5. **NMDS 收敛**：NMDS 为迭代算法，`stress > 0.2` 时结果可靠性下降，建议换用 PCoA
