# abundance_processing 使用说明

微生物组丰度数据处理流水线：**聚合 → 过滤 → 标准化 → 差异分析**。
四个函数各自独立（可单独调用），也可以按顺序串联，覆盖从 OTU/ASV 原始丰度表
到差异分析结果的完整链路。

| 函数 | 作用 | 文件 |
|---|---|---|
| `aggregate_taxa()` | 按分类学水平聚合特征 (OTU→Genus 等) | `aggregate_taxa.R` |
| `filter_abundance()` | 按流行率/丰度/方差等过滤低质量特征 | `filter_abundance.R` |
| `norm_abundance()` | 标准化/转换 (TSS、CLR、rarefy 等 23 种方法) | `norm_abundance.R` |
| `diff_abundance()` | 组间差异分析 (含效应量) (15 种方法) | `diff_abundance.R` |

---

## 加载函数

```r
source("R/r_functions/lib/abundance_processing/aggregate_taxa.R")
source("R/r_functions/lib/abundance_processing/filter_abundance.R")
source("R/r_functions/lib/abundance_processing/norm_abundance.R")
source("R/r_functions/lib/abundance_processing/diff_abundance.R")
```

四个函数均为纯 base R 实现，无强制依赖；仅部分 `method` 选项需要按需安装对应包
(见各函数文末的依赖表，或直接看函数内的 `requireNamespace()` 报错提示)。

---

## 数据约定

**所有函数统一约定：行 = 特征 (OTU/ASV/Gene/Taxon)，列 = 样本。**

```
              S1    S2    S3   ...
OTU1          10    5     0
OTU2          3     8     12
...
```

- `abund_table` 可以是 `matrix` 或 `data.frame`，行名必须是 Feature ID（`aggregate_taxa`
  要求行名必存在；其他三个函数在行名缺失时会自动补 `Feature_N`，但那样就没法再和
  分类表匹配，因此聚合步骤务必保留原始行名）。
- `group`/`condition` 等分组参数长度必须与样本数（列数）一致，且顺序一一对应。

---

## 流水线总览

```
原始 OTU/ASV 丰度表 (行=Feature ID, 列=样本)
        │
        ▼
aggregate_taxa()     ← 可选。按 Genus/Species 等水平聚合，减少特征数、
        │               去除分类噪音；不需要聚合可跳过此步
        ▼
filter_abundance()   ← 按流行率/丰度/方差过滤低质量特征，减少多重检验负担
        │
        ▼
norm_abundance()     ← 按下游分析目标选择标准化方法 (见下方选择表)
        │
        ▼
diff_abundance()     ← 组间差异分析，输出统一的 EffectSize + P 值格式
```

> 不是每一步都必须做。例如已经是 Genus 水平的数据可以跳过 `aggregate_taxa()`；
> 用 DESeq2/edgeR/ALDEx2 做差异分析时它们自带标准化，`norm_abundance()` 也可以跳过
> （直接把 `filter_abundance()` 输出的**原始计数**传给 `diff_abundance()`）。

---

## 完整链式示例

```r
source("R/r_functions/lib/abundance_processing/aggregate_taxa.R")
source("R/r_functions/lib/abundance_processing/filter_abundance.R")
source("R/r_functions/lib/abundance_processing/norm_abundance.R")
source("R/r_functions/lib/abundance_processing/diff_abundance.R")

# otu_table: 行=OTU ID, 列=样本；tax_table: 行=OTU ID, 列=Kingdom...Species
# group: 与 otu_table 列一一对应的分组向量

# 1. 聚合到 Genus 水平
genus_table <- aggregate_taxa(
  otu_table, tax_table,
  target_rank = "Genus",
  fill_na_with_rank = "Family"   # Genus 缺失时用 Family 名智能填充
)

# 2. 过滤低质量特征 (差异分析专用预设：保留未分类特征，因为它们也可能是真实差异信号)
filt <- filter_for_difftest(genus_table)
genus_clean <- filt$filtered

# 3a. 非参数检验路线：标准化用相对丰度即可，wilcox 本身对深度差异有一定鲁棒性
rel <- norm_abundance(genus_clean, method = "TSS")
res_wilcox <- diff_abundance(rel, group, method = "wilcox")

# 3b. 参数检验路线：CLR 处理组成性约束后再做 t 检验/ANOVA
clr <- norm_abundance(genus_clean, method = "clr")
res_ttest <- diff_abundance(clr, group, method = "t.test")

# 3c. 计数模型路线：DESeq2 自带标准化，直接传原始计数 (跳过 norm_abundance)
res_deseq2 <- diff_abundance(genus_clean, group, method = "DESeq2")

# 4. 查看显著特征
sig <- res_wilcox[res_wilcox$Significance != "ns", ]
head(sig[, c("Feature", "EffectSize", "EffectSize_type", "P.adj", "Significance")])

# 5. 可视化
plot_effect_bar(res_wilcox, top_n = 20, sig_only = TRUE)
plot_volcano(res_deseq2, log2fc_cutoff = 1)   # 仅适用于 log2FC 类效应量的方法
```

---

## 分步用法

### 1. `aggregate_taxa()` — 分类学聚合

```r
aggregate_taxa(
  abund_table,                # 丰度矩阵，行=Feature ID，列=样本
  tax_table,                  # 分类注释表，行=Feature ID，列=分类等级
  target_rank,                # 目标聚合水平，如 "Phylum"/"Genus"/"Species"
  na_action = "unclassified",  # "unclassified"=保留缺失值(重命名) | "remove"=丢弃
  fill_na_with_rank = NULL,    # 用上一级分类名智能填充缺失，如 target="Species", fill="Genus"
  keep_all_abundance = TRUE,   # TRUE=保留丰度表全部特征 | FALSE=只保留与 tax_table 的交集
  verbose = TRUE
)
```

**返回**：`data.frame`，行=聚合后的分类名，列=样本，总丰度守恒（不会漏计任何一个原始特征，
除非显式选择 `na_action = "remove"` 或 `keep_all_abundance = FALSE`）。

```r
genus_table <- aggregate_taxa(otu_table, tax_table, target_rank = "Genus")
```

> 仅在**已知需要跨分类水平合并**（如从 ASV 聚合到 Genus 做属水平比较）时使用；
> 如果测序流程输出的表本身已经是目标水平（比如 16S 分析直接给了 Genus 丰度表），
> 跳过这一步。

---

### 2. `filter_abundance()` — 过滤

```r
filter_abundance(
  abund_table,
  method = "combined",         # combined(推荐) | prevalence | abundance | variance | liberal | custom
  prev_threshold = 0.1,        # 流行率阈值 (至少 10% 样本检出)
  abund_threshold = 1e-5,      # 平均相对丰度阈值
  min_count = 2,                # 最小计数 (仅对计数数据生效)
  filter_unclassified = TRUE,   # 是否移除 "unclassified"/"uncultured" 等未分类特征
  group = NULL,                 # 配合 max_one_group_zero 使用
  verbose = TRUE
)
```

**返回**：`list(filtered, removed, stats, summary)`，实际过滤后的矩阵在 `$filtered`。

**预设封装**（按分析目的选择，等价于调整默认参数组合）：

| 预设函数 | 适用场景 |
|---|---|
| `filter_strict()` | 发表级/列表级分析，流行率≥20%、丰度≥0.01% |
| `filter_lenient()` | 探索性分析，尽量保留特征 |
| `filter_for_difftest()` | 差异分析前置过滤，**保留未分类特征**（可能是真实信号） |
| `filter_for_network()` | 网络分析，流行率≥30%、额外做方差过滤 |

```r
filt <- filter_for_difftest(genus_table)
genus_clean <- filt$filtered

# 检查保留了多少 reads / 特征
filt$summary$pct_features_kept
filt$summary$pct_reads_kept
```

> 用 `plot_filter_threshold()` / `plot_filter_sensitivity()` 可视化辅助选阈值，
> 经验上 `pct_reads_kept` 低于 90% 说明阈值可能过严。

---

### 3. `norm_abundance()` — 标准化/转换

```r
norm_abundance(abund_table, method = "TSS", ...)
```

23 种方法分四大类，**按下游分析目标选择**（详见函数 Roxygen 文档或 `自定义函数.md`
里的完整决策表）：

| 下游分析 | 推荐标准化 |
|---|---|
| α 多样性 (Shannon 等) | `rarefy` / `SRS` |
| β 多样性 (Bray-Curtis) | `TSS` / `rarefy` |
| β 多样性 (Euclidean/PCA) | `clr` / `hellinger` |
| 差异分析 (Wilcoxon/KW) | `TSS` 或不标准化 |
| 差异分析 (t.test/ANOVA) | `clr` / `AST` |
| 差异分析 (DESeq2/edgeR/ALDEx2) | **不需要**，工具自带标准化 |
| 网络分析 (相关性) | `clr` |

```r
rel <- norm_abundance(genus_clean, method = "TSS")
clr <- norm_abundance(genus_clean, method = "clr")
```

`norm_batch()` 可一次跑多种方法便于对比：

```r
results <- norm_batch(genus_clean, methods = c("TSS", "clr", "hellinger"))
```

---

### 4. `diff_abundance()` — 差异分析

```r
diff_abundance(
  abund_table, group,
  method = "wilcox",           # 15 种方法，见下表
  p_adjust_method = "fdr",
  paired = FALSE,               # 仅 wilcox/t.test
  ...
)
```

**方法选择速查**：

| 数据要求 | 方法 |
|---|---|
| 相对丰度/计数均可 | wilcox, t.test, KW, KW_dunn, anova, lm, lefse |
| 相对丰度 [0,1] | betareg |
| 必须是整数计数 | DESeq2, edgeR, ALDEx2_t, ALDEx2_kw, metagenomeSeq, ancombc2, linda |

只支持两组的方法（wilcox, t.test, DESeq2, edgeR, metagenomeSeq, ALDEx2_t）在
`group` 有 ≥3 组时会自动生成所有两两组合分别检验。

**返回**：统一格式的 `data.frame`，核心列 `Feature`/`Comparison`/`Group`/`EffectSize`/
`EffectSize_type`/`P.unadj`/`P.adj`/`Significance`，已按 `|EffectSize|` 降序排列。

```r
res <- diff_abundance(rel, group, method = "wilcox")
sig <- res[res$Significance != "ns", ]

# 辅助函数
abund_summary(genus_clean, group)         # 各组 Mean/SD/Median/SE
plot_effect_bar(res, top_n = 20)          # 效应量条形图
plot_volcano(res_deseq2, log2fc_cutoff=1) # 火山图 (仅 log2FC 类方法)
```

---

## 常见问题

**Q: 一定要按顺序四步都做吗？**
不需要。最常见的精简路径是 `filter_abundance()` → `diff_abundance()`（跳过聚合和标准化，
直接用 DESeq2/edgeR 等自带标准化的方法）。

**Q: `aggregate_taxa()` 输出可以直接喂给 `filter_abundance()`/`norm_abundance()` 吗？**
可以，三者的行列约定完全一致（行=特征，列=样本），无需转换。

**Q: 为什么 `norm_abundance()` 标准化之后还要再跑 `diff_abundance()` 里 DESeq2？**
不要这样做。DESeq2/edgeR/ALDEx2/metagenomeSeq/ancombc2/linda 这几个方法要求**原始整数计数**，
且自带标准化逻辑，如果先用 `norm_abundance()` 转换过（尤其是 CLR、TSS 这类非整数输出），
会破坏这些工具的内部假设。这几个方法应直接接在 `filter_abundance()` 之后。
