# ClusterGVis 完整学习指南

## 目录
1. [包概述](#1-包概述)
2. [核心工作流程](#2-核心工作流程)
3. [函数详解](#3-函数详解)
4. [基础使用教程](#4-基础使用教程)
5. [高级应用场景](#5-高级应用场景)
6. [源码关键设计](#6-源码关键设计)
7. [最佳实践](#7-最佳实践)
8. [常见问题与技巧](#8-常见问题与技巧)

---

## 1. 包概述

### 1.1 功能定位

**ClusterGVis** 是一个R语言包，专门用于基因表达数据的聚类分析和可视化。它的核心优势是：

- 🔬 **多种聚类方法**：支持kmeans、mfuzz、TCseq、WGCNA四种主流方法
- 📊 **丰富的可视化**：基于ComplexHeatmap和ggplot2，提供出版级图表
- 🧬 **生物学注释**：无缝整合GO/KEGG富集分析
- 🔗 **广泛兼容性**：支持bulk RNA-seq、单细胞、WGCNA、拟时序等多种数据
- 🎨 **高度可定制**：提供极其丰富的参数控制和自定义扩展能力

### 1.2 适用场景

- ✅ 时序基因表达数据分析（发育、分化、处理时间序列）
- ✅ 单细胞marker基因可视化
- ✅ WGCNA共表达网络结果展示
- ✅ 拟时序分析结果可视化
- ✅ 需要整合聚类和富集分析的任何场景

### 1.3 安装

```r
# 从Bioconductor安装
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ClusterGVis")

# 或从GitHub安装开发版
devtools::install_github("junjunlab/ClusterGVis")
```

---

## 2. 核心工作流程

ClusterGVis的标准分析流程分为4个步骤：

```r
# 步骤1: 确定最优聚类数
getClusters(obj = exps)

# 步骤2: 执行聚类分析
ck <- clusterData(obj = exps, 
                  clusterMethod = "kmeans", 
                  clusterNum = 8)

# 步骤3: 富集分析（可选）
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05)

# 步骤4: 可视化
visCluster(object = ck, 
           plotType = "both",
           annoTermData = enrich)
```

**数据流示意**：
```
原始表达矩阵 
    ↓ getClusters()
最优k值选择
    ↓ clusterData()
聚类结果对象 (wide.res + long.res + cluster.list)
    ↓ enrichCluster()
富集注释数据
    ↓ visCluster()
最终可视化图表
```

---

## 3. 函数详解

### 3.1 getClusters() - 确定最优聚类数

#### 功能
使用**Elbow方法**（肘部法则）确定最优的聚类数k值。

#### 参数
```r
getClusters(obj = NULL, ...)
```

- `obj`: 输入数据，支持以下类型：
  - `cell_data_set` (Monocle3对象)
  - `matrix` 或 `data.frame` (基因×样本矩阵)
  - `SummarizedExperiment` 对象
- `...`: 传递给预处理函数的额外参数

#### 返回值
一个**ggplot对象**，显示不同k值对应的WSS（组内平方和）。最优k值通常在曲线的"肘部"位置。

#### 实现原理
```r
# 核心实现
factoextra::fviz_nbclust(exp, stats::kmeans, method = "wss") +
  ggplot2::labs(subtitle = "Elbow method")
```

#### 使用示例
```r
library(ClusterGVis)
data(exps)

# 查看最优聚类数
p <- getClusters(obj = exps)
print(p)

# 通过观察图形，假设选择k=8
```

#### 如何判断最优k值
- 寻找WSS下降速度明显放缓的点（肘部）
- WSS在该点之后趋于平缓
- 通常在3-15之间

---

### 3.2 clusterData() - 核心聚类函数

#### 功能
对基因表达数据执行聚类分析，支持4种方法：mfuzz、TCseq、kmeans、wgcna。

#### 完整参数
```r
clusterData(
  obj = NULL,                    # 输入数据
  scaleData = TRUE,              # 是否标准化（z-score）
  clusterMethod = c("mfuzz", "TCseq", "kmeans", "wgcna"),
  TCseqParamsList = list(),      # TCseq额外参数
  kmeansParamsList = list(),     # kmeans额外参数
  object = NULL,                 # WGCNA预计算对象
  minStd = 0,                    # 最小标准差过滤
  clusterNum = NULL,             # 聚类数
  subcluster = NULL,             # 选择特定簇
  ...                            # 额外参数
)
```

#### 四种聚类方法详解

##### 1. mfuzz (模糊C均值聚类)

**适用场景**：时序表达数据，基因可能属于多个模式

**特点**：
- 软聚类：每个基因有membership值（0-1），表示属于某簇的隶属度
- 适合识别重叠的表达模式
- 基于fuzzy c-means算法

**使用示例**：
```r
cm <- clusterData(
  obj = exps,
  clusterMethod = "mfuzz",
  clusterNum = 8,
  minStd = 0.1  # 过滤低变异基因
)

# 查看结果结构
str(cm)
# $wide.res: 包含gene, 各样本表达值, cluster, membership
# $long.res: 长格式，用于可视化
# $cluster.list: 按簇分组的基因列表
# $type: "mfuzz"
```

**membership解释**：
- 值接近1：基因强烈属于该簇
- 值在0.5左右：基因特征模糊，可能介于两个模式之间
- 可以设置membership阈值筛选核心基因

##### 2. TCseq (时序聚类)

**适用场景**：时序表达数据，强调时间依赖性

**特点**：
- 也是基于fuzzy c-means
- 专门为时序数据优化
- 同样提供membership信息

**使用示例**：
```r
ct <- clusterData(
  obj = exps,
  clusterMethod = "TCseq",
  clusterNum = 8,
  TCseqParamsList = list(
    algo = "cm",           # 聚类算法
    standardize = TRUE     # 是否标准化
  )
)
```

**与mfuzz的区别**：
- TCseq更强调时间序列的连续性
- 对时间点顺序更敏感
- 通常用于有明确时序关系的数据

##### 3. kmeans (K均值聚类)

**适用场景**：一般表达数据，不限于时序

**特点**：
- 硬聚类：每个基因只属于一个簇
- 快速、稳定
- 无membership信息
- 最常用的聚类方法

**使用示例**：
```r
ck <- clusterData(
  obj = exps,
  clusterMethod = "kmeans",
  clusterNum = 8,
  kmeansParamsList = list(
    nstart = 10,      # 随机起始点数量
    iter.max = 100    # 最大迭代次数
  )
)
```

**优点**：
- 计算速度快
- 结果易于解释
- 适合大数据集

##### 4. wgcna (加权基因共表达网络)

**适用场景**：已经完成WGCNA分析，需要可视化模块

**特点**：
- 不重新聚类，而是解析WGCNA结果
- 继承模块颜色信息
- 适合大规模共表达分析

**使用示例**：
```r
# 前提：已完成WGCNA分析
library(WGCNA)
net <- blockwiseModules(
  datExpr,
  power = 6,
  TOMType = "unsigned",
  minModuleSize = 30,
  ...
)

# 使用ClusterGVis可视化
ck <- clusterData(
  obj = datExpr,        # 注意：这是过滤后的表达矩阵
  clusterMethod = "wgcna",
  object = net          # WGCNA结果对象
)
```

**重要提示**：
- `obj`参数应该是WGCNA使用的`datExpr`（过滤后的转置矩阵）
- `object`参数是`blockwiseModules()`的输出
- 簇名称会包含WGCNA的模块颜色（如"cluster 1 (200 turquoise)"）

#### 返回值结构

所有方法都返回一个**list**，包含以下元素：

```r
list(
  wide.res = ...,      # 宽格式数据框
  long.res = ...,      # 长格式数据框
  cluster.list = ...,  # 基因列表（按簇分组）
  type = ...,          # 聚类方法名称
  geneMode = "none",   # 预留
  geneType = "none"    # 预留
)
```

**wide.res结构**（以kmeans为例）：
```
         gene    sample1  sample2  ...  cluster
1       Oog4      1.313    1.237   ...    1
2      Psmd9      1.092    1.316   ...    1
...
```

对于mfuzz/TCseq，还包含`membership`列。

**long.res结构**：
```
   cluster  gene      cell_type  norm_value  cluster_name
1    1      Oog4      sample1    1.313       cluster 1 (381)
2    1      Oog4      sample2    1.237       cluster 1 (381)
...
```

**cluster.list结构**：
```r
$C1
[1] "Oog4"   "Psmd9"  "Sephs2" ...

$C2
[1] "Gene1"  "Gene2"  "Gene3"  ...
```

#### 数据预处理选项

##### scaleData参数
```r
# 标准化（推荐）
ck <- clusterData(obj = exps, scaleData = TRUE)  # z-score标准化

# 不标准化（数据已经标准化或使用原始值）
ck <- clusterData(obj = exps, scaleData = FALSE)
```

##### minStd参数
```r
# 过滤低变异基因（减少噪音）
ck <- clusterData(obj = exps, minStd = 0.1)
# 只保留标准差 > 0.1 的基因
```

##### subcluster参数
```r
# 只关注特定的簇
ck <- clusterData(obj = exps, clusterNum = 8, subcluster = c(1, 3, 5))
# 只保留簇1、3、5的结果
```

---

### 3.3 enrichCluster() - 富集分析函数

#### 功能
对聚类结果的每个簇分别进行GO或KEGG富集分析。

#### 完整参数
```r
enrichCluster(
  object = NULL,              # clusterData对象
  type = c("BP", "MF", "CC", "KEGG", "ownSet"),
  TERM2GENE = NULL,           # 自定义集合：term到gene映射
  TERM2NAME = NULL,           # 自定义集合：term到name映射
  OrgDb = NULL,               # 物种注释数据库
  idTrans = TRUE,             # 是否进行ID转换
  fromType = "SYMBOL",        # 输入ID类型
  toType = "ENTREZID",        # 目标ID类型
  readable = TRUE,            # 结果是否转回可读格式
  organism = "hsa",           # KEGG物种代码
  pvalueCutoff = 0.05,        # p值阈值
  topn = 5,                   # 每簇提取的top条目数
  addGene = FALSE,            # 是否添加基因列表
  useInternalData = FALSE,    # KEGG使用本地还是在线数据
  ...
)
```

#### 五种富集类型

##### 1. GO富集 (type = "BP"/"MF"/"CC")

```r
library(org.Mm.eg.db)  # 小鼠
# library(org.Hs.eg.db)  # 人类

enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",              # 生物过程
  # type = "MF",            # 分子功能
  # type = "CC",            # 细胞组分
  pvalueCutoff = 0.05,
  topn = 5
)

head(enrich)
#   group                Description       pvalue    ratio
# 1    C1  lymphocyte differentiation 4.26e-09  42.10526
# 2    C1       T cell differentiation 1.28e-08  36.84211
```

##### 2. KEGG通路 (type = "KEGG")

```r
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "KEGG",
  organism = "mmu",         # 小鼠：mmu，人类：hsa
  pvalueCutoff = 0.05,
  topn = 5,
  useInternalData = FALSE   # TRUE使用本地KEGG.db，FALSE在线查询
)
```

**常用organism代码**：
- `hsa`: Homo sapiens (人)
- `mmu`: Mus musculus (小鼠)
- `rno`: Rattus norvegicus (大鼠)
- `dme`: Drosophila melanogaster (果蝇)
- `cel`: Caenorhabditis elegans (线虫)
- `sce`: Saccharomyces cerevisiae (酵母)

##### 3. 自定义基因集 (type = "ownSet")

```r
# 准备自定义基因集
TERM2GENE <- data.frame(
  term = c("MySet1", "MySet1", "MySet2", "MySet2"),
  gene = c("Gene1", "Gene2", "Gene3", "Gene4")
)

TERM2NAME <- data.frame(
  term = c("MySet1", "MySet2"),
  name = c("Custom Set 1", "Custom Set 2")
)

enrich <- enrichCluster(
  object = ck,
  type = "ownSet",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  idTrans = FALSE,          # 自定义集合通常不需要ID转换
  pvalueCutoff = 0.05
)
```

#### ID转换详解

富集分析通常需要ENTREZID，但基因通常是SYMBOL格式。

```r
# 自动转换（推荐）
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  idTrans = TRUE,
  fromType = "SYMBOL",      # 输入是SYMBOL
  toType = "ENTREZID",      # 转换为ENTREZID
  readable = TRUE           # 结果转回SYMBOL显示
)

# 手动转换（如果数据已经是ENTREZID）
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  idTrans = FALSE,
  fromType = "ENTREZID",    # 已经是ENTREZID
  readable = FALSE
)
```

**常用ID类型**：
- `SYMBOL`: 基因符号（如"TP53"）
- `ENTREZID`: Entrez ID（如"7157"）
- `ENSEMBL`: Ensembl ID（如"ENSG00000141510"）
- `REFSEQ`: RefSeq ID（如"NM_000546"）

#### 返回值结构

```r
head(enrich)
#   group                Description       pvalue    ratio  geneID
# 1    C1  lymphocyte differentiation 4.26e-09  42.10526  Cd3d/Lck/...
# 2    C1       T cell differentiation 1.28e-08  36.84211  Cd3d/Lef1/...
```

- `group`: 簇ID（C1, C2, ...）
- `Description`: 通路/功能描述
- `pvalue`: 富集显著性p值
- `ratio`: 富集比例（基因比例）
- `geneID`: 涉及的基因（如果addGene=TRUE）

#### topn参数灵活使用

```r
# 方式1：所有簇使用相同的topn
enrich <- enrichCluster(object = ck, topn = 5)

# 方式2：每个簇使用不同的topn
enrich <- enrichCluster(object = ck, topn = c(5, 3, 7, 5, 10, 4, 6, 8))
# 簇1取5条，簇2取3条，依此类推
```

#### 非模式生物富集

对于非模式生物，没有OrgDb注释包，可以：

**方法1：使用自定义基因集**
```r
# 从其他来源获取基因-GO映射
# 例如通过BLAST比对到模式生物
TERM2GENE <- read.csv("my_species_GO_mapping.csv")

enrich <- enrichCluster(
  object = ck,
  type = "ownSet",
  TERM2GENE = TERM2GENE,
  idTrans = FALSE
)
```

**方法2：构建自己的OrgDb**
```r
# 使用AnnotationForge包
library(AnnotationForge)
makeOrgPackage(...)
```

---

### 3.4 visCluster() - 可视化核心函数

这是ClusterGVis最复杂也是最强大的函数，有**1582行源码**，提供极其丰富的可视化选项。

#### 三种绘图模式

##### 模式1: 折线图 (plotType = "line")

**适用场景**：初步探索表达趋势

```r
visCluster(
  object = cm,
  plotType = "line",
  lineSize = 0.1,
  lineCol = "grey90",
  addMline = TRUE,        # 添加中位线
  mlineSize = 2,
  mlineCol = "#CC3333",
  ncol = 4                # facet列数
)
```

**特点**：
- 基于ggplot2
- 每个基因一条线
- 支持facet分面展示
- mfuzz/TCseq可以按membership着色

**示例：改变颜色**
```r
# mfuzz结果的membership梯度颜色
visCluster(
  object = cm,
  plotType = "line",
  msCol = c("green", "orange", "red")  # 低-中-高membership
)
```

##### 模式2: 热图 (plotType = "heatmap")

**适用场景**：详细展示表达模式，添加注释

```r
visCluster(
  object = ck,
  plotType = "heatmap",
  column_names_rot = 45,
  show_row_dend = FALSE,
  border = TRUE
)
```

**特点**：
- 基于ComplexHeatmap
- 支持丰富的注释系统
- 可标记基因
- 可添加GO/KEGG注释

**核心参数**：
```r
# 热图颜色控制
htColList = list(
  col_range = c(-2, 0, 2),              # 数值范围
  col_color = c("#08519C", "white", "#A50F15")  # 对应颜色
)

# 簇注释条颜色
ctAnnoCol = ggsci::pal_npg()(8)

# 标记基因
markGenes = c("Oog4", "Psmd9", "Sephs2")
markGenesSide = "right"  # 或 "left"
genesGp = c('italic', 10, NA)  # 字体样式、大小、颜色
```

##### 模式3: 组合图 (plotType = "both")

**适用场景**：最常用，展示热图+趋势

```r
pdf('result.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "right",     # 趋势图位置
  column_names_rot = 45
)
dev.off()
```

**特点**：
- 热图在主体
- 趋势图（折线/箱线图）在左侧或右侧
- 可同时显示多种注释

**布局控制**：
```r
# 趋势在右侧（默认）
visCluster(object = ck, plotType = "both", lineSide = "right")

# 趋势在左侧
visCluster(object = ck, plotType = "both", lineSide = "left")

# 面板大小控制
panelArg = c(
  2,      # 面板大小（单位）
  0.25,   # 间隙
  4,      # 宽度
  "grey90",  # 填充色
  NA      # 边框色
)
```

#### 样本注释系统

##### 单个样本分组
```r
visCluster(
  object = ck,
  plotType = "heatmap",
  sampleGroup = rep(c("A", "B", "C"), each = 2),
  sampleCol = c("A" = "red", "B" = "blue", "C" = "green")
)
```

##### 多个样本注释
```r
library(ComplexHeatmap)

# 准备多个注释
mg1 = rep(c("D1", "D2"), each = 3)
names(mg1) <- colnames(exps)

mg2 = rep(c("E1", "E2", "E3"), each = 2)
names(mg2) <- colnames(exps)

# 创建ComplexHeatmap注释对象
topAnnotation <- HeatmapAnnotation(
  mg1 = mg1,
  mg2 = mg2,
  col = list(
    mg1 = c("D1" = "#C147E9", "D2" = "#FF7000"),
    mg2 = c("E1" = "#54B435", "E2" = "#31C6D4", "E3" = "#D9CB50")
  ),
  gp = gpar(col = "white")
)

# 应用注释
visCluster(
  object = ck,
  plotType = "heatmap",
  heatmapAnnotation = topAnnotation
)
```

##### 样本排序
```r
# 自定义样本顺序
visCluster(
  object = ck,
  plotType = "heatmap",
  sampleOrder = c("t8.cell", "tmorula", "blastocyst", 
                  "zygote", "t2.cell", "t4.cell")
)
```

##### 分割列
```r
visCluster(
  object = ck,
  plotType = "heatmap",
  heatmapAnnotation = topAnnotation,
  columnSplit = rep(c(1, 2, 3), each = 2)  # 分成3组
)
```

#### 行注释系统

行注释在热图主体左侧或右侧添加额外信息。

##### 基本用法
```r
# 准备注释数据（必须与聚类后的基因顺序一致）
df <- ck$wide.res %>% arrange(cluster)

# 创建命名向量
anno_vec <- rep(LETTERS[1:3], c(100, 150, 200))
names(anno_vec) <- df$gene

# 方式1：使用ComplexHeatmap对象（当lineSide="right"时）
library(ComplexHeatmap)
rowAnno <- rowAnnotation(
  anno = anno_vec,
  text = anno_text(sample(letters, nrow(df), replace=TRUE)),
  bar = anno_barplot(runif(nrow(df)))
)

visCluster(
  object = ck,
  plotType = "heatmap",
  rowAnnotationObj = rowAnno
)
```

```r
# 方式2：使用命名列表（当lineSide="left"时）
rowAnnoList <- list(
  anno = anno_vec,
  text = anno_text(sample(letters, nrow(df), replace=TRUE), which = "row"),
  bar = anno_barplot(runif(nrow(df)), which = "row")
)

visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",
  rowAnnotationObj = rowAnnoList
)
```

#### GO/KEGG富集注释

这是ClusterGVis的特色功能，可以直接在热图上展示富集结果。

##### GO注释（基础）
```r
# 先进行富集
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",
  topn = 5
)

# 添加到可视化
pdf('with_GO.pdf', height = 10, width = 10)
visCluster(
  object = ck,
  plotType = "both",
  annoTermData = enrich,
  annoTermMside = "right",    # GO注释在右侧
  column_names_rot = 45
)
dev.off()
```

##### GO注释（高级）
```r
pdf('GO_advanced.pdf', height = 10, width = 12)
visCluster(
  object = ck,
  plotType = "both",
  annoTermData = enrich,
  annoTermMside = "right",
  
  # GO文本样式
  goCol = rep(ggsci::pal_d3()(8), each = 3),  # 文本颜色
  goSize = "pval",                             # 大小基于p值
  # goSize = 3,                                # 或固定大小
  
  # GO框样式
  byGo = "anno_link",                          # 连线样式
  # byGo = "anno_block",                       # 块状样式
  
  # 文本控制
  wordWrap = TRUE,                             # 自动换行
  addNewLine = TRUE,                           # 长文本换行
  termTextLimit = c(10, 18),                   # 文本长度限制
  
  # 框面板样式
  termAnnoArg = c("grey95", "grey50")          # 填充色、边框色
)
dev.off()
```

##### KEGG注释
```r
# 富集KEGG
kegg_enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "KEGG",
  organism = "mmu",
  topn = 3
)

# 添加KEGG注释
pdf('with_KEGG.pdf', height = 10, width = 12)
visCluster(
  object = ck,
  plotType = "both",
  annoKeggData = kegg_enrich,
  annoKeggMside = "right",
  keggCol = rep(ggsci::pal_lancet()(8), each = 3),
  keggSize = "pval",
  byKegg = "anno_link",
  keggAnnoArg = c("grey95", "grey50")
)
dev.off()
```

##### 同时添加GO和KEGG
```r
pdf('GO_and_KEGG.pdf', height = 10, width = 14)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",              # 趋势在左
  
  annoTermData = go_enrich,       # GO在右
  annoTermMside = "right",
  goCol = rep(ggsci::pal_d3()(8), each = 3),
  
  annoKeggData = kegg_enrich,     # KEGG也在右（会并排）
  annoKeggMside = "right",
  keggCol = rep(ggsci::pal_lancet()(8), each = 3)
)
dev.off()
```

##### 为特定簇添加注释
```r
# 只为簇1、3、5添加GO注释
enrich_subset <- enrich %>% filter(group %in% c("C1", "C3", "C5"))

visCluster(
  object = ck,
  plotType = "both",
  annoTermData = enrich_subset,
  subgroupAnno = c(1, 3, 5)  # 只在这些簇添加
)
```

##### 添加条形图展示
```r
pdf('GO_with_bar.pdf', height = 10, width = 12)
visCluster(
  object = ck,
  plotType = "both",
  annoTermData = enrich,
  addBar = TRUE,          # 添加条形图
  barWidth = 8,           # 条宽度
  textbarPos = c(0.8, 0.8)  # 文本位置
)
dev.off()
```

#### 箱线图和点图

```r
# 添加箱线图
pdf('with_boxplot.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "both",
  addBox = TRUE,
  boxCol = ggsci::pal_npg()(8),        # 箱体颜色
  boxArg = c(0.1, "grey50"),           # 宽度、边框色
  addLine = FALSE                      # 移除折线
)
dev.off()

# 添加点图
pdf('with_points.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "both",
  addBox = TRUE,
  addPoint = TRUE,
  pointArg = c(19, "orange", "orange", 1),  # 形状、填充、边框、大小
  addLine = FALSE
)
dev.off()
```

#### 多组趋势线

适用于比较不同实验条件的趋势。

```r
# 假设有6个样本，分成2组，每组3个重复
pdf('multi_group.pdf', height = 10, width = 11)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",
  mulGroup = c(3, 3),                  # 2组，每组3个样本
  mlineCol = c("#E64B35", "#4DBBD5"),  # 每组的颜色
  lgdLabel = c("Control", "Treatment") # 图例标签
)
dev.off()

# 3组的例子
visCluster(
  object = ck,
  plotType = "both",
  mulGroup = c(2, 2, 2),               # 3组
  mlineCol = c("#E64B35", "#4DBBD5", "#00A087"),
  lgdLabel = c("Group A", "Group B", "Group C")
)
```

#### 簇排序

```r
# 自定义簇顺序
visCluster(
  object = ck,
  plotType = "both",
  clusterOrder = c(3, 1, 5, 2, 7, 4, 8, 6)  # 新的簇顺序
)
```

#### 自定义ggplot图形（第7章核心功能）

这是ClusterGVis最强大的扩展功能，允许为每个簇添加任意ggplot2图形。

##### 基本用法
```r
# 为每个簇创建自定义图
gglist <- lapply(1:8, function(i) {
  # 这里可以是任意ggplot2代码
  ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) +
    geom_line() +
    ggtitle(paste("Cluster", i, "custom plot"))
})

names(gglist) <- paste0("C", 1:8)

# 插入到可视化
pdf('with_custom_plots.pdf', height = 10, width = 14)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",
  gglist = gglist,
  ggplotPanelArg = c(2, 0.25, 4, "grey90", NA)  # 控制ggplot面板
)
dev.off()
```

##### 单细胞应用：添加Feature Plot
```r
library(Seurat)

# 准备单细胞数据
markers <- FindAllMarkers(pbmc)
st.data <- prepareDataFromscRNA(pbmc, markers, showAverage = TRUE)

# 为每个细胞类型创建Feature Plot
gglist <- lapply(unique(markers$cluster), function(x) {
  genes <- markers %>% filter(cluster == x) %>% pull(gene)
  FeaturePlot(pbmc, features = genes, ncol = 4)
})
names(gglist) <- paste0("C", unique(markers$cluster))

# 可视化
pdf('scRNA_with_featurePlot.pdf', height = 12, width = 16)
visCluster(
  object = st.data,
  plotType = "both",
  lineSide = "left",
  gglist = gglist,
  column_names_rot = 45
)
dev.off()
```

##### 使用scRNAtoolVis的cornerAxes
```r
library(scRNAtoolVis)

gglist <- lapply(unique(markers$cluster), function(x) {
  genes <- markers %>% filter(cluster == x) %>% pull(gene)
  featureCornerAxes(
    object = pbmc,
    reduction = 'umap',
    features = genes,
    relLength = 0.65,
    relDist = 0.05
  )
})
names(gglist) <- paste0("C", unique(markers$cluster))

visCluster(
  object = st.data,
  plotType = "both",
  gglist = gglist
)
```

#### ComplexHeatmap参数传递

visCluster可以接收ComplexHeatmap::Heatmap()的所有参数。

```r
visCluster(
  object = ck,
  plotType = "heatmap",
  
  # ComplexHeatmap参数
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,
  cluster_columns = FALSE,  # 或 clusterColumns = FALSE
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  heatmap_legend_param = list(
    title = "Expression",
    legend_direction = "vertical",
    legend_width = unit(4, "cm")
  ),
  use_raster = TRUE,        # 大数据集使用光栅化
  raster_quality = 2
)
```

#### 完整示例

```r
# 完整的可视化流程
library(ClusterGVis)
library(org.Mm.eg.db)
library(ComplexHeatmap)

# 1. 数据和聚类
data(exps)
ck <- clusterData(obj = exps, clusterMethod = "kmeans", clusterNum = 8)

# 2. 富集分析
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",
  pvalueCutoff = 0.05,
  topn = 5
)

# 3. 准备标记基因
markGenes <- sample(rownames(exps), 30)

# 4. 准备样本注释
sampleAnno <- HeatmapAnnotation(
  Group = rep(c("A", "B", "C"), each = 2),
  col = list(Group = c("A" = "red", "B" = "blue", "C" = "green"))
)

# 5. 完整可视化
pdf('complete_visualization.pdf', height = 12, width = 14)
visCluster(
  object = ck,
  plotType = "both",
  
  # 热图设置
  htColList = list(
    col_range = c(-2, 0, 2),
    col_color = c("#08519C", "white", "#A50F15")
  ),
  column_names_rot = 45,
  
  # 簇注释
  ctAnnoCol = ggsci::pal_npg()(8),
  
  # 样本注释
  heatmapAnnotation = sampleAnno,
  
  # 标记基因
  markGenes = markGenes,
  markGenesSide = "left",
  genesGp = c('italic', 10, 'black'),
  
  # 趋势图
  lineSide = "left",
  addBox = TRUE,
  boxCol = ggsci::pal_npg()(8),
  addLine = FALSE,
  
  # GO注释
  annoTermData = enrich,
  annoTermMside = "right",
  goCol = rep(ggsci::pal_d3()(8), each = 3),
  goSize = "pval",
  addBar = TRUE,
  
  # 面板设置
  panelArg = c(2, 0.25, 4, "grey90", NA),
  
  # ComplexHeatmap参数
  show_row_dend = FALSE,
  border = TRUE
)
dev.off()
```

---

## 4. 基础使用教程

### 4.1 输入数据要求

ClusterGVis接受多种类型的输入数据：

```r
# 1. 标准表达矩阵（最常用）
# 行 = 基因，列 = 样本
exps <- read.csv("expression_matrix.csv", row.names = 1)
head(exps, 3)
#           sample1  sample2  sample3  sample4
# Gene1     1.313    1.237    1.326    1.262
# Gene2     1.092    1.316    1.174    1.065
# Gene3     0.986    1.201    1.123    1.085

# 2. SummarizedExperiment对象
library(SummarizedExperiment)
se <- SummarizedExperiment(assays = list(counts = exps))
# 函数会自动提取assay(se)

# 3. Monocle的cell_data_set对象（拟时序）
# 见第6章

# 4. WGCNA的结果对象
# 见第4章
```

**数据预处理建议**：
- 归一化：TPM、FPKM、RPKM、RPM等
- 行名必须是基因名称
- 去除全为0或低表达的基因
- 可选：log转换

### 4.2 标准工作流程

#### 步骤1：加载示例数据
```r
library(ClusterGVis)
data(exps)

# 查看数据
dim(exps)  # 3767 genes × 6 samples
head(exps, 3)
#           zygote  t2.cell  t4.cell  t8.cell   tmorula blastocyst
# Oog4   1.3132282 1.237078 1.325978 1.262073 0.6549312  0.2067114
# Psmd9  1.0917337 1.315989 1.174417 1.064756 0.8685598  0.4845448
# Sephs2 0.9859232 1.201026 1.123076 1.084673 0.8878931  0.7174088
```

这是小鼠胚胎发育前着床期的蛋白质表达数据，包含6个时间点。

#### 步骤2：确定聚类数
```r
getClusters(obj = exps)
# 观察Elbow图，选择拐点处的k值
# 假设选择k = 8
```

#### 步骤3：执行聚类
```r
# 使用mfuzz
cm <- clusterData(
  obj = exps,
  clusterMethod = "mfuzz",
  clusterNum = 8
)

# 或使用kmeans
ck <- clusterData(
  obj = exps,
  clusterMethod = "kmeans",
  clusterNum = 8
)
```

#### 步骤4：初步可视化
```r
# 折线图查看趋势
visCluster(object = cm, plotType = "line")

# 热图查看详细模式
visCluster(object = ck, plotType = "heatmap")
```

#### 步骤5：富集分析
```r
library(org.Mm.eg.db)

enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",
  pvalueCutoff = 0.05,
  topn = 5
)

head(enrich)
```

#### 步骤6：完整可视化
```r
pdf('final_result.pdf', height = 10, width = 12)
visCluster(
  object = ck,
  plotType = "both",
  annoTermData = enrich,
  column_names_rot = 45
)
dev.off()
```

### 4.3 自定义颜色方案

```r
# 热图颜色
visCluster(
  object = ck,
  plotType = "heatmap",
  htColList = list(
    col_range = c(-2, 0, 2),
    col_color = c("blue", "white", "red")  # 冷暖色
  )
)

# 簇注释条颜色
library(ggsci)
visCluster(
  object = ck,
  plotType = "heatmap",
  ctAnnoCol = pal_npg()(8)      # Nature配色
  # ctAnnoCol = pal_lancet()(8) # Lancet配色
  # ctAnnoCol = pal_jco()(8)    # JCO配色
)

# 折线图颜色（mfuzz）
visCluster(
  object = cm,
  plotType = "line",
  msCol = c("green", "yellow", "red")  # membership梯度
)
```

### 4.4 标记感兴趣的基因

```r
# 方式1：手动指定
markGenes <- c("Oog4", "Psmd9", "Sephs2", "Lck", "Cd3d")

# 方式2：随机采样
markGenes <- sample(rownames(exps), 30)

# 方式3：基于差异倍数
# 假设有差异分析结果
de_genes <- read.csv("DE_results.csv")
markGenes <- de_genes %>% 
  filter(abs(log2FC) > 2, padj < 0.05) %>% 
  pull(gene)

# 添加到图中
pdf('marked_genes.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "heatmap",
  markGenes = markGenes,
  markGenesSide = "right",  # 或 "left"
  genesGp = c('italic', 10, 'darkred'),
  show_row_dend = FALSE
)
dev.off()
```

### 4.5 样本分组和排序

```r
# 添加分组
pdf('sample_groups.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "heatmap",
  sampleGroup = rep(c("Early", "Mid", "Late"), each = 2),
  sampleCol = c("Early" = "green", "Mid" = "orange", "Late" = "red")
)
dev.off()

# 改变样本顺序
visCluster(
  object = ck,
  plotType = "heatmap",
  sampleOrder = c("blastocyst", "tmorula", "t8.cell", 
                  "t4.cell", "t2.cell", "zygote")  # 反向时序
)

# 分组+排序+多注释
library(ComplexHeatmap)

group1 <- rep(c("A", "B"), each = 3)
names(group1) <- colnames(exps)

group2 <- rep(c("X", "Y", "Z"), each = 2)
names(group2) <- colnames(exps)

topAnno <- HeatmapAnnotation(
  Stage = group1,
  Type = group2,
  col = list(
    Stage = c("A" = "red", "B" = "blue"),
    Type = c("X" = "green", "Y" = "orange", "Z" = "purple")
  )
)

pdf('multi_anno.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "heatmap",
  heatmapAnnotation = topAnno,
  columnSplit = rep(c(1, 2, 3), each = 2)  # 分成3组
)
dev.off()
```

### 4.6 添加行注释

```r
# 准备行注释数据
df <- ck$wide.res %>% arrange(cluster)

# 创建一些注释
pathway_anno <- sample(c("PathwayA", "PathwayB", "PathwayC"), 
                       nrow(df), replace = TRUE)
names(pathway_anno) <- df$gene

expression_level <- sample(c("High", "Low"), 
                           nrow(df), replace = TRUE)
names(expression_level) <- df$gene

# 创建rowAnnotation对象
library(ComplexHeatmap)
leftAnno <- rowAnnotation(
  Pathway = pathway_anno,
  Level = expression_level,
  col = list(
    Pathway = c("PathwayA" = "red", "PathwayB" = "blue", "PathwayC" = "green"),
    Level = c("High" = "darkred", "Low" = "lightblue")
  )
)

# 添加到图中
pdf('with_row_anno.pdf', height = 10, width = 7)
visCluster(
  object = ck,
  plotType = "heatmap",
  rowAnnotationObj = leftAnno,
  column_names_rot = 45
)
dev.off()
```

### 4.7 趋势线注释

```r
# 基础趋势线
pdf('with_trend.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "right"
)
dev.off()

# 箱线图代替趋势线
pdf('with_boxplot.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "both",
  addBox = TRUE,
  addLine = FALSE,
  boxCol = ggsci::pal_npg()(8)
)
dev.off()

# 箱线图+点图
pdf('box_and_point.pdf', height = 10, width = 6)
visCluster(
  object = ck,
  plotType = "both",
  addBox = TRUE,
  addPoint = TRUE,
  addLine = FALSE,
  boxCol = ggsci::pal_npg()(8),
  pointArg = c(19, "orange", "orange", 1)
)
dev.off()
```

### 4.8 多组比较

```r
# 假设前3个样本是对照，后3个是处理
pdf('multi_condition.pdf', height = 10, width = 11)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",
  mulGroup = c(3, 3),
  mlineCol = c("blue", "red"),
  lgdLabel = c("Control", "Treatment"),
  annoTermData = enrich,
  annoTermMside = "right"
)
dev.off()
```

---

## 5. 高级应用场景

### 5.1 WGCNA结果可视化

WGCNA（Weighted Gene Co-expression Network Analysis）是识别共表达模块的强大工具。ClusterGVis可以很好地展示WGCNA结果。

#### 完整WGCNA流程
```r
library(WGCNA)
library(ClusterGVis)

# ===== 步骤1: 准备数据 =====
data(exps)
datExpr0 <- as.data.frame(t(exps))  # WGCNA需要转置（样本×基因）
rownames(datExpr0) <- colnames(exps)

# ===== 步骤2: 数据质控 =====
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# ===== 步骤3: 样本聚类，检查离群值 =====
sampleTree <- hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "")
abline(h = 15, col = "red")  # 设置阈值

clust <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]

# ===== 步骤4: 选择软阈值 =====
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 查看推荐的power
sft$powerEstimate

# ===== 步骤5: 一步法构建网络和识别模块 =====
net <- blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# 查看模块
table(net$colors)

# ===== 步骤6: 使用ClusterGVis可视化 =====
# 准备表达矩阵（要用WGCNA过滤后的）
exps_filtered <- t(datExpr)  # 转回基因×样本

# 聚类
ck <- clusterData(
  obj = exps_filtered,
  clusterMethod = "wgcna",
  object = net
)

# 富集分析
library(org.Mm.eg.db)
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",
  pvalueCutoff = 0.05,
  topn = 3
)

# 可视化
pdf('WGCNA_visualization.pdf', height = 12, width = 14)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",
  annoTermData = enrich,
  annoTermMside = "right",
  goCol = rep(ggsci::pal_d3()(length(unique(net$colors))), each = 3),
  column_names_rot = 45
)
dev.off()
```

**注意事项**：
1. `obj`参数必须是过滤后的表达矩阵（基因×样本）
2. `object`参数是`blockwiseModules()`的输出
3. 簇名称会包含WGCNA的模块颜色

#### 提取特定模块的基因
```r
# 提取模块1的基因
module1_genes <- ck$cluster.list$C1

# 导出所有模块
library(writexl)
write_xlsx(ck$cluster.list, "WGCNA_modules.xlsx")
```

### 5.2 单细胞数据可视化

ClusterGVis可以展示单细胞数据的marker基因和细胞类型特征。

#### 完整单细胞流程
```r
library(Seurat)
library(ClusterGVis)
library(SeuratData)

# ===== 步骤1: 准备Seurat对象 =====
data("pbmc_small")  # 示例数据
pbmc <- pbmc_small

# 或加载自己的数据
# pbmc <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# pbmc <- CreateSeuratObject(counts = pbmc, project = "pbmc3k", 
#                            min.cells = 3, min.features = 200)

# ===== 步骤2: 标准Seurat流程 =====
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# ===== 步骤3: 找marker基因 =====
# 所有marker
pbmc.markers.all <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 每个细胞类型的top10 marker
pbmc.markers <- pbmc.markers.all %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# ===== 步骤4: 准备ClusterGVis数据 =====
# 方式1: 显示平均表达
st.data.avg <- prepareDataFromscRNA(
  object = pbmc,
  diffData = pbmc.markers,
  showAverage = TRUE  # 每个细胞类型的平均值
)

# 方式2: 显示所有单细胞
st.data.all <- prepareDataFromscRNA(
  object = pbmc,
  diffData = pbmc.markers,
  showAverage = FALSE  # 所有细胞
)

# ===== 步骤5: 富集分析 =====
library(org.Hs.eg.db)
enrich <- enrichCluster(
  object = st.data.avg,
  OrgDb = org.Hs.eg.db,
  type = "BP",
  organism = "hsa",
  pvalueCutoff = 0.05,
  topn = 5
)

# ===== 步骤6: 基础可视化 =====
# 折线图
visCluster(
  object = st.data.avg,
  plotType = "line",
  ncol = 3
)

# 热图
pdf('scRNA_heatmap.pdf', height = 10, width = 8)
visCluster(
  object = st.data.avg,
  plotType = "heatmap",
  column_names_rot = 45,
  markGenes = head(pbmc.markers$gene, 20)
)
dev.off()

# 组合图+富集
pdf('scRNA_complete.pdf', height = 12, width = 14)
visCluster(
  object = st.data.avg,
  plotType = "both",
  lineSide = "left",
  annoTermData = enrich,
  annoTermMside = "right",
  column_names_rot = 45
)
dev.off()

# ===== 步骤7: 自定义排序 =====
# 按细胞类型重新排序
celltype_order <- c("CD8 T", "CD4 T", "NK", "B", "Mono", "DC")

pdf('scRNA_ordered.pdf', height = 10, width = 10)
visCluster(
  object = st.data.avg,
  plotType = "both",
  sampleOrder = celltype_order,
  column_names_rot = 45
)
dev.off()

# ===== 步骤8: 添加Feature Plot（第7章功能）=====
# 为每个cluster创建Feature Plot
gglist <- lapply(unique(pbmc.markers$cluster), function(x) {
  genes <- pbmc.markers %>% filter(cluster == x) %>% pull(gene)
  FeaturePlot(object = pbmc, features = genes, ncol = 4)
})
names(gglist) <- paste0("C", unique(pbmc.markers$cluster))

# 组合展示
pdf('scRNA_with_featurePlot.pdf', height = 12, width = 18)
visCluster(
  object = st.data.avg,
  plotType = "both",
  lineSide = "left",
  gglist = gglist,
  ggplotPanelArg = c(2, 0.25, 6, "grey90", NA),
  column_names_rot = 45
)
dev.off()
```

#### 显示所有单细胞（不平均）
```r
# 使用showAverage = FALSE
st.data.all <- prepareDataFromscRNA(
  object = pbmc,
  diffData = pbmc.markers,
  showAverage = FALSE
)

# 可视化时需要指定细胞类型顺序
pdf('scRNA_all_cells.pdf', height = 12, width = 15)
visCluster(
  object = st.data.all,
  plotType = "both",
  sampleCellOrder = levels(pbmc$seurat_clusters),  # 细胞类型顺序
  column_names_rot = 90,
  show_column_names = FALSE  # 太多细胞，不显示名称
)
dev.off()
```

### 5.3 拟时序数据可视化（Monocle）

ClusterGVis支持Monocle2和Monocle3的拟时序分析结果。

#### Monocle2流程
```r
library(monocle)
library(ClusterGVis)

# ===== 准备Monocle对象（这里省略标准流程）=====
# cds <- ...  # 假设已有Monocle2 CellDataSet对象

# ===== 使用ClusterGVis的Monocle2函数 =====
# 方式1: plot_pseudotime_heatmap2
pt_data <- plot_pseudotime_heatmap2(
  cds_subset = cds[marker_genes, ],
  cluster_rows = TRUE,
  num_clusters = 6,
  show_rownames = FALSE
)

# 方式2: plot_genes_branched_heatmap2（有分支）
branch_data <- plot_genes_branched_heatmap2(
  cds_subset = cds[marker_genes, ],
  branch_point = 1,
  num_clusters = 4
)

# ===== 可视化 =====
pdf('pseudotime.pdf', height = 10, width = 8)
visCluster(
  object = pt_data,
  plotType = "both",
  lineSide = "right"
)
dev.off()
```

#### Monocle3流程
```r
library(monocle3)
library(ClusterGVis)

# ===== 准备Monocle3对象 =====
# cds <- ...  # 假设已有Monocle3 cell_data_set对象

# ===== 直接使用clusterData =====
# ClusterGVis可以直接处理cell_data_set对象
pt_cluster <- clusterData(
  obj = cds,
  clusterMethod = "kmeans",
  clusterNum = 6
)

# 可视化
visCluster(
  object = pt_cluster,
  plotType = "both",
  pseudotimeCol = c("blue", "red")  # 拟时序颜色
)
```

### 5.4 批量富集分析和比较

有时需要对多个聚类结果进行富集分析并比较。

```r
library(ClusterGVis)
library(org.Mm.eg.db)

# ===== 不同聚类数的比较 =====
cluster_numbers <- c(4, 6, 8, 10)

results <- lapply(cluster_numbers, function(k) {
  # 聚类
  ck <- clusterData(
    obj = exps,
    clusterMethod = "kmeans",
    clusterNum = k
  )
  
  # 富集
  enrich <- enrichCluster(
    object = ck,
    OrgDb = org.Mm.eg.db,
    type = "BP",
    pvalueCutoff = 0.05,
    topn = 5
  )
  
  # 可视化
  pdf(paste0('cluster_k', k, '.pdf'), height = 10, width = 12)
  p <- visCluster(
    object = ck,
    plotType = "both",
    annoTermData = enrich
  )
  print(p)
  dev.off()
  
  list(cluster = ck, enrich = enrich)
})

names(results) <- paste0("k", cluster_numbers)

# ===== 比较不同k值的富集结果 =====
# 提取所有GO terms
all_terms <- lapply(results, function(x) {
  x$enrich %>% pull(Description) %>% unique()
})

# 找共同的GO terms
common_terms <- Reduce(intersect, all_terms)
print(paste("共同富集的GO terms:", length(common_terms)))
```

### 5.5 自定义基因集富集

```r
# ===== 准备自定义基因集 =====
# 例如：从MSigDB或其他来源
custom_genesets <- list(
  "MyPathway1" = c("Gene1", "Gene2", "Gene3", "Gene4"),
  "MyPathway2" = c("Gene5", "Gene6", "Gene7"),
  "MyPathway3" = c("Gene8", "Gene9", "Gene10", "Gene11", "Gene12")
)

# 转换为TERM2GENE格式
TERM2GENE <- data.frame(
  term = rep(names(custom_genesets), sapply(custom_genesets, length)),
  gene = unlist(custom_genesets)
)

TERM2NAME <- data.frame(
  term = names(custom_genesets),
  name = paste("Custom", names(custom_genesets))
)

# ===== 富集分析 =====
enrich <- enrichCluster(
  object = ck,
  type = "ownSet",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  idTrans = FALSE,
  pvalueCutoff = 0.05,
  topn = 3
)

# 可视化
pdf('custom_geneset.pdf', height = 10, width = 12)
visCluster(
  object = ck,
  plotType = "both",
  annoTermData = enrich,
  goCol = rep(c("red", "blue", "green"), each = 3)
)
dev.off()
```

### 5.6 整合多种注释

```r
library(ComplexHeatmap)
library(ClusterGVis)

# ===== 准备数据 =====
data(exps)
ck <- clusterData(obj = exps, clusterMethod = "kmeans", clusterNum = 8)

# GO富集
go_enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",
  topn = 3
)

# KEGG富集
kegg_enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "KEGG",
  organism = "mmu",
  topn = 2
)

# 准备样本注释
stage <- rep(c("Early", "Mid", "Late"), each = 2)
names(stage) <- colnames(exps)

treatment <- rep(c("Control", "Treatment"), 3)
names(treatment) <- colnames(exps)

topAnno <- HeatmapAnnotation(
  Stage = stage,
  Treatment = treatment,
  col = list(
    Stage = c("Early" = "green", "Mid" = "orange", "Late" = "red"),
    Treatment = c("Control" = "blue", "Treatment" = "purple")
  )
)

# 准备行注释
df <- ck$wide.res %>% arrange(cluster)
pathway <- sample(c("P1", "P2", "P3"), nrow(df), replace = TRUE)
names(pathway) <- df$gene

leftAnno <- rowAnnotation(
  Pathway = pathway,
  col = list(Pathway = c("P1" = "pink", "P2" = "cyan", "P3" = "yellow"))
)

# 标记基因
markGenes <- sample(rownames(exps), 30)

# ===== 整合所有注释 =====
pdf('integrated_annotations.pdf', height = 14, width = 18)
visCluster(
  object = ck,
  plotType = "both",
  
  # 热图设置
  htColList = list(
    col_range = c(-2, 0, 2),
    col_color = c("#08519C", "white", "#A50F15")
  ),
  
  # 列注释
  heatmapAnnotation = topAnno,
  columnSplit = rep(c(1, 2, 3), each = 2),
  column_names_rot = 45,
  
  # 行注释
  rowAnnotationObj = leftAnno,
  
  # 标记基因
  markGenes = markGenes,
  markGenesSide = "left",
  genesGp = c('italic', 8, 'darkred'),
  
  # 趋势线
  lineSide = "left",
  addBox = TRUE,
  boxCol = ggsci::pal_npg()(8),
  addLine = FALSE,
  
  # GO注释
  annoTermData = go_enrich,
  annoTermMside = "right",
  goCol = rep(ggsci::pal_d3()(8), each = 3),
  goSize = "pval",
  addBar = TRUE,
  barWidth = 6,
  
  # KEGG注释
  annoKeggData = kegg_enrich,
  annoKeggMside = "right",
  keggCol = rep(ggsci::pal_lancet()(8), each = 2),
  keggSize = 2.5,
  
  # 面板设置
  panelArg = c(2, 0.25, 4, "grey90", NA),
  ggplotPanelArg = c(2, 0.25, 5, "grey90", NA),
  
  # ComplexHeatmap参数
  show_row_dend = FALSE,
  border = TRUE
)
dev.off()
```

---

## 6. 源码关键设计

### 6.1 数据结构设计

#### clusterData返回对象
```r
# 返回一个list，包含5个元素
list(
  wide.res = data.frame(...),   # 宽格式：基因×样本+cluster
  long.res = data.frame(...),   # 长格式：用于ggplot2
  cluster.list = list(...),     # 基因列表，按簇分组
  type = "kmeans",              # 聚类方法
  geneMode = "none",            # 预留字段
  geneType = "none"             # 预留字段
)
```

**wide.res结构**：
- 适合数据操作和导出
- 每行一个基因，包含所有样本的表达值
- 额外列：cluster（簇ID）、membership（如果是mfuzz/TCseq）

**long.res结构**：
- 适合ggplot2绘图
- 每行一个基因-样本组合
- 列：cluster、gene、cell_type、norm_value、cluster_name、membership（如果有）

**cluster.list结构**：
- 便于提取特定簇的基因
- 键是簇ID（C1, C2, ...），值是基因向量

#### enrichCluster返回对象
```r
# 返回一个data.frame
data.frame(
  group = c("C1", "C1", "C2", ...),        # 簇ID
  Description = c("term1", "term2", ...),  # GO/KEGG描述
  pvalue = c(0.001, 0.005, ...),           # p值
  ratio = c(42.1, 36.8, ...),              # 富集比例
  geneID = c("Gene1/Gene2", ...)           # 基因列表（如果addGene=TRUE）
)
```

### 6.2 关键算法实现

#### getClusters - Elbow方法
```r
# 核心：计算不同k值的WSS（组内平方和）
factoextra::fviz_nbclust(exp, stats::kmeans, method = "wss")

# WSS计算原理：
# WSS = Σ(每个簇内所有点到簇中心的距离平方和)
# k越大，WSS越小
# 最优k在WSS下降速度明显放缓的点（肘部）
```

#### clusterData - 四种聚类算法

##### 1. mfuzz实现
```r
# 创建ExpressionSet对象
myset <- Biobase::ExpressionSet(assayData = as.matrix(exp))

# 过滤低变异基因
myset <- filter.std(myset, minStd = minStd, visu = FALSE)

# 标准化
if (scaleData) myset <- standardise(myset)

# 估计模糊参数m
m <- mestimate(myset)

# 执行模糊C均值聚类
mfuzz_res <- e1071::cmeans(
  Biobase::exprs(myset),
  centers = clusterNum,
  method = "cmeans",
  m = m
)

# mfuzz_res包含：
# - cluster: 每个基因的主要簇
# - membership: 隶属度矩阵（基因×簇）
```

**membership解释**：
- 值在0-1之间
- 对于基因i和簇j，membership[i,j]表示基因i属于簇j的程度
- Σ(membership[i,]) = 1（每个基因的membership和为1）

##### 2. kmeans实现
```r
# 数据标准化
if (scaleData) {
  hclust_matrix <- exp %>% t() %>% scale() %>% t()
} else {
  hclust_matrix <- exp
}

# 执行kmeans
km <- stats::kmeans(
  x = hclust_matrix,
  centers = clusterNum,
  nstart = 10  # 多次随机初始化，取最好结果
)

# km包含：
# - cluster: 每个基因的簇分配
# - centers: 每个簇的中心
# - withinss: 每个簇的WSS
# - tot.withinss: 总WSS
```

##### 3. WGCNA解析
```r
# 提取WGCNA结果
cinfo <- data.frame(
  cluster = net$colors + 1,  # WGCNA的colors是0-based
  modulecol = WGCNA::labels2colors(net$colors)
)

# 标准化表达矩阵
expm <- data.frame(t(scale(exp)))

# 合并信息
final.res <- cbind(expm, cinfo)
```

**WGCNA特点**：
- 不重新聚类，只是重新组织数据
- 保留模块颜色信息
- 簇名称格式："cluster 1 (200 turquoise)"

#### enrichCluster - 富集分析流程

```r
# 对每个簇循环
purrr::map_df(seq_len(length(unique(cluster_ids))), function(x) {
  # 1. 提取该簇的基因
  genes <- enrich.data %>% filter(cluster == x) %>% pull(gene)
  
  # 2. ID转换（如果需要）
  if (idTrans) {
    gene.ent <- clusterProfiler::bitr(
      genes,
      fromType = fromType,
      toType = toType,
      OrgDb = OrgDb
    )
    target_genes <- gene.ent[[toType]]
  } else {
    target_genes <- genes
  }
  
  # 3. 执行富集
  if (type %in% c("BP", "MF", "CC")) {
    ego <- clusterProfiler::enrichGO(
      gene = target_genes,
      keyType = toType,
      OrgDb = OrgDb,
      ont = type,
      pvalueCutoff = 1,  # 先不过滤
      readable = readable
    )
  } else if (type == "KEGG") {
    ego <- clusterProfiler::enrichKEGG(
      gene = target_genes,
      organism = organism,
      pvalueCutoff = 1
    )
  }
  
  # 4. 提取top N
  result <- ego@result %>%
    filter(pvalue < pvalueCutoff) %>%
    arrange(pvalue) %>%
    head(topn[x])
  
  # 5. 添加簇信息
  result$group <- paste0("C", x)
  
  return(result)
})
```

### 6.3 可视化架构

#### visCluster架构
```r
visCluster <- function(...) {
  # 1. 参数检查和预处理
  plotType <- match.arg(plotType)
  
  # 2. 根据plotType分支
  if (plotType == "line") {
    # 使用ggplot2绘制折线图
    ggplot(data, aes(x = cell_type, y = norm_value)) +
      geom_line(aes(group = gene)) +
      facet_wrap(~cluster_name)
      
  } else if (plotType == "heatmap") {
    # 使用ComplexHeatmap绘制热图
    Heatmap(
      matrix,
      col = col_fun,
      top_annotation = ...,
      left_annotation = ...,
      right_annotation = ...
    )
    
  } else if (plotType == "both") {
    # 组合：热图 + 趋势注释
    # 创建主热图
    ht <- Heatmap(matrix, ...)
    
    # 根据lineSide添加左侧或右侧注释
    if (lineSide == "left") {
      left_anno <- create_trend_annotation()
      ht <- left_anno + ht
    } else {
      right_anno <- create_trend_annotation()
      ht <- ht + right_anno
    }
    
    # 添加GO/KEGG注释
    if (!is.null(annoTermData)) {
      term_anno <- create_term_annotation()
      ht <- ht + term_anno
    }
    
    # 绘制
    draw(ht)
  }
}
```

#### 注释系统层次
```
主热图 (Heatmap)
├── 顶部注释 (top_annotation)
│   ├── 样本分组
│   └── 其他列注释
├── 左侧注释 (left_annotation)
│   ├── 行注释（当lineSide="right"时）
│   └── 趋势图（当lineSide="left"时）
└── 右侧注释 (right_annotation)
    ├── 簇注释条
    ├── 基因标签（当markGenesSide="right"时）
    ├── 趋势图（当lineSide="right"时）
    ├── GO注释
    ├── KEGG注释
    └── 自定义ggplot（gglist）
```

### 6.4 颜色系统

#### 热图颜色映射
```r
# 使用circlize包创建颜色映射函数
col_fun <- circlize::colorRamp2(
  breaks = col_range,   # 如 c(-2, 0, 2)
  colors = col_color    # 如 c("blue", "white", "red")
)

# 应用到热图
Heatmap(matrix, col = col_fun)
```

#### membership梯度颜色（mfuzz）
```r
# 将membership值映射到颜色
# membership从0到1，对应msCol的颜色梯度
ggplot(data, aes(x = cell_type, y = norm_value)) +
  geom_line(aes(color = membership, group = gene)) +
  scale_color_gradientn(colors = msCol)
```

### 6.5 数据流转换

#### 宽格式 → 长格式
```r
# 使用reshape2::melt
long_data <- reshape2::melt(
  wide_data,
  id.vars = c("cluster", "gene", "membership"),
  variable.name = "cell_type",
  value.name = "norm_value"
)

# 添加簇名称（包含基因数）
long_data$cluster_name <- paste0(
  "cluster ", long_data$cluster,
  " (", table(wide_data$cluster)[long_data$cluster], ")"
)
```

#### 长格式 → 矩阵（用于热图）
```r
# 使用dplyr + tidyr
matrix_data <- long_data %>%
  select(gene, cell_type, norm_value) %>%
  pivot_wider(names_from = cell_type, values_from = norm_value) %>%
  column_to_rownames("gene") %>%
  as.matrix()
```

### 6.6 性能优化技巧

#### 1. 大数据集处理
```r
# 使用光栅化
visCluster(
  object = ck,
  plotType = "heatmap",
  use_raster = TRUE,
  raster_quality = 2
)
```

#### 2. 减少基因数
```r
# 预先过滤
ck <- clusterData(
  obj = exps,
  minStd = 0.5,  # 提高阈值
  clusterNum = 8
)
```

#### 3. 批量处理
```r
# 使用lapply而不是for循环
results <- lapply(cluster_numbers, function(k) {
  clusterData(obj = exps, clusterNum = k)
})
```

---

## 7. 最佳实践

### 7.1 聚类数选择策略

#### 方法1：Elbow方法（推荐）
```r
# 使用getClusters
p <- getClusters(obj = exps)
print(p)

# 观察肘部位置
# WSS下降明显放缓的点就是最优k
```

#### 方法2：Silhouette系数
```r
library(cluster)
library(factoextra)

# 计算silhouette
fviz_nbclust(exps, kmeans, method = "silhouette")

# silhouette系数越高越好
# 范围：-1到1，>0.5表示结构良好
```

#### 方法3：Gap统计量
```r
library(cluster)

# 计算gap统计量
gap_stat <- clusGap(as.matrix(exps), 
                    FUN = kmeans, 
                    K.max = 15, 
                    B = 50)

fviz_gap_stat(gap_stat)

# 选择gap最大的k
```

#### 方法4：生物学意义
```r
# 尝试多个k值，比较生物学解释性
k_values <- c(4, 6, 8, 10)

for (k in k_values) {
  ck <- clusterData(obj = exps, clusterNum = k)
  enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, type = "BP")
  
  # 评估：
  # - 每个簇是否有清晰的生物学功能？
  # - 簇之间是否有明显差异？
  # - 是否符合预期的生物学过程？
}
```

### 7.2 数据预处理建议

#### 归一化选择
```r
# TPM/FPKM/RPKM已归一化，直接使用
exps_tpm <- read.csv("tpm_matrix.csv")

# count数据需要归一化
library(edgeR)
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge)
exps_cpm <- cpm(dge, log = FALSE)

# 或使用DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(count_matrix, colData, ~condition)
dds <- estimateSizeFactors(dds)
exps_norm <- counts(dds, normalized = TRUE)
```

#### 过滤低表达基因
```r
# 过滤低计数基因
keep <- rowSums(count_matrix >= 10) >= 3  # 至少3个样本 >= 10 counts
exps_filtered <- count_matrix[keep, ]

# 过滤低变异基因
gene_var <- apply(exps, 1, var)
keep <- gene_var > quantile(gene_var, 0.5)  # 保留前50%
exps_filtered <- exps[keep, ]
```

#### log转换
```r
# 对于count类数据，建议log转换
exps_log <- log2(exps + 1)

# 对于已经归一化的数据，视情况而定
# TPM/FPKM通常不需要log转换
```

### 7.3 可视化最佳实践

#### 使用pdf保存复杂图形
```r
# ❌ 错误：直接在RStudio查看
visCluster(object = ck, plotType = "both")
# 布局可能错位

# ✅ 正确：保存为PDF
pdf('result.pdf', height = 10, width = 12)
visCluster(object = ck, plotType = "both")
dev.off()
```

#### 图形尺寸选择
```r
# 基础热图
pdf('heatmap.pdf', height = 10, width = 6)

# 添加GO注释
pdf('with_GO.pdf', height = 10, width = 12)

# 添加自定义ggplot
pdf('with_custom.pdf', height = 12, width = 16)

# 大数据集
pdf('large_dataset.pdf', height = 15, width = 10)
```

#### 颜色方案选择
```r
# 学术期刊配色
library(ggsci)

# Nature系列
ctAnnoCol = pal_npg()(8)

# Science系列
ctAnnoCol = pal_aaas()(8)

# Lancet系列
ctAnnoCol = pal_lancet()(8)

# 自定义渐变
htColList = list(
  col_range = c(-2, 0, 2),
  col_color = c("#2166AC", "white", "#B2182B")  # RdBu配色
)
```

#### 字体大小调整
```r
visCluster(
  object = ck,
  plotType = "both",
  
  # 列名
  column_names_gp = gpar(fontsize = 10),
  
  # 基因标签
  genesGp = c('italic', 8, 'black'),
  
  # GO文本
  termTextLimit = c(8, 15),  # 最小和最大字号
  
  # 簇标注文本
  annnoblockGp = c("white", 8)
)
```

### 7.4 富集分析技巧

#### 多层次富集
```r
# 1. 先做BP
bp_enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, type = "BP", topn = 5)

# 2. 再做MF
mf_enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, type = "MF", topn = 3)

# 3. 最后做KEGG
kegg_enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, 
                             type = "KEGG", organism = "mmu", topn = 3)

# 合并展示
combined_enrich <- rbind(
  bp_enrich %>% mutate(Type = "BP"),
  mf_enrich %>% mutate(Type = "MF"),
  kegg_enrich %>% mutate(Type = "KEGG")
)
```

#### 控制富集term数量
```r
# 每个簇不同的topn
enrich <- enrichCluster(
  object = ck,
  OrgDb = org.Mm.eg.db,
  type = "BP",
  topn = c(5, 3, 7, 4, 6, 5, 4, 8)  # 8个簇，各自的topn
)
```

#### 富集结果后处理
```r
# 去除冗余term
library(rrvgo)

# 计算相似性矩阵
simMatrix <- calculateSimMatrix(
  enrich$Description,
  orgdb = "org.Mm.eg.db",
  ont = "BP",
  method = "Rel"
)

# 移除冗余
scores <- setNames(-log10(enrich$pvalue), enrich$Description)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7)

# 过滤enrich
enrich_reduced <- enrich %>%
  filter(Description %in% reducedTerms$parent)
```

### 7.5 数据导出

#### 导出聚类结果
```r
library(writexl)

# 方式1：导出宽格式
write_xlsx(ck$wide.res, "cluster_results_wide.xlsx")

# 方式2：按簇分sheet导出
write_xlsx(ck$cluster.list, "cluster_results_by_cluster.xlsx")

# 方式3：导出长格式
write_xlsx(ck$long.res, "cluster_results_long.xlsx")
```

#### 导出富集结果
```r
# 单个Excel文件
write_xlsx(enrich, "enrichment_results.xlsx")

# 按簇分sheet
enrich_list <- split(enrich, enrich$group)
write_xlsx(enrich_list, "enrichment_by_cluster.xlsx")
```

#### 导出图形
```r
# PDF（矢量图，推荐）
pdf('figure.pdf', height = 10, width = 12)
visCluster(...)
dev.off()

# PNG（位图，用于演示）
png('figure.png', height = 1000, width = 1200, res = 150)
visCluster(...)
dev.off()

# TIFF（高分辨率，用于投稿）
tiff('figure.tiff', height = 3000, width = 3600, res = 300, compression = "lzw")
visCluster(...)
dev.off()
```

### 7.6 报告生成

#### 使用RMarkdown生成完整报告
```r
# report.Rmd
---
title: "ClusterGVis Analysis Report"
author: "Your Name"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    code_folding: hide
---

## Data Loading
```{r}
library(ClusterGVis)
data(exps)
```

## Determine Optimal Clusters
```{r}
getClusters(obj = exps)
```

## Clustering
```{r}
ck <- clusterData(obj = exps, clusterMethod = "kmeans", clusterNum = 8)
```

## Enrichment
```{r}
library(org.Mm.eg.db)
enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, type = "BP")
knitr::kable(head(enrich, 20))
```

## Visualization
```{r fig.height=10, fig.width=12}
visCluster(object = ck, plotType = "both", annoTermData = enrich)
```
```

渲染报告：
```r
rmarkdown::render("report.Rmd")
```

---

## 8. 常见问题与技巧

### 8.1 常见错误及解决

#### 问题1：图形布局错位
```r
# ❌ 问题
visCluster(object = ck, plotType = "both")
# 在RStudio中查看，趋势线和热图对不齐

# ✅ 解决
pdf('correct.pdf', height = 10, width = 6)
visCluster(object = ck, plotType = "both")
dev.off()
# 用PDF查看器打开，布局正确
```

#### 问题2：基因标签顺序不匹配
```r
# ❌ 问题
markGenes <- c("Gene1", "Gene2", "Gene3")
visCluster(object = ck, markGenes = markGenes)
# 标签位置可能不对

# ✅ 解决
# 确保基因顺序与聚类后的顺序一致
df <- ck$wide.res %>% arrange(cluster)
markGenes <- df$gene[c(1, 50, 100)]  # 使用聚类后的基因名
visCluster(object = ck, markGenes = markGenes)
```

#### 问题3：ID转换失败
```r
# ❌ 问题
enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, 
                       fromType = "ENSEMBL", toType = "ENTREZID")
# Error: 很多ID无法转换

# ✅ 解决方案1：检查ID格式
head(rownames(ck$wide.res))
# 如果已经是SYMBOL，使用fromType = "SYMBOL"

# ✅ 解决方案2：预先转换
library(clusterProfiler)
id_map <- bitr(rownames(ck$wide.res), 
              fromType = "ENSEMBL",
              toType = "ENTREZID",
              OrgDb = org.Mm.eg.db)
# 检查转换率
nrow(id_map) / nrow(ck$wide.res)
```

#### 问题4：内存不足
```r
# ❌ 问题
# 大数据集（>10000基因）导致内存溢出

# ✅ 解决
# 方案1：过滤基因
ck <- clusterData(obj = exps, minStd = 0.5)  # 提高阈值

# 方案2：使用光栅化
visCluster(ck, plotType = "heatmap", use_raster = TRUE)

# 方案3：分批处理
clusters <- split(ck$wide.res$gene, ck$wide.res$cluster)
for (i in seq_along(clusters)) {
  subset_data <- exps[clusters[[i]], ]
  # 处理每个簇
}
```

#### 问题5：富集无结果
```r
# ❌ 问题
enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, 
                       type = "BP", pvalueCutoff = 0.05)
# 返回空或很少结果

# ✅ 解决
# 方案1：放宽p值阈值
enrich <- enrichCluster(ck, pvalueCutoff = 0.1, topn = 10)

# 方案2：检查基因数
table(ck$wide.res$cluster)
# 如果某些簇基因太少（<10），考虑合并簇

# 方案3：换富集类型
enrich_mf <- enrichCluster(ck, type = "MF")  # 试试MF
enrich_kegg <- enrichCluster(ck, type = "KEGG", organism = "mmu")
```

### 8.2 性能优化技巧

#### 技巧1：预先过滤
```r
# 在聚类前过滤
# 1. 去除低表达
keep_expr <- rowMeans(exps) > 1
exps_filtered <- exps[keep_expr, ]

# 2. 去除低变异
gene_var <- apply(exps_filtered, 1, var)
keep_var <- gene_var > quantile(gene_var, 0.25)  # 保留前75%
exps_final <- exps_filtered[keep_var, ]

# 现在聚类会快很多
ck <- clusterData(obj = exps_final, clusterNum = 8)
```

#### 技巧2：并行处理
```r
library(parallel)
library(foreach)
library(doParallel)

# 注册并行后端
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# 并行处理多个k值
results <- foreach(k = 4:10, .packages = "ClusterGVis") %dopar% {
  ck <- clusterData(obj = exps, clusterNum = k)
  enrich <- enrichCluster(ck, OrgDb = org.Mm.eg.db, type = "BP")
  list(cluster = ck, enrich = enrich)
}

stopCluster(cl)
```

#### 技巧3：缓存中间结果
```r
# 使用缓存避免重复计算
cache_file <- "cluster_k8.rds"

if (file.exists(cache_file)) {
  ck <- readRDS(cache_file)
} else {
  ck <- clusterData(obj = exps, clusterNum = 8)
  saveRDS(ck, cache_file)
}
```

### 8.3 高级技巧

#### 技巧1：自定义聚类
```r
# 如果有自己的聚类结果，可以手动构造ClusterGVis对象
my_clusters <- kmeans(t(scale(t(exps))), centers = 8)

# 构造wide.res
wide.res <- data.frame(
  t(scale(t(exps))),
  gene = rownames(exps),
  cluster = my_clusters$cluster,
  check.names = FALSE
) %>% arrange(cluster)

# 构造long.res
long.res <- wide.res %>%
  pivot_longer(
    cols = -c(gene, cluster),
    names_to = "cell_type",
    values_to = "norm_value"
  ) %>%
  mutate(cluster_name = paste0("cluster ", cluster, 
                               " (", table(wide.res$cluster)[cluster], ")"))

# 构造cluster.list
cluster.list <- split(wide.res$gene, paste0("C", wide.res$cluster))

# 创建对象
ck <- list(
  wide.res = wide.res,
  long.res = long.res,
  cluster.list = cluster.list,
  type = "custom",
  geneMode = "none",
  geneType = "none"
)

# 现在可以使用ClusterGVis的可视化功能
visCluster(object = ck, plotType = "both")
```

#### 技巧2：提取ComplexHeatmap对象
```r
# 获取热图对象而不直接绘制
ht <- visCluster(object = ck, plotType = "heatmap")

# 现在可以进一步操作
# 添加额外的注释
library(ComplexHeatmap)

extra_anno <- rowAnnotation(
  foo = anno_text(sample(LETTERS, nrow(ck$wide.res), replace = TRUE))
)

# 重新绘制
pdf('modified.pdf', height = 10, width = 8)
draw(ht + extra_anno)
dev.off()
```

#### 技巧3：批量生成图形
```r
# 为每个簇生成单独的详细图
library(dplyr)

for (i in 1:8) {
  # 提取该簇的数据
  cluster_genes <- ck$cluster.list[[paste0("C", i)]]
  cluster_data <- exps[cluster_genes, ]
  
  # 重新聚类（更细分）
  sub_ck <- clusterData(
    obj = cluster_data,
    clusterMethod = "kmeans",
    clusterNum = 3
  )
  
  # 富集
  sub_enrich <- enrichCluster(
    object = sub_ck,
    OrgDb = org.Mm.eg.db,
    type = "BP",
    topn = 5
  )
  
  # 可视化
  pdf(paste0('cluster_', i, '_detail.pdf'), height = 8, width = 10)
  visCluster(
    object = sub_ck,
    plotType = "both",
    annoTermData = sub_enrich
  )
  dev.off()
}
```

#### 技巧4：交互式可视化
```r
library(plotly)

# 将折线图转为交互式
p <- visCluster(object = ck, plotType = "line")
ggplotly(p)

# 或使用heatmaply
library(heatmaply)
heatmaply(
  ck$wide.res %>% select(-gene, -cluster) %>% as.matrix(),
  Rowv = FALSE,
  Colv = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(100)
)
```

### 8.4 调试技巧

#### 检查数据结构
```r
# 检查聚类结果
str(ck)
summary(ck$wide.res)
table(ck$wide.res$cluster)

# 检查富集结果
head(enrich)
table(enrich$group)
summary(enrich$pvalue)
```

#### 逐步可视化
```r
# 从简单开始
visCluster(object = ck, plotType = "line")  # 先看折线

visCluster(object = ck, plotType = "heatmap")  # 再看热图

# 逐步添加元素
pdf('step1.pdf', height = 10, width = 6)
visCluster(object = ck, plotType = "both")  # 基础组合
dev.off()

pdf('step2.pdf', height = 10, width = 8)
visCluster(
  object = ck, 
  plotType = "both",
  markGenes = head(rownames(exps), 20)  # 添加基因标签
)
dev.off()

pdf('step3.pdf', height = 10, width = 12)
visCluster(
  object = ck,
  plotType = "both",
  markGenes = head(rownames(exps), 20),
  annoTermData = enrich  # 添加GO注释
)
dev.off()
```

#### 使用traceback
```r
# 如果出错
options(error = recover)  # 进入调试模式

# 或使用browser()
visCluster <- function(...) {
  browser()  # 在这里暂停
  # ... 原始代码
}
```

### 8.5 文档和帮助

#### 查看函数文档
```r
?getClusters
?clusterData
?enrichCluster
?visCluster
```

#### 查看示例
```r
example(visCluster)
```

#### 在线资源
- 官方文档：https://github.com/junjunlab/ClusterGVis
- 教程：https://junjunlab.github.io/ClusterGVis/
- 问题反馈：https://github.com/junjunlab/ClusterGVis/issues

---

## 附录：快速参考

### A. 核心函数参数速查

#### getClusters()
```r
getClusters(
  obj = NULL,  # 表达矩阵/CellDataSet/SummarizedExperiment
  ...          # 额外参数
)
```

#### clusterData()
```r
clusterData(
  obj = NULL,
  scaleData = TRUE,
  clusterMethod = c("mfuzz", "TCseq", "kmeans", "wgcna"),
  clusterNum = NULL,
  minStd = 0,
  subcluster = NULL
)
```

#### enrichCluster()
```r
enrichCluster(
  object = NULL,
  type = c("BP", "MF", "CC", "KEGG", "ownSet"),
  OrgDb = NULL,
  idTrans = TRUE,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  organism = "hsa",
  pvalueCutoff = 0.05,
  topn = 5
)
```

#### visCluster() - 核心参数
```r
visCluster(
  object = NULL,
  plotType = c("line", "heatmap", "both"),
  
  # 热图颜色
  htColList = list(col_range = c(-2, 0, 2),
                   col_color = c("#08519C", "white", "#A50F15")),
  
  # 折线图
  lineSize = 0.1,
  lineCol = "grey90",
  msCol = c("#0099CC", "grey90", "#CC3333"),  # mfuzz
  
  # 趋势注释
  lineSide = "right",
  addBox = FALSE,
  addLine = TRUE,
  
  # 标记基因
  markGenes = NULL,
  markGenesSide = "right",
  
  # GO注释
  annoTermData = NULL,
  annoTermMside = "right",
  goCol = NULL,
  goSize = NULL,
  
  # 样本注释
  sampleGroup = NULL,
  sampleOrder = NULL,
  heatmapAnnotation = NULL,
  
  # 行注释
  rowAnnotationObj = NULL,
  
  # 自定义ggplot
  gglist = NULL,
  
  # ComplexHeatmap参数
  column_names_rot = 45,
  show_row_dend = FALSE,
  ...
)
```

### B. 常用配色方案

```r
library(ggsci)

# Nature系列
pal_npg()(10)

# Science系列
pal_aaas()(10)

# Lancet系列
pal_lancet()(10)

# JCO系列
pal_jco()(10)

# 自定义渐变
colorRampPalette(c("blue", "white", "red"))(100)
```

### C. 数据格式示例

#### 输入表达矩阵
```
         sample1  sample2  sample3  sample4
Gene1     1.313    1.237    1.326    1.262
Gene2     1.092    1.316    1.174    1.065
Gene3     0.986    1.201    1.123    1.085
```

#### clusterData输出（wide.res）
```
   gene    sample1  sample2  sample3  cluster  membership
1  Gene1    0.255    0.608    0.708      1        0.765
2  Gene2   -1.193    0.347    0.113      1        0.437
```

#### enrichCluster输出
```
  group              Description       pvalue    ratio
1    C1  lymphocyte differentiation  4.26e-09  42.10526
2    C1       T cell differentiation  1.28e-08  36.84211
```

---

## 结语

ClusterGVis是一个功能强大、高度灵活的R包，特别适合基因表达聚类分析和可视化。通过本文档，你应该能够：

1. ✅ 理解ClusterGVis的核心工作流程
2. ✅ 掌握四个核心函数的使用
3. ✅ 应用于bulk RNA-seq、单细胞、WGCNA等多种场景
4. ✅ 创建出版级别的可视化图表
5. ✅ 解决常见问题并优化性能

**学习建议**：
- 从基础教程（第4章）开始实践
- 尝试不同的聚类方法和参数
- 逐步添加注释元素
- 参考高级应用场景（第5章）
- 遇到问题查阅第8章

**进一步学习资源**：
- 官方GitHub：https://github.com/junjunlab/ClusterGVis
- ComplexHeatmap文档：https://jokergoo.github.io/ComplexHeatmap-reference/book/
- clusterProfiler教程：https://yulab-smu.top/biomedical-knowledge-mining-book/

祝你分析顺利！🎉
