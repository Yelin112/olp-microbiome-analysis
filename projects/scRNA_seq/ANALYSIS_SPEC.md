# OLP scRNA-seq 免疫/T 细胞亚群分析规范

> 本文档面向执行本分析的 Claude Code 实例（服务器端）。执行前请先完整阅读，
> 并对照仓库根目录 `CLAUDE.md` 的项目约定（路径引导、目录结构、配色/主题系统、
> 分工原则）。有任何与本文档冲突或本文档未覆盖的情况，以 `CLAUDE.md` 为准。

## 1. 背景与目标

OLP（口腔扁平苔藓）vs 健康对照的 scRNA-seq 数据，目标是回答：

1. OLP 状态下免疫细胞（尤其 T 细胞）的组成是否发生变化；
2. 哪些 T 细胞亚群的比例变化最明显；
3. 这些亚群是否存在特定信号通路的激活，激活程度在不同亚群间是否有差异。

最终交付：细胞类型/亚群注释结果、组间比例差异、亚群水平差异表达与通路富集结果，
以及支撑结论的图表（中间产出，非最终期刊图）。

## 2. 数据

- 服务器路径：`~/OLP_scRNA_seq/`（21 个样本子目录，每个目录是 10x Cell Ranger
  `filtered_feature_bc_matrix` 标准输出：`barcodes.tsv.gz` / `features.tsv.gz` /
  `matrix.mtx.gz`）。
- 样本清单：
  - Healthy: `Healthy1`–`Healthy5`（n=5）
  - OLP: `OLP1`–`OLP16`（n=16）
- **已确认**：仅有基因表达（GEX）矩阵，无 TCR/VDJ 数据（已向用户核实，
  不会补传）。本次分析**不做克隆型/克隆扩增分析**，通路激活的结论完全基于
  转录组特征（差异表达 + 通路打分）推断，需在最终报告的局限性中注明这一点。
- 原始数据不入库（遵循 `CLAUDE.md` 的 data 不入 git 约定），中间产出（RDS/CSV）
  也不入库，只有脚本入库。

## 3. 产出约定（不限制具体目录结构）

服务器上的目录组织由执行者（服务器端 Claude Code）自行决定，不强制对齐本地
`projects/<name>/data|scripts|output` 的三段式——本文档只约束内容规范：

- 脚本开头统一用项目标准引导（不要写死路径、不要用 `sys.frame(1)$ofile`）：
  ```r
  source(file.path(Sys.getenv("OLP_ROOT"), "R/init.R"))
  load_utils(theme = TRUE)
  ```
  运行前确认 `~/.Renviron` 中 `OLP_ROOT` 已设置，且跑过
  `Rscript R/check_deps.R` 确认依赖齐全（见第 4 节）。
- 每个阶段结束后落一份可重新加载的中间结果（.rds 或等价形式），支持断点续跑，
  具体存放路径由执行者自定，但要在脚本/日志里说明清楚。
- **分工原则（沿用 CLAUDE.md）**：服务器只做计算和导出中间表/RDS/PDF/PNG；
  不要在服务器上做 PPTX 导出或最终期刊级图表微调，那部分留到本地做。
- 图表文字一律英文（沿用 `R/PLOTTING_CONVENTIONS.md`），配色/主题必须走
  `get_colors()` / `scale_fill_pub_d()` / `theme_pub_*()`，**不要在脚本里
  自定义调色板或另起一套主题**。

## 4. 环境依赖

预计新增以下 R 包依赖，执行前检查是否已安装，缺失的登记进
`R/check_deps.R` 的 `DEPS` 表（只报告缺失，不自动装，装包需征得用户同意）：

- 核心：`Seurat` (>=5.0)、`SeuratObject`、`Matrix`
- 整合：`harmony`
- 双细胞检测：`scDblFinder`（Bioconductor）
- QC 离群值判定：不需要额外依赖，`median()`/`mad()`（base R）即可实现，
  想省事也可以用 `scuttle::isOutlier()`（可选）
- 差异表达/伪批量：`DESeq2`、`edgeR`
- 通路富集：`clusterProfiler`、`msigdbr`、`fgsea`
- 单细胞打分：`UCell` 或 Seurat 自带 `AddModuleScore`
- 组成差异检验：非必需，默认用按病人聚合比例 + t-test/ANOVA（见第 6 节）；
  `speckle`（`propeller`）/`miloR` 仅作可选加强项
- 细胞通讯：`CellChat`（第 10 阶段，纳入核心流程）
- 加速：`presto`（加速 Seurat 的 Wilcoxon 检验，仅用于探索性排序，不用于最终统计结论）

缺 `officer`/`rvg` 属于预期内（PPTX 本就不在服务器做），不影响本分析。

## 5. 分析流程（分阶段，每阶段结束落一个 checkpoint .rds）

### 阶段 0：数据核查
- 确认 21 个样本矩阵可正常读入，行数/列数符合预期。
- 核查是否存在 VDJ 数据（见第 2 节）。
- 输出：`output/checkpoints/00_raw_list.rds`（各样本原始 Seurat 对象列表）。

### 阶段 1：逐样本 QC
- 计算 UMI 数、基因数、线粒体基因比例（`^MT-`，注意大小写，人类数据）、
  红细胞/血红蛋白基因比例（口腔活检组织，排查血液污染）。
- 阈值用 MAD 自适应，逐样本单独定，不用全局硬阈值；实现上很轻量，
  对每个 QC 指标算 `median(x) ± k * mad(x)`（k 常用 3），几行 base R
  就能做完，Seurat 对象的 QC 列（`nCount_RNA`/`nFeature_RNA`/`percent.mt`等）
  本身就是普通数值向量，没有额外集成障碍。
- 记录每个样本过滤前后细胞数，检查 Healthy vs OLP 组间是否有系统性差异
  （写入 `output/tables/01_qc_summary.csv`）。
- 输出：`output/checkpoints/01_qc_filtered_list.rds`。

### 阶段 2：双细胞检测
- 逐样本单独跑（不要合并后跑），标记但先不急着删除，记录比例。
- 输出：`output/checkpoints/02_doublet_annotated_list.rds`。

### 阶段 3：整合（Harmony）
- 合并 21 个样本，按样本 ID 做 Harmony 整合（仅整合技术性批次，不整合
  Healthy/OLP 分组本身）。
- 输出整合前后 UMAP（按样本着色 vs 按 marker 基因着色）到
  `output/figures/03_integration_check.pdf`，用于人工判断整合是否合理。
- 输出：`output/checkpoints/03_integrated.rds`。

### 阶段 4：全局聚类与细胞大类注释
- 标准 Seurat 流程聚类，用经典 marker 划分：上皮、免疫（PTPRC+）、
  成纤维、内皮、肌细胞等大类。
- 纯 marker 基因人工注释，**不使用 label transfer / 参考图谱**（评估过
  CELLxGENE 的 Human Oral and Craniofacial Cell Atlas，组织类型虽然匹配，
  但最终决定注释工作由用户人工复核完成，不引入自动标签）。
- **人工复核断点**：本阶段结束后**暂停，不自动往下跑**。落好 checkpoint
  （细胞大类聚类结果 + UMAP + marker dotplot/热图 + 各 cluster 的
  marker 基因排名表）后，交给用户人工核对/修正大类标签，确认无误后再
  手动触发第 5 阶段。执行者不要自作主张替用户拍板细胞类型。
- 输出：`output/checkpoints/04_annotated_global.rds`，
  UMAP + marker dotplot 存入 `output/figures/`。

### 阶段 5：免疫细胞大类组成比较
- 抽取 CD45+ 免疫细胞，划分 T/B/髓系/NK/mast 等大类。
- 按病人计算各大类比例（**统计单元是病人，不是细胞**），默认用
  t-test/ANOVA 比较 OLP vs Healthy（与项目常规统计惯例一致）；
  `propeller`/`miloR` 作为可选加强项，非必需。
- 输出：`output/tables/05_immune_composition.csv` + 对应箱线图。

### 阶段 6：T 细胞二次聚类与亚群注释
- 单独抽取 T 细胞重新降维聚类（更高分辨率）。
- 亚群目标清单（至少覆盖）：CD4 naive / Th1 / Th17 / Tfh / Treg(FOXP3+)；
  CD8 naive / 效应细胞毒性(GZMB/GZMK/PRF1) / 耗竭样(PDCD1/HAVCR2/LAG3/TOX) /
  组织驻留记忆(CD69/ITGAE)；γδT；MAIT。
- 纯 marker 基因人工判定亚群身份，**不使用 label transfer**（同阶段 4）；
  命名尽量参照已发表口腔黏膜/皮肤扁平苔藓单细胞图谱交叉核对。
- **人工复核断点**：同阶段 4，本阶段结束后暂停，落好 checkpoint（T 细胞
  二次聚类结果 + UMAP + marker dotplot/热图 + 各 cluster marker 排名表）
  交给用户人工核对/修正亚群标签，确认无误后再手动触发第 7 阶段。
- 输出：`output/checkpoints/06_tcell_subclustered.rds`，
  UMAP + marker dotplot/热图。

### 阶段 7：T 细胞亚群比例组间差异
- 同阶段 5 的统计方法（默认 t-test/ANOVA），按病人计算各 T 细胞亚群
  占总 T 细胞比例，比较 OLP vs Healthy，做多重检验校正。
- 输出：`output/tables/07_tcell_subset_composition.csv` + 图。
- **这一步的结果直接回答"哪些 T 细胞亚群变化最明显"。**

### 阶段 8：亚群水平差异表达（pseudobulk）
- 对每个 T 细胞亚群，按"病人 × 亚群"聚合成伪批量样本，用
  DESeq2/edgeR 做 OLP vs Healthy 差异表达。**不要用 Seurat 默认的
  细胞级 Wilcoxon 检验作为最终统计结论**（伪重复问题，只能用于
  探索性排序）。
- 输出：每个亚群一份 DE 结果表，存 `output/tables/08_de_<subset>.csv`。

### 阶段 9：通路富集与单细胞通路打分
- 用阶段 8 的 DE 结果做 GSEA/超几何富集（Hallmark/KEGG/Reactome/GO），
  重点关注：IFN-γ/IFN-α 应答、TNF 信号、IL-17 信号、JAK-STAT、
  TCR 信号、细胞毒性/颗粒酶通路、耗竭特征基因集。
- 用 `AddModuleScore`/`UCell` 给每个细胞打分（细胞毒性、耗竭、Th17、
  IFN 应答等特征），按亚群、按组比较分布，UMAP 展示打分空间分布。
- 输出：`output/tables/09_pathway_enrichment_<subset>.csv` +
  `output/figures/09_module_scores.pdf`。
- **两层结果（pseudobulk 富集 + 单细胞打分）互相印证，是回答
  "是否存在特定通路激活、哪个亚群更明显"的核心证据。**

### 阶段 10：细胞间通讯
- CellChat 比较 OLP vs Healthy 的配体-受体信号强度，重点看变化明显的
  T 细胞亚群与上皮细胞/髓系细胞之间的通讯变化，作为"T 细胞攻击上皮"
  这一 OLP 核心病理机制的补充证据。
- 纳入核心流程，跑完阶段 9 后直接执行。

### 阶段 11（可选）：拟时序
- 若耗竭/效应分化相关亚群变化明显，用 Monocle3/Slingshot 做补充验证。
- 非必需，视前序结果决定。

## 6. 统计规范（贯穿全流程，必须遵守）

1. **病人是统计单元，细胞不是**——所有组间比例比较、差异表达都要
   在病人层面聚合或用能处理重复测量/伪重复的方法，不能直接把细胞当
   独立样本做检验。
2. **多重检验校正**——T 细胞亚群比较、通路富集都涉及多重假设检验，
   必须做校正（如 BH/FDR）。**注意**：项目里"vs-Control 默认
   t-test/ANOVA、不做多重校正"的惯例是针对常规少量组间图表比较，
   *不适用于*本分析里的大规模亚群/通路检验场景，这里必须单独校正。
3. **QC/整合阶段的系统性偏差**要在阶段 1、3 提前排查，避免技术性
   伪影被误读为生物学差异。

## 7. 执行方式提示

- 全程非交互式执行（`Rscript`），每阶段结束保存 checkpoint，支持从任意
  阶段断点续跑，不要一次性跑完不落中间结果。
- **阶段 4 和阶段 6 结束后必须真正停下来，等用户人工复核完注释结果再继续**
  ——不要自动接着跑后续阶段，也不要自己代替用户判断细胞类型是否正确。
- 每个阶段脚本单独成文件，按 `01_qc.R`、`02_doublet.R` ... 顺序命名，
  存入 `projects/scRNA_seq/scripts/`。
- 计算耗时较长的阶段（整合、pseudobulk DE、GSEA）建议 `nohup` /
  后台方式运行并记录日志。
- 脚本入库前对照 `CLAUDE.md` 里"新函数"的约定（如果过程中沉淀出可复用
  函数，按 `R_FUNCTION_DEV_WORKFLOW.md` 流程走 dev → lib，不要把可复用
  逻辑散落在一次性脚本里）。

## 8. 待确认事项（开始前需要用户拍板）

- [x] VDJ/TCR 数据：确认不存在，不做克隆型分析（见第 2 节）。
- [x] label transfer：确定不使用，纯 marker 基因人工注释；阶段 4（全局
      大类）和阶段 6（T 细胞亚群）结束后各设一个**人工复核断点**，
      由用户确认后再继续（见第 4、6、7 节）。
- [ ] 阶段 11（拟时序）：视阶段 6–9 结果是否支持相关结论再决定是否执行；
      阶段 10（细胞通讯）已纳入核心流程，不再单独确认。
