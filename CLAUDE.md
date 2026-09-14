# CLAUDE.md — 项目导航与工作约定

## 项目背景

OLP（口腔扁平苔藓）微生物组研究。实验涉及流式细胞术、16S rRNA 扩增子、宏基因组、qPCR、FISH、IHC、裂解酶活性、细菌生长曲线等。测序原始数据在服务器上分析，本地只放仪器导出的数据（Excel/CSV）和服务器下载的处理后结果。

## 目录结构

```
分析/
├── R/                              # 所有可复用 R 代码（核心）
│   ├── utils/                      # 通用工具函数（my_themes, compare, filter_abundance 等）
│   ├── r_functions/                # 正式函数库 dev → lib → archive 开发流程
│   │   ├── dev/                    # 开发中的函数
│   │   ├── lib/                    # 稳定函数（项目脚本的 source 目标）
│   │   └── archive/                # 旧版本备份
│   ├── scratch/                    # 一次性/实验性脚本
│   └── R_FUNCTION_DEV_WORKFLOW.md  # 函数开发流程详细文档
├── projects/                       # 14个实验/分析项目，每个含 data/ + scripts/ + output/
│   ├── flow_cytometry/             # 流式细胞术
│   ├── 16S/                        # 16S rRNA
│   ├── metagenome/                 # 宏基因组
│   ├── lefse/                      # LEfSe 差异分析
│   ├── enrichment_replot/          # 富集重画
│   ├── fish/ ihc/ qpcr/            # FISH / 免疫组化 / qPCR
│   ├── lysin_enzyme/               # 裂解酶实验
│   ├── growth_curves/              # 生长曲线
│   ├── clinical_tables/            # 临床 Tableone
│   ├── phage_host/                 # 噬菌体-宿主互作
│   ├── venn_diagram/               # Venn 图
│   └── ACS_ADS_Global/             # ACS/ADS 全局分析
├── collaborators/                  # 帮其他人做的分析（LHY/WY/YZX）
├── data/metadata/                  # 跨项目共用的样本信息
├── archive/                        # 废弃/旧文件（main.r 等）
├── python/                         # Python 脚本预留
├── ref/                            # 参考文献
└── output/                         # 跨项目输出
```

## 命名约定

- **文件夹**：英文 snake_case 小写，如 `flow_cytometry/`
- **R 脚本**：snake_case，`.R` 扩展名，如 `compare_plot.R`
- **数据文件**：可保留中文原名，如 `裂解酶0416.xlsx`
- **归档版本**：加日期后缀，如 `function_v1_20250405.R`

## 日常工作模式

### 模式 1：数据分析（最常见）
用户拿着仪器导出的 Excel/CSV 做分析和画图。流程：
1. 数据放 `projects/<项目>/data/`
2. 脚本放 `projects/<项目>/scripts/`
3. 输出放 `projects/<项目>/output/`
4. 在脚本中引用库函数：
```r
source("../../R/utils/my_themes.R")
source("../../R/r_functions/lib/compare_plot/compare_plot_optimized.R")
```

### 模式 2：帮合作者分析
在 `collaborators/<姓名>/` 下独立管理，每个人有独立的 data/scripts/figures。

### 模式 3：R 函数开发
遵循 `R_FUNCTION_DEV_WORKFLOW.md` 中的流程：
1. 在 `R/r_functions/dev/<function_name>/` 创建 function.R + test.R + notes.md
2. 开发测试通过后，归档到 `R/r_functions/lib/<category>/`
3. 旧版本备份到 `R/r_functions/archive/`

## 技术栈

- **主语言**：R（ggplot2 画图，tidyverse 数据处理）
- **辅助**：Python（偶尔，在 `python/` 下）
- **版本控制**：Git，已关联 GitHub remote origin
- **AI 工具**：Claude Code（`.claude/settings.local.json`）

## 跨项目共享知识库

`~/projects/research-commons`（GitHub: `Yelin112/research-commons`）是跨项目共享的知识库，收纳函数库、流程模板、方法笔记、工具手册。做分析前先去那里查有没有现成的，避免重复摸索；本项目里定型的可复用资产（函数、流程参数组合、方法论笔记）也应该按 `anthropic-skills:commons-harvest` skill 的规范沉淀过去，而不是只留在本项目里。

与本项目直接相关的资料：
- `references/biobakery/bioBakery-overview.md` —— bioBakery工具全景手册（HUMAnN/MetaPhlAn/StrainPhlAn/MaAsLin等），含 OLP口腔微生物组相关性速查、以及 Assembly workflow(sgb_pipeline) vs 手工 MEGAHIT+Prodigal+CD-HIT+eggNOG-mapper 流程的取舍结论。做宏基因组物种/功能谱分析、菌株分析、差异丰度检验前应先查这份手册。

云端 Claude Code 会话默认看不到 research-commons，需要显式 attach（`add_repo` owner=Yelin112 repo=research-commons）；本地/SSH 会话直接用文件路径访问。

## .gitignore 规则

脚本纳入 Git；PDF/PNG/SVG 输出不入库；data/metadata 不入库；RData/Rhistory 不入库。

## 注意事项

- 不要修改文件内容，只动目录结构（约定）
- 新分析优先在 `projects/` 下建独立项目，不要往根目录扔脚本
- 复用已有库函数，避免重复造轮子
- 函数开发参考 `R_FUNCTION_DEV_WORKFLOW.md` 中的模板
