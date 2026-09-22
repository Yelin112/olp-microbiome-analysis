# OLP Microbiome 分析工作目录

OLP (口腔扁平苔藓) 微生物组研究相关的数据分析工作目录。

## 目录结构

```
分析/
├── R/                     # 所有可复用的 R 代码
│   ├── utils/             # 通用工具函数 (主题、丰度处理、比较分析等)
│   ├── r_functions/       # 正式函数库 (dev → lib → archive 开发流程)
│   ├── scratch/           # 一次性/实验性脚本
│   └── R_FUNCTION_DEV_WORKFLOW.md  # 函数开发流程文档
├── python/                # Python 脚本 (预留)
├── projects/              # 各实验/分析项目
│   ├── flow_cytometry/    # 流式细胞术 (合并 流式/ + Control_vs_R_Analysis)
│   ├── 16S/               # 16S rRNA 扩增子分析
│   ├── fish/              # FISH 荧光原位杂交
│   ├── ihc/               # 免疫组化
│   ├── qpcr/              # 定量 PCR
│   ├── lefse/             # LEfSe 差异分析
│   ├── metagenome/        # 宏基因组学
│   ├── enrichment_replot/ # 富集重画
│   ├── lysin_enzyme/      # 裂解酶活性实验
│   ├── ACS_ADS_Global/    # ACS/ADS 全局分析
│   ├── growth_curves/     # 细菌生长曲线
│   ├── clinical_tables/   # 临床 Tableone
│   ├── phage_host/        # 噬菌体-宿主互作网络
│   └── venn_diagram/      # Venn 图分析
├── collaborators/         # 协作者分析 (LHY/WY/YZX)
├── data/                  # 跨项目共享数据 (大部分数据跟随项目)
├── output/                # 跨项目共享输出
├── archive/               # 废弃/旧文件归档
├── ref/                   # 参考文献和方法论文档
└── .gitignore             # Git 忽略规则
```

## 命名规范

| 类别 | 规范 | 示例 |
|------|------|------|
| 文件夹 | 英文 snake_case | `flow_cytometry/` |
| R 脚本 | snake_case, `.R` | `compare_plot.R` |
| 数据文件 | 可保留中文原名 | `裂解酶0416.xlsx` |
| 归档版本 | 加日期后缀 | `function_v1_20250405.R` |

## 常用路径

```r
# 在项目脚本中引用库函数
source("../../R/r_functions/lib/compare_plot/compare_plot_optimized.R")
source("../../R/utils/helpers.R")
source("../../R/utils/palette_system.R")    # 颜色系统
source("../../R/utils/theme_system.R")     # 主题系统（可选）
```

## 工作流程

1. **数据分析**：在 `projects/<项目>/scripts/` 编写分析脚本
2. **函数开发**：在 `R/r_functions/dev/` 开发，成熟后移至 `R/r_functions/lib/`
3. **协作者**：在 `collaborators/<姓名>/` 下独立管理

## 版本控制

- R/Python 脚本纳入 Git 版本控制
- 输出文件 (PDF/PNG/SVG) 不纳入版本控制
- 原始数据不纳入版本控制，定期备份
- 测序原始数据在服务器上分析，不在本地存放
