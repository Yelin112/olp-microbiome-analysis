# CLAUDE.md — 项目导航与工作约定

## 项目背景

OLP（口腔扁平苔藓）微生物组研究。实验涉及流式细胞术、16S rRNA 扩增子、宏基因组、qPCR、FISH、IHC、裂解酶活性、细菌生长曲线等。测序原始数据在服务器上分析，本地只放仪器导出的数据（Excel/CSV）和服务器下载的处理后结果。

## 两仓库架构

通用 R 工具（主题系统、配色系统、compare_plot/pcoa_plot 等可复用函数）已拆到独立仓库
**[r-pub-toolkit](https://github.com/Yelin112/r-pub-toolkit)**，不再绑定本课题——它是
跨项目通用的画图/分析工具，任何课题（包括帮合作者做的分析）都能引用，不必依赖本仓库。

- **本仓库（olp-microbiome-analysis）**：OLP 课题自己的项目、数据、脚本、合作者分析
- **r-pub-toolkit**：主题系统、配色系统、`compare_plot`/`pcoa_plot`/`abundance_processing`/
  `network_analysis` 等可复用函数库，`init.R`/`check_deps.R`/`check_conventions.R` 等基础设施

两个仓库都通过 GitHub 在三台机器（本地 Windows + 2 台 Linux 服务器）上同步，各自独立
clone、独立 `git pull`/`push`（不是 submodule，是两个平行仓库）。

## 目录结构

```
分析/                                # 本仓库 = olp-microbiome-analysis
├── R/
│   └── scratch/                    # 一次性/实验性脚本（OLP 专属，不进 r-pub-toolkit）
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
│   ├── scRNA_seq/                  # 单细胞 RNA-seq
│   └── ACS_ADS_Global/             # ACS/ADS 全局分析
├── collaborators/                  # 帮其他人做的分析（LHY/WY/YZX/Gilly/ZSY 等）
├── data/metadata/                  # 跨项目共用的样本信息
├── archive/                        # 废弃/旧文件（main.r 等）
├── python/                         # Python 脚本预留
├── ref/                            # 参考文献
└── output/                         # 跨项目输出

r-pub-toolkit/                      # 独立仓库，与本仓库平级 clone，不在本仓库内
├── init.R                          # 统一入口：定位工具库根 + load_utils() / load_lib()
├── check_deps.R                    # 检查当前机器缺哪些 R 包
├── check_conventions.R             # 检查新函数是否符合规范
├── R_FUNCTION_DEV_WORKFLOW.md      # 函数开发流程详细文档
├── PLOTTING_CONVENTIONS.md         # 分析图表产出规范（导出格式/英文标签/统计设计）
├── utils/                          # 通用工具函数（主题系统、颜色系统、面板固定等）
│   ├── helpers.R                    # 零依赖基础工具 (%||%, %ni%)
│   ├── palette_system.R             # 颜色唯一真相来源 (调色板库 + 解析引擎 + scale_xxx_pub_*)
│   ├── theme_system.R               # 主题工厂 (theme_pub_base → 期刊预设)
│   ├── panel_fix.R                  # 固定面板尺寸工具
│   ├── FUNCTION_TEMPLATE.R          # 新可视化函数开发模板
│   └── README.md                    # 主题与配色系统完整使用手册 ⭐
└── r_functions/                    # 正式函数库 dev → lib → archive 开发流程
    ├── dev/                        # 开发中的函数
    ├── lib/                        # 稳定函数（项目脚本的 source 目标）
    └── archive/                    # 旧版本备份
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
4. 在脚本中引用 r-pub-toolkit 的库函数：
```r
source(file.path(Sys.getenv("R_TOOLKIT_ROOT"), "init.R"))   # 定位 r-pub-toolkit 仓库根目录
load_utils(theme = TRUE)                                      # helpers + palette（theme=TRUE 加主题；pptx=TRUE 加 PPTX 导出）
load_lib("compare_plot")                                      # 加载 lib 下的分类目录或单个函数
```

不要再写 `../../R/...` 这类依赖目录深度的相对路径，也不要用 `sys.frame(1)$ofile`（Rscript 下失效）。
`R_TOOLKIT_ROOT` 指向 r-pub-toolkit 仓库（跨仓库调用，必须显式设置环境变量，见下文"多机协作"）。

### 模式 2：帮合作者分析
在 `collaborators/<姓名>/` 下独立管理，每个人有独立的 data/scripts/figures，同样通过
`R_TOOLKIT_ROOT` 引用 r-pub-toolkit（合作者项目和 OLP 课题本身共用同一套工具库，
这正是拆出独立仓库的原因）。

### 模式 3：R 函数开发
函数开发在 **r-pub-toolkit 仓库**（不是本仓库）里进行，完整流程、编码约定、各阶段
检查要求见它自己的 [CLAUDE.md](https://github.com/Yelin112/r-pub-toolkit/blob/master/CLAUDE.md)
和 `R_FUNCTION_DEV_WORKFLOW.md`——这里不重复，避免两份文档各自漂移。

本仓库的 `R/scratch/` 保留，放 OLP 专属的一次性/探索性脚本（不是通用工具，不进
r-pub-toolkit；一次性脚本里如果沉淀出值得复用的逻辑，再搬到 r-pub-toolkit 的
`dev/` 走正式开发流程）。

## 多机协作（本地 Windows + 2 台 Linux 服务器）

**两个仓库**（本仓库 + r-pub-toolkit）在三台机器上各有一份 clone，通过 GitHub 独立同步；
数据不入库，各机自己放。三台机器上都会用 Claude Code 直接改代码，所以：

- 开始工作前先对**两个仓库**分别 `git pull`，结束后分别 commit + push，不要带着未推送的改动离开
- 函数开发（r-pub-toolkit 的 `dev/`）在本地进行；服务器上只改 `lib/` 里的 bug 和小调整，改完立刻 push
- 同一个函数不要在两台机器上同时改
- 每台机器在 `~/.Renviron` 中设置 `R_TOOLKIT_ROOT=<r-pub-toolkit 仓库路径>`；未设置时
  `init.R` 会回退到 git / 向上查找（仅当调用脚本本身在 r-pub-toolkit 仓库内时才有效，
  跨仓库调用——也就是本仓库里几乎所有脚本——必须显式设置这个环境变量）
- 已有分析脚本开头是 4 行"引导"（`R_TOOLKIT_ROOT` 优先，回退 `git rev-parse`，再
  `source(init.R)`），新脚本可直接用上面的一行写法；两种都不要写死 `E:/...` 绝对路径
- 服务器 RStudio Server 同样适用：`~/.Renviron` 改动后需重启 R 会话（Session → Restart R）
  才生效；RStudio 默认工作目录是家目录，不在任何仓库内，此时必须依赖 `R_TOOLKIT_ROOT`。
  RStudio 用的 R 可能与命令行 `Rscript` 不是同一个，需在 RStudio 控制台里也跑一次
  `check_deps()`
- 环境为系统 R；换机或首次 clone 后在 r-pub-toolkit 仓库里运行 `Rscript check_deps.R`
  检查缺失的包（只报告，不自动安装）。新函数引入新依赖时，同步更新 `check_deps.R` 的
  `DEPS` 表
- 服务器缺 `officer`/`rvg` 时，PPTX 导出会警告并跳过，PDF/PNG 照常输出，不要为此中断脚本
- 分工：服务器做重计算（标准化、差异分析、beta 多样性等），只导出中间表（CSV/RDS）；出图微调和 PPTX 回本地做
- 服务器直连 GitHub 可能超时，需要走本地代理（如 mihomo，常见监听 `127.0.0.1:7890`）；
  给对应仓库配 `git config http(s).proxy http://127.0.0.1:7890` 后 `git pull`/`push` 才通

## 技术栈

- **主语言**：R（ggplot2 画图，tidyverse 数据处理）
- **辅助**：Python（偶尔，在 `python/` 下）
- **版本控制**：Git，已关联 GitHub remote origin；另有独立的 r-pub-toolkit 仓库
- **AI 工具**：Claude Code（`.claude/settings.local.json`）

## 跨项目共享知识库

`~/projects/research-commons`（GitHub: `Yelin112/research-commons`）是跨项目共享的知识库，收纳函数库、流程模板、方法笔记、工具手册。做分析前先去那里查有没有现成的，避免重复摸索；本项目里定型的可复用资产（函数、流程参数组合、方法论笔记）也应该按 `anthropic-skills:commons-harvest` skill 的规范沉淀过去，而不是只留在本项目里。

与本项目直接相关的资料：
- `references/biobakery/bioBakery-overview.md` —— bioBakery工具全景手册（HUMAnN/MetaPhlAn/StrainPhlAn/MaAsLin等），含 OLP口腔微生物组相关性速查、以及 Assembly workflow(sgb_pipeline) vs 手工 MEGAHIT+Prodigal+CD-HIT+eggNOG-mapper 流程的取舍结论。做宏基因组物种/功能谱分析、菌株分析、差异丰度检验前应先查这份手册。

云端 Claude Code 会话默认看不到 research-commons，需要显式 attach（`add_repo` owner=Yelin112 repo=research-commons）；本地/SSH 会话直接用文件路径访问。r-pub-toolkit 同理，云端会话需要显式 attach 才能看到。

## .gitignore 规则

脚本纳入 Git；PDF/PNG/SVG 输出不入库；data/metadata 不入库；RData/Rhistory 不入库。

## 主题与配色系统

所有 ggplot2 绘图的颜色和主题通过 **r-pub-toolkit** 的 `utils/` 下的统一系统管理。**在任何脚本或函数内不定义调色板列表**——颜色永远从 `palette_system.R` 获取。

### 日常绘图快速启动

```r
source(file.path(Sys.getenv("R_TOOLKIT_ROOT"), "init.R"))
load_utils()

ggplot(df, aes(x, y, fill = group)) +
  geom_boxplot() +
  scale_fill_pub_d("NPG")       # 离散 fill 配色（也支持颜色向量）
```

### 配色 API (pub = publication)

- `scale_fill_pub_d("NPG")` — 离散 fill
- `scale_color_pub_d("NPG")` — 离散 color
- `scale_fill_pub_c("RWB")` — 连续 fill
- `scale_color_pub_c("RWB")` — 连续 color
- `get_colors("NPG", n = 5)` — 直接取颜色向量
- `set_palette("name", colors)` — 运行时注册配色
- `preview_palette("NPG")` — 预览

### 主题 API

- `theme_pub_base()` — 基石主题
- `theme_pub_stat()` — 统计图
- `theme_nature()` / `theme_cell()` / `theme_jama()` / `theme_science()` — 期刊预设

### 开发新绘图函数

在 r-pub-toolkit 仓库里复制 `utils/FUNCTION_TEMPLATE.R`，函数签名必须包含 `palette` 和 `theme_use` 参数。
配色用 `scale_fill_pub_d(palette)` 对接，主题用 `if (is.function(theme_use)) p + theme_use()` 对接。

### 详细文档

完整使用手册、自定义指南、速查卡片 → r-pub-toolkit 仓库的 `utils/README.md`

### 分析图表产出规范

导出格式（PDF/PNG/PPTX）、图内文字用英文、组间比较统计设计惯例 → r-pub-toolkit 仓库的 `PLOTTING_CONVENTIONS.md`

## 注意事项

- 不要修改文件内容，只动目录结构（约定）
- 新分析优先在 `projects/` 下建独立项目，不要往根目录扔脚本
- 复用已有库函数（r-pub-toolkit），避免重复造轮子
- 函数开发在 r-pub-toolkit 仓库进行，参考其 `R_FUNCTION_DEV_WORKFLOW.md` 中的模板
- **不要在脚本或函数内定义调色板列表**，用 `palette_system.R` 的 `get_colors()` 或 `scale_fill_pub_d()`
- 本仓库不再包含 `R/utils`、`R/r_functions`、`R/init.R` 等——那些都在 r-pub-toolkit 仓库；
  本仓库 `R/` 下只剩 `scratch/`（OLP 专属一次性脚本）
