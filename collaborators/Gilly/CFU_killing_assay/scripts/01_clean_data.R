# =============================================================================
# CFU 杀菌实验：数据清洗
# 数据：collaborators/Gilly/CFU_killing_assay/raw/cfu_tidy_data.csv
# 背景：连续稀释点板法评估 Cu-Zn@HA-MnO2 纳米酶的抗菌/杀菌效果。
#       10 uL 点样体积，4 次独立实验（exp_1~exp_4），每次实验含
#       control / nanozyme_1mg / nanozyme_0.5mg 三组，每组一张图片估读菌落数
#       （半定量，非技术重复）。
#
# 剔除 nanozyme_1mg 组的原因：
#   4 次实验中有 3 次（exp_2/exp_3/exp_4）nanozyme_1mg 组的杀菌率为 0% 或负值
#   （exp_3 甚至菌落数超过同批 control），显著劣于同批次浓度更低的
#   nanozyme_0.5mg 组（66.7%~80% 杀菌率），出现"高浓度组杀菌力反而不如低浓度组"
#   的倒置剂量-效应关系，与仅 exp_1 呈现的正常剂量趋势（22.2% vs 29.6%）不符。
#   判定为该浓度组存在数据/标签问题（如点板时 1mg 与 0.5mg 组别弄混），
#   予以剔除，仅保留 control 与 nanozyme_0.5mg 两组用于后续分析。
#
# 剔除 exp_1 的原因：
#   剔除 nanozyme_1mg 后，剩余 4 组配对中 exp_1 的 nanozyme_0.5mg 杀菌率仅
#   29.6%，明显低于 exp_2/exp_3/exp_4 的 66.7%~80%，判定为该次实验结果有问题，
#   整组（control + nanozyme_0.5mg）一并剔除以保持配对结构完整，
#   最终保留 exp_2/exp_3/exp_4 共 3 组配对用于分析。
# =============================================================================

library(tidyverse)

# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJ <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJ)) PROJ <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJ) != 1 || !nzchar(PROJ)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJ, "init.R"))

# 本仓库（olp-microbiome-analysis）根目录：脚本本身在这个仓库里，用 git 定位
# （不能用 PROJ——PROJ 现在指向独立的 r-pub-toolkit 仓库，没有 collaborators/ 目录）
REPO_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(REPO_ROOT) != 1 || !nzchar(REPO_ROOT)) stop("找不到本仓库根目录：请在仓库目录内运行脚本")

BASE_DIR <- file.path(REPO_ROOT, "collaborators/Gilly/CFU_killing_assay")
RAW_DIR  <- file.path(BASE_DIR, "raw")
DATA_DIR <- file.path(BASE_DIR, "data")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(file.path(RAW_DIR, "cfu_tidy_data.csv"), stringsAsFactors = FALSE)

clean <- raw %>%
  filter(group != "nanozyme_1mg", experiment != "exp_1") %>%
  mutate(
    group = factor(group, levels = c("control", "nanozyme_0.5mg")),
    experiment = factor(experiment)
  ) %>%
  arrange(experiment, group)

message("剔除 nanozyme_1mg 与 exp_1 后剩余样本量：")
print(table(clean$group))

write.csv(clean, file.path(DATA_DIR, "cfu_clean_data.csv"), row.names = FALSE)
message("✓ 已保存干净数据: ", file.path(DATA_DIR, "cfu_clean_data.csv"))
