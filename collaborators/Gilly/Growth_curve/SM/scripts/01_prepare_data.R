# =============================================================================
# 生长曲线（OD600）数据准备 —— 单批次
#
# 用法：
#   1) 改下面的 BATCH_ID 为要处理的批次文件夹名（对应 raw/<BATCH_ID>/），然后整个
#      source 运行；或者
#   2) 命令行：Rscript 01_prepare_data.R <batch_id>（不改脚本，参数会覆盖 BATCH_ID）
#
# 要求 raw/<BATCH_ID>/ 下有 growth_curve_tidy_A01_F03.csv 和
# growth_curve_summary_by_group.csv 两张表（长格式 tidy + 按组按时间点汇总），
# 列结构须与 batch1（见 raw/batch1/positron_claude_handoff.md）一致：
#   tidy: read_index, time_sec, time_h, group, well, od600, blank_mean, od600_blank_corrected
#   summary: group, read_index, time_h, od600_mean, od600_sd,
#            od600_blank_corrected_mean, od600_blank_corrected_sd, n
# 若某批次拿到的是原始仪器宽表（wide format，如 raw/batch2/growth_curve.tsv），需要
# 先转换成上述 tidy/summary 两张表（G01-G03 视为空白孔背景扣除），再指向这里。
#
# 背景：96 孔板 OD600 生长曲线，评估同一纳米酶材料对细菌生长的抑制作用
#       （与 collaborators/Gilly/CCK8 是同一材料体系的配套实验）。
#       每 30 min 读数一次，共 28 个时间点（0~13.5 h）。
#       A~F 六组，每组 3 复孔（01~03），G01~03 为空白孔（仅培养基，无细菌），
#       已在各批次数据中用于背景扣除：od600_blank_corrected = od600 - blank_mean(G组)。
#
# 分组浓度映射（用户 2026-08-28 确认，所有批次共用同一套材料浓度体系，定义在
# _common.R 的 conc_map 里）：A = 阳性对照，B~F 为同一材料的 2 倍稀释系列
# （200 / 100 / 50 / 25 / 12.5 ug/mL）。
# =============================================================================

# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJ <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJ)) PROJ <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJ) != 1 || !nzchar(PROJ)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJ, "init.R"))
source(file.path(PROJ, "collaborators/Gilly/Growth_curve/SM/scripts/_common.R"))

BATCH_ID <- "batch3" # <-- 手动改成要处理的批次，例如 "batch2" / "batch3"
BATCH_ID <- resolve_batch_id(BATCH_ID)

raw_dir <- file.path(RAW_DIR, BATCH_ID)
data_dir <- file.path(DATA_DIR, BATCH_ID)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "raw/<BATCH_ID>/growth_curve_tidy_A01_F03.csv 不存在，请检查 BATCH_ID 或先准备好该批次的 tidy 表" = file.exists(file.path(
    raw_dir,
    "growth_curve_tidy_A01_F03.csv"
  )),
  "raw/<BATCH_ID>/growth_curve_summary_by_group.csv 不存在，请检查 BATCH_ID" = file.exists(file.path(
    raw_dir,
    "growth_curve_summary_by_group.csv"
  ))
)

tidy_raw <- read.csv(
  file.path(raw_dir, "growth_curve_tidy_A01_F03.csv"),
  stringsAsFactors = FALSE
)

tidy_data <- tidy_raw %>%
  left_join(conc_map, by = "group") %>%
  mutate(
    group = factor(group, levels = conc_map$group, labels = conc_map$label),
    well = factor(well)
  )
write.csv(
  tidy_data,
  file.path(data_dir, "growth_curve_tidy.csv"),
  row.names = FALSE
)

summary_raw <- read.csv(
  file.path(raw_dir, "growth_curve_summary_by_group.csv"),
  stringsAsFactors = FALSE
)
summary_data <- summary_raw %>%
  left_join(conc_map, by = "group") %>%
  mutate(
    group = factor(group, levels = conc_map$group, labels = conc_map$label)
  )
write.csv(
  summary_data,
  file.path(data_dir, "growth_curve_summary.csv"),
  row.names = FALSE
)

# --- 逐孔终点 OD600（末次读数）+ 逐孔 AUC（梯形积分），用于组间统计检验（每组 n=3） ---
endpoint_well <- tidy_raw %>%
  filter(time_h == max(time_h)) %>%
  transmute(
    group,
    well,
    final_od600 = od600,
    final_od600_blank_corrected = od600_blank_corrected
  )

auc_well <- tidy_raw %>%
  arrange(group, well, time_h) %>%
  group_by(group, well) %>%
  summarise(
    auc_od600 = trapz(time_h, od600),
    auc_od600_blank_corrected = trapz(time_h, od600_blank_corrected),
    .groups = "drop"
  )

replicate_metrics <- endpoint_well %>%
  left_join(auc_well, by = c("group", "well")) %>%
  left_join(conc_map, by = "group") %>%
  mutate(
    group = factor(group, levels = conc_map$group, labels = conc_map$label),
    well = factor(well)
  ) %>%
  arrange(group, well)
write.csv(
  replicate_metrics,
  file.path(data_dir, "replicate_endpoint_auc.csv"),
  row.names = FALSE
)

message("✓ [", BATCH_ID, "] 数据准备完成，输出至: ", data_dir)
message("  - growth_curve_tidy.csv (", nrow(tidy_data), " 行)")
message("  - growth_curve_summary.csv (", nrow(summary_data), " 行)")
message(
  "  - replicate_endpoint_auc.csv (",
  nrow(replicate_metrics),
  " 行, 每组 n=3)"
)
message(
  "下一步：Rscript 02_analyze_batch.R ",
  BATCH_ID,
  "（或改 02_analyze_batch.R 里的 BATCH_ID 后 source）"
)
