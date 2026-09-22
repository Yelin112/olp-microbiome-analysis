# 01_read_explore.R
# 目标：读取 S1/ 下所有 FCS 文件，建立通道名 → 抗体标记映射表，识别重复文件

library(flowCore)
library(tidyverse)

# ── 路径（相对于项目根目录 分析/）─────────────────────────────────────────────
DATA_DIR <- "R/scratch/flow_cyto/data/S1"
OUTPUT_DIR <- "R/scratch/flow_cyto/output"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fcs_files <- list.files(DATA_DIR, pattern = "\\.fcs$", full.names = TRUE)

# ── 从单个 FCS 提取通道信息 ───────────────────────────────────────────────────
extract_channel_info <- function(path) {
  fcs <- read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE)
  params <- pData(parameters(fcs))
  expr <- exprs(fcs)

  tibble(
    file = basename(path),
    n_events = nrow(expr),
    channel = params$name, # 仪器通道，如 "PE-A"
    label = params$desc, # 抗体标记，如 "CD3"；可能为 NA
    val_min = apply(expr, 2, min),
    val_max = apply(expr, 2, max)
  )
}

# ── 批量读取（进度只打到控制台，不写入结果文件）─────────────────────────────
cat("正在读取 FCS 文件...\n")
channel_df <- map_dfr(fcs_files, \(p) {
  cat(sprintf("  %s\n", basename(p)))
  extract_channel_info(p)
})
cat(sprintf("读取完成，共 %d 个文件。\n\n", length(fcs_files)))

# ── 计算衍生表 ────────────────────────────────────────────────────────────────
SCATTER_CHANNELS <- c(
  "FSC-A",
  "FSC-H",
  "FSC-W",
  "SSC-A",
  "SSC-H",
  "SSC-W",
  "Time"
)

# 重复文件检测
all_names <- basename(fcs_files)
stripped <- sub("_\\d{8}_\\d{6}\\.fcs$", ".fcs", all_names)
dup_table <- tibble(file = all_names, stripped = stripped) |>
  filter(duplicated(stripped) | duplicated(stripped, fromLast = TRUE)) |>
  arrange(stripped)

# 荧光通道汇总（每文件×通道）
panel_summary <- channel_df |>
  filter(!channel %in% SCATTER_CHANNELS) |>
  distinct(file, channel, label) |>
  arrange(file, channel)

# 跨文件通道频次
channel_freq <- channel_df |>
  filter(!channel %in% SCATTER_CHANNELS) |>
  count(channel, label, name = "n_files") |>
  arrange(desc(n_files), channel)

# ── 写结果文件 ────────────────────────────────────────────────────────────────

# 1. 文字报告（channel_report.txt）
report_path <- file.path(OUTPUT_DIR, "channel_report.txt")
sink(report_path)

cat(sprintf(
  "FCS 通道探索报告\n生成时间: %s\n数据目录: %s\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  DATA_DIR
))

cat(
  "\n── 各 FCS 文件通道详情 ────────────────────────────────────────────────\n"
)
channel_df |>
  group_split(file) |>
  walk(\(df) {
    cat(sprintf(
      "\n[%s]  事件数 = %s\n",
      df$file[1],
      format(df$n_events[1], big.mark = ",")
    ))
    df |> select(channel, label, val_min, val_max) |> print(n = Inf)
  })

cat(
  "\n── 重复文件检查 ────────────────────────────────────────────────────────\n"
)
if (nrow(dup_table) > 0) {
  cat("发现可能重复的文件组:\n")
  print(dup_table)
  cat("建议：确认后只保留一个版本，将原始 FCS 整理到 data/S1/raw/\n")
} else {
  cat("未发现重复文件\n")
}

cat("\n── 各通道出现频次（n_files = 出现在几个 FCS 文件中）─────────────────\n")
print(channel_freq, n = Inf)

sink()

# 2. 荧光通道映射表（panel_summary.csv）——用于填写 config
write_csv(panel_summary, file.path(OUTPUT_DIR, "panel_summary.csv"))

# 3. 完整通道数据（channel_detail.csv）——包含 val_min/val_max
write_csv(channel_df, file.path(OUTPUT_DIR, "channel_detail.csv"))

# ── 控制台简报 ────────────────────────────────────────────────────────────────
cat(sprintf(
  "结果已写入 %s/\n  channel_report.txt  — 文字报告\n  panel_summary.csv   — 荧光通道映射表\n  channel_detail.csv  — 完整通道数据\n",
  OUTPUT_DIR
))
if (nrow(dup_table) > 0) {
  cat(sprintf(
    "\n注意：发现 %d 个可能重复的文件，详见 channel_report.txt\n",
    nrow(dup_table)
  ))
}
