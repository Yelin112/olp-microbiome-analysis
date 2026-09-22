# ============================================================================
# HE 上皮厚度 — 组间比较分析
# OXA 动物模型，分部位（上唇/舌/下唇）比较 4 组上皮厚度
# 使用 compare_plot (stable) + pub_ 主题配色系统
# ============================================================================
# 创建日期: 2025-06-24
# 更新日期: 2026-06-27 — 迁移至 compare_plot + pub_ 系统，增加统计检验 + 固定面板
# ============================================================================

# === 加载依赖 ===
library(ggplot2)
library(grid) # convertWidth / convertHeight (fix_panel_size 需要)
library(dplyr)
library(officer)
library(rvg)

# === 加载工具模块（绝对路径，无需 setwd）===
# 定位仓库根目录（本地/服务器通用）：优先 R_TOOLKIT_ROOT 环境变量，否则用 git
PROJ <- Sys.getenv("R_TOOLKIT_ROOT", unset = "")
if (!nzchar(PROJ)) PROJ <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(PROJ) != 1 || !nzchar(PROJ)) stop("找不到仓库根目录：请设置 R_TOOLKIT_ROOT 环境变量，或在仓库目录内运行")
source(file.path(PROJ, "init.R"))
setwd(
  'E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\projects\\animal_model\\OXA\\HE上皮厚度统计\\scripts'
)
source(file.path(PROJ, "utils/helpers.R")) # %||%, %ni%
source(file.path(PROJ, "utils/palette_system.R")) # get_colors(), scale_fill_pub_d()
source(file.path(PROJ, "utils/theme_system.R")) # theme_pub_base(), theme_nature()
source(file.path(PROJ, "utils/panel_fix.R")) # fix_panel_size(), render_fixed()
source(file.path(
  PROJ,
  "r_functions/lib/compare_plot/compare_plot_optimized.R"
)) # compare_plot()

# 本仓库（olp-microbiome-analysis）根目录：脚本本身在这个仓库里，用 git 定位
# （不能用 PROJ——PROJ 现在指向独立的 r-pub-toolkit 仓库，没有 projects/ 目录）
REPO_ROOT <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE), error = function(e) "")
if (length(REPO_ROOT) != 1 || !nzchar(REPO_ROOT)) stop("找不到本仓库根目录：请在仓库目录内运行脚本")

# === 读取数据 ===
raw_data <- read.csv(
  file.path(
    REPO_ROOT,
    "projects/animal_model/OXA/HE上皮厚度统计/data/上皮厚度_标准长表.csv"
  ),
  stringsAsFactors = FALSE
)

# === 数据预处理 ===
# 统一部位编码
site_labels <- c(
  "sc" = "上唇 (sc)",
  "she" = "舌 (she)",
  "xc" = "下唇 (xc)"
)
raw_data$site_label <- site_labels[raw_data$site_code]

cat("=== 数据概览 ===\n")
cat("总行数:", nrow(raw_data), "\n")
cat("分组:", paste(unique(raw_data$group), collapse = ", "), "\n")
cat("部位:", paste(unique(raw_data$site_label), collapse = ", "), "\n\n")

# 各组各部位样本量
print(table(raw_data$group, raw_data$site_label))
cat("\n")

# === 定义比较对：各处理组 vs Control ===
comparison_pairs <- list(
  c("Control", "OXA"),
  c("Control", "OXA_MN"),
  c("Control", "OXA_RM")
)

# === 分部位分析 ===
sites <- unique(raw_data$site_label)

for (site in sites) {
  cat("========================================\n")
  cat("  分析部位:", site, "\n")
  cat("========================================\n")

  site_data <- raw_data[raw_data$site_label == site, ]

  # 每组样本量
  site_n <- table(site_data$group)
  for (g in names(site_n)) {
    cat(sprintf("    %s: n = %d\n", g, site_n[g]))
  }
  cat("\n")

  # === 方案 A: 整体 ANOVA + 箱线图 + 散点 ===
  p_anova <- compare_plot(
    data = site_data,
    value.var = "thickness",
    group.by = "group",
    strategy = "auto", # 根据样本量自适应选择图型
    palette = "NPG",
    add_stat = "anova", # ★ 整体 ANOVA p 值
    stat_label = "p.format",
    hide_ns = FALSE,
    y_expand = 0.25, # 为 p 值标签留空间
    title = paste("HE 上皮厚度 —", site, "(ANOVA)"),
    xlab = "分组",
    ylab = "上皮厚度 (μm)",
    theme_use = theme_pub_base # ★ pub_ 主题系统
  )

  # === 方案 B: 成对 Wilcoxon + 箱线图 + 散点 ===
  p_wilcox <- compare_plot(
    data = site_data,
    value.var = "thickness",
    group.by = "group",
    strategy = "auto",
    palette = "NPG",
    add_stat = "wilcox.test", # ★ 成对 Wilcoxon 检验
    comparisons = comparison_pairs, # ★ Control vs 各处理组
    stat_label = "p.signif", # 星号标注
    hide_ns = TRUE, # 隐藏不显著的
    step_increase = 0.15,
    y_expand = 0.35, # 三组比较需要更大空间
    title = paste("HE 上皮厚度 —", site, "(Pairwise)"),
    xlab = "分组",
    ylab = "上皮厚度 (μm)",
    theme_use = theme_pub_base
  )

  # === 输出文件名（去括号和空格）===
  site_file <- gsub("[ ()]", "_", site)
  site_file <- gsub("_+", "_", site_file)
  site_file <- gsub("_$", "", site_file)

  # 固定面板尺寸 (cm) — ★ 所有图统一 panel 物理尺寸
  PANEL_W <- 9
  PANEL_H <- 8

  # ---- 导出函数：固定 panel → 计算总 gtable 尺寸 → ggsave ----
  save_fixed <- function(
    plot,
    file_base,
    panel_w = PANEL_W,
    panel_h = PANEL_H
  ) {
    # 1. 固定 panel 尺寸 → 返回 gtable
    g <- fix_panel_size(plot, width = panel_w, height = panel_h)

    # 2. 计算 gtable 实际需要的总宽高（panel + 标签 + 图例 + 标题）
    total_w <- convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
    total_h <- convertHeight(sum(g$heights), "cm", valueOnly = TRUE)

    # 3. PDF
    ggsave(
      paste0(file_base, ".pdf"),
      g,
      width = total_w,
      height = total_h,
      device = "pdf"
    )

    # 4. PNG
    ggsave(
      paste0(file_base, ".png"),
      g,
      width = total_w,
      height = total_h,
      dpi = 300,
      device = "png"
    )

    # 5. PPTX — 用原始 ggplot（dml 需要 ggplot 对象，且 PPTX 内可自由缩放）
    pptx <- read_pptx()
    pptx <- add_slide(
      pptx,
      layout = "Title and Content",
      master = "Office Theme"
    )
    pptx <- ph_with(pptx, dml(ggobj = plot), location = ph_location_fullsize())
    print(pptx, target = paste0(file_base, ".pptx"))

    cat(sprintf(
      "  ✓ %s.{pdf, png, pptx}  (panel: %.0f×%.0f cm)\n",
      file_base,
      panel_w,
      panel_h
    ))
  }

  # ---- 导出 ----
  save_fixed(
    p_anova,
    sprintf("../output/epithelial_thickness_%s_anova", site_file)
  )
  save_fixed(
    p_wilcox,
    sprintf("../output/epithelial_thickness_%s_wilcox", site_file)
  )
  cat("\n")
}

cat("========================================\n")
cat("  全部完成！输出文件在 ../output/\n")
cat("========================================\n")
