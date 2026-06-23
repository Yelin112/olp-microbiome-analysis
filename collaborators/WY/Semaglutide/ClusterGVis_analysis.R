# ================================================================
# ClusterGVis 基因表达趋势聚类分析
# 项目：Semaglutide 干预糖尿病小鼠多部位转录组
# 比较组：WT_Control vs DB_Control vs DB_Treated
# 部位：Lung, BO, Oral
# 方法：kmeans 聚类 + GO BP 富集分析
# ================================================================

# ────────────────────────────────────────────────────────────────
# 0. 包加载
# ────────────────────────────────────────────────────────────────

required_pkgs <- c(
  "ClusterGVis",
  "org.Mm.eg.db",
  "dplyr",
  "stringr",
  "ggsci",
  "ComplexHeatmap",
  "grid",
  "ggplot2",
  "writexl"
)

not_installed <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]
if (length(not_installed) > 0) {
  message("请先安装以下包：", paste(not_installed, collapse = ", "))
  if (
    any(not_installed %in% c("ClusterGVis", "org.Mm.eg.db", "ComplexHeatmap"))
  ) {
    cat(
      'BiocManager::install(c("ClusterGVis","org.Mm.eg.db","ComplexHeatmap"))\n'
    )
  }
  stop("缺少必要包，请安装后重新运行", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ClusterGVis)
  library(org.Mm.eg.db)
  library(dplyr)
  library(stringr)
  library(ggsci)
  library(ComplexHeatmap)
  library(grid)
  library(ggplot2)
  library(writexl)
})

set.seed(42)

# ────────────────────────────────────────────────────────────────
# 1. 路径与参数
# ────────────────────────────────────────────────────────────────

DATA_DIR <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/data/RNA/1.count"
OUT_BASE <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/analysis_results/RNA"

SITES <- c("Lung", "BO", "Oral")

# 聚类参数 —— 请先查看 elbow 图再决定是否调整
K_CLUSTERS <- 6 # 默认聚类数（可根据各部位 elbow 图单独调整）
MIN_STD <- 0.1 # clusterData 内部低方差基因过滤阈值
FPKM_CUTOFF <- 1 # 平均 FPKM 过滤阈值
ENRICH_TOPN <- 5 # 每个簇富集展示的 top term 数

# 分组对应关系（样本前缀 → 分组名）
PREFIX_TO_GROUP <- c(
  "WTN" = "WT_Control",
  "dbN" = "DB_Control",
  "dbS" = "DB_Treated"
)

# 配色方案（与主分析脚本一致）
GROUP_COLORS <- c(
  "WT_Control" = "#4DBBD5",
  "DB_Control" = "#E64B35",
  "DB_Treated" = "#F39B7F"
)

# 创建输出目录
for (site in SITES) {
  dir.create(
    file.path(OUT_BASE, site, "ClusterGVis"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ────────────────────────────────────────────────────────────────
# 2. 读取 FPKM 数据
# ────────────────────────────────────────────────────────────────

message("读取 gene_fpkm.xls ...")

fpkm_raw <- read.table(
  file.path(DATA_DIR, "gene_fpkm.xls"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1 # 第一列 gene_id 作为行名
)

# 分离样本列与注释列
is_sample_col <- grepl("_(Lung|BO|Oral)$", colnames(fpkm_raw))
fpkm_all <- fpkm_raw[, is_sample_col, drop = FALSE] # 样本矩阵
gene_anno <- fpkm_raw[, !is_sample_col, drop = FALSE] # 注释列

message(sprintf(
  "  总基因数: %d  |  总样本数: %d",
  nrow(fpkm_all),
  ncol(fpkm_all)
))

# ────────────────────────────────────────────────────────────────
# 3. 辅助函数
# ────────────────────────────────────────────────────────────────

# 保存 ComplexHeatmap / visCluster 图为 PDF
save_ht_pdf <- function(draw_func, filepath, width = 14, height = 12) {
  pdf(filepath, width = width, height = height)
  tryCatch(draw_func(), error = function(e) {
    message("  [PDF 保存警告] ", conditionMessage(e))
  })
  dev.off()
  message("  已保存: ", filepath)
}

# 保存 ComplexHeatmap / visCluster 图为 PNG
save_ht_png <- function(
  draw_func,
  filepath,
  width = 4200,
  height = 3600,
  res = 300
) {
  png(filepath, width = width, height = height, res = res)
  tryCatch(draw_func(), error = function(e) {
    message("  [PNG 保存警告] ", conditionMessage(e))
  })
  dev.off()
  message("  已保存: ", filepath)
}

# 保存 ggplot 对象为 PDF + PNG
save_gg <- function(p, base_path, width = 10, height = 8) {
  ggsave(paste0(base_path, ".pdf"), plot = p, width = width, height = height)
  ggsave(
    paste0(base_path, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
  message("  已保存: ", base_path, ".pdf / .png")
}

# ────────────────────────────────────────────────────────────────
# 4. 按部位循环分析
# ────────────────────────────────────────────────────────────────

for (site in SITES) {
  cat(sprintf(
    "\n\n══════════════════════════════════════════\n  分析部位: %s\n══════════════════════════════════════════\n",
    site
  ))

  out_dir <- file.path(OUT_BASE, site, "ClusterGVis")

  # ── 4.1 选取该部位的 3 组样本 ────────────────────────────────
  # 目标前缀：WTN（WT_Control）、dbN（DB_Control）、dbS（DB_Treated）

  site_cols <- grep(paste0("_", site, "$"), colnames(fpkm_all), value = TRUE)
  sel_cols <- site_cols[grepl("^(WTN|dbN|dbS)", site_cols)]

  if (length(sel_cols) == 0) {
    warning("  [跳过] 未找到部位 ", site, " 的样本列")
    next
  }

  # 按 WT_Control → DB_Control → DB_Treated 排序
  col_order <- c(
    sort(grep("^WTN", sel_cols, value = TRUE)),
    sort(grep("^dbN", sel_cols, value = TRUE)),
    sort(grep("^dbS", sel_cols, value = TRUE))
  )
  fpkm_site <- fpkm_all[, col_order, drop = FALSE]

  # 重命名列：去掉部位后缀，前缀统一为可读名称
  new_names <- col_order
  new_names <- gsub(paste0("_", site, "$"), "", new_names)
  new_names <- gsub("^WTN(\\d+)$", "WT_Ctrl\\1", new_names)
  new_names <- gsub("^dbN(\\d+)$", "DB_Ctrl\\1", new_names)
  new_names <- gsub("^dbS(\\d+)$", "DB_Trt\\1", new_names)
  colnames(fpkm_site) <- new_names

  n_groups <- 3 # WT_Control / DB_Control / DB_Treated
  n_reps <- length(col_order) / n_groups # 每组重复数
  group_vec <- rep(names(GROUP_COLORS), each = n_reps)

  cat(sprintf(
    "  选取样本: %d（每组 %d 个重复）\n  列名: %s\n",
    length(col_order),
    n_reps,
    paste(new_names, collapse = ", ")
  ))

  # ── 4.2 基因过滤 ────────────────────────────────────────────
  # 蛋白编码基因
  if ("gene_biotype" %in% colnames(gene_anno)) {
    coding_ids <- rownames(gene_anno)[
      !is.na(gene_anno$gene_biotype) &
        gene_anno$gene_biotype == "protein_coding"
    ]
    fpkm_site <- fpkm_site[rownames(fpkm_site) %in% coding_ids, , drop = FALSE]
    cat(sprintf("  蛋白编码基因: %d\n", nrow(fpkm_site)))
  }

  # 低表达过滤（任意组平均 FPKM > 阈值）
  group_means <- sapply(1:n_groups, function(g) {
    idx <- ((g - 1) * n_reps + 1):(g * n_reps)
    rowMeans(fpkm_site[, idx, drop = FALSE])
  })
  keep_expr <- apply(group_means, 1, max) > FPKM_CUTOFF
  fpkm_site <- fpkm_site[keep_expr, , drop = FALSE]
  cat(sprintf("  平均FPKM > %g 后保留基因: %d\n", FPKM_CUTOFF, nrow(fpkm_site)))

  # ── 4.3 log2 转换 + 替换行名为基因符号 ──────────────────────
  fpkm_log <- log2(fpkm_site + 1)

  gene_sym <- gene_anno[rownames(fpkm_log), "gene_name"]
  valid <- !is.na(gene_sym) & gene_sym != "" & gene_sym != "-" & gene_sym != "."
  gene_sym[!valid] <- rownames(fpkm_log)[!valid]
  rownames(fpkm_log) <- make.unique(as.character(gene_sym))

  cat(sprintf("  进入聚类的基因数: %d\n", nrow(fpkm_log)))

  # ── 4.3b 组内平均：将每组 n_reps 个重复合并为 1 列 ──────────
  # 最终矩阵为 3 列：WT_Control | DB_Control | DB_Treated
  fpkm_avg <- data.frame(
    WT_Control = rowMeans(fpkm_log[, seq_len(n_reps), drop = FALSE]),
    DB_Control = rowMeans(fpkm_log[, n_reps + seq_len(n_reps), drop = FALSE]),
    DB_Treated = rowMeans(fpkm_log[,
      2 * n_reps + seq_len(n_reps),
      drop = FALSE
    ]),
    check.names = FALSE
  )
  cat(sprintf("  组内平均后矩阵：%d 基因 × 3 组\n", nrow(fpkm_avg)))

  # ── 4.4 Elbow 图（确定最优聚类数）──────────────────────────
  cat("  计算 Elbow 图 ...\n")
  elbow_base <- file.path(out_dir, paste0(site, "_01_elbow"))

  p_elbow <- getClusters(obj = fpkm_avg)
  save_gg(p_elbow, elbow_base, width = 7, height = 5)

  # ── 4.5 clusterData（kmeans 聚类）──────────────────────────
  cat(sprintf("  执行 kmeans 聚类（k = %d）...\n", K_CLUSTERS))
  ck <- clusterData(
    obj = fpkm_avg,
    clusterMethod = "kmeans",
    clusterNum = K_CLUSTERS,
    minStd = MIN_STD,
    scaleData = TRUE
  )

  cat("  各簇基因数：\n")
  print(table(ck$wide.res$cluster))

  # ── 4.6 GO BP 富集分析 ──────────────────────────────────────
  cat("  GO BP 富集分析 ...\n")
  enrich <- tryCatch(
    enrichCluster(
      object = ck,
      OrgDb = org.Mm.eg.db,
      type = "BP",
      pvalueCutoff = 0.05,
      topn = ENRICH_TOPN,
      readable = TRUE
    ),
    error = function(e) {
      message("  [富集警告] ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(enrich) && nrow(enrich) > 0) {
    cat(sprintf("  富集 term 数: %d\n", nrow(enrich)))
    write_xlsx(
      enrich,
      file.path(out_dir, paste0(site, "_GO_BP_enrichment.xlsx"))
    )
  } else {
    cat("  富集无显著结果（将跳过 GO 注释）\n")
    enrich <- NULL
  }

  # ── 4.7 导出聚类结果表 ────────────────────────────────────
  cluster_export <- list(
    cluster_summary = ck$wide.res %>%
      dplyr::select(gene, cluster) %>%
      left_join(
        gene_anno %>%
          tibble::rownames_to_column("gene_id") %>%
          dplyr::select(gene_id, gene_name, gene_biotype, gene_description) %>%
          mutate(gene_description = substr(gene_description, 1, 500)),
        by = c("gene" = "gene_name")
      ),
    cluster_gene_list = do.call(
      rbind,
      lapply(names(ck$cluster.list), function(nm) {
        data.frame(cluster = nm, gene = ck$cluster.list[[nm]])
      })
    )
  )
  write_xlsx(
    cluster_export,
    file.path(out_dir, paste0(site, "_cluster_gene_list.xlsx"))
  )

  # ── 4.8 可视化：折线趋势图 ──────────────────────────────────
  cat("  绘制折线趋势图 ...\n")

  line_base <- file.path(out_dir, paste0(site, "_02_trend_line"))

  p_line <- visCluster(
    object = ck,
    plotType = "line",
    ncol = 3,
    lineSize = 0.08,
    lineCol = "grey80",
    addMline = TRUE,
    mlineSize = 1.5,
    mlineCol = "#333333" # 每簇单条均值线（3列已是组均值）
  )
  save_gg(p_line, line_base, width = 14, height = 10)

  # ── 4.9 可视化：热图 + GO 注释（both 模式）──────────────────
  cat("  绘制热图 + 趋势图组合（both 模式）...\n")

  # 顶部样本分组注释（3 列：每组 1 列均值）
  top_anno <- HeatmapAnnotation(
    Group = names(GROUP_COLORS),
    col = list(Group = GROUP_COLORS),
    gp = gpar(col = "white"),
    annotation_name_gp = gpar(fontsize = 9)
  )

  # 绘图参数 - 用于多次调用
  # 注意：show_row_names 由 visCluster 内部控制，不在此传入以避免重复参数报错
  vis_args_base <- list(
    object = ck,
    plotType = "both",
    lineSide = "left",
    # 热图颜色
    htColList = list(
      col_range = c(-2, 0, 2),
      col_color = c("#2166AC", "white", "#B2182B")
    ),
    # 顶部注释
    heatmapAnnotation = top_anno,
    # 趋势图：3列已是组均值，直接显示均值线即可
    mlineCol = "#333333",
    mlineSize = 1.2,
    addBox = FALSE,
    addLine = TRUE,
    # 簇颜色
    ctAnnoCol = pal_npg()(K_CLUSTERS),
    # ComplexHeatmap 参数
    column_names_rot = 45,
    show_row_dend = FALSE,
    border = TRUE
  )

  # 无 GO 注释版本
  both_base_noGO <- file.path(out_dir, paste0(site, "_03_both_noGO"))

  save_ht_pdf(
    function() do.call(visCluster, vis_args_base),
    paste0(both_base_noGO, ".pdf"),
    width = 12,
    height = 12
  )
  save_ht_png(
    function() do.call(visCluster, vis_args_base),
    paste0(both_base_noGO, ".png"),
    width = 3600,
    height = 3600,
    res = 300
  )

  # 有 GO 注释版本
  if (!is.null(enrich)) {
    cat("  绘制带 GO 注释的组合图 ...\n")

    vis_args_go <- c(
      vis_args_base,
      list(
        annoTermData = enrich,
        annoTermMside = "right",
        goCol = rep(pal_d3()(K_CLUSTERS), each = ENRICH_TOPN),
        goSize = "pval",
        byGo = "anno_link",
        wordWrap = TRUE,
        termTextLimit = c(10, 20),
        termAnnoArg = c("grey95", "grey50"),
        addBar = TRUE,
        barWidth = 6
      )
    )

    both_base_GO <- file.path(out_dir, paste0(site, "_04_both_withGO"))

    save_ht_pdf(
      function() do.call(visCluster, vis_args_go),
      paste0(both_base_GO, ".pdf"),
      width = 18,
      height = 12
    )
    save_ht_png(
      function() do.call(visCluster, vis_args_go),
      paste0(both_base_GO, ".png"),
      width = 5400,
      height = 3600,
      res = 300
    )
  }

  # ── 4.10 可视化：仅热图（便于论文配图）──────────────────────
  cat("  绘制单独热图 ...\n")

  ht_args <- list(
    object = ck,
    plotType = "heatmap",
    htColList = list(
      col_range = c(-2, 0, 2),
      col_color = c("#2166AC", "white", "#B2182B")
    ),
    heatmapAnnotation = top_anno,
    ctAnnoCol = pal_npg()(K_CLUSTERS),
    column_names_rot = 45,
    show_row_dend = FALSE,
    border = TRUE
  )

  if (!is.null(enrich)) {
    ht_args <- c(
      ht_args,
      list(
        annoTermData = enrich,
        annoTermMside = "right",
        goCol = rep(pal_d3()(K_CLUSTERS), each = ENRICH_TOPN),
        goSize = "pval"
      )
    )
  }

  ht_base <- file.path(out_dir, paste0(site, "_05_heatmap"))

  save_ht_pdf(
    function() do.call(visCluster, ht_args),
    paste0(ht_base, ".pdf"),
    width = 10,
    height = 12
  )
  save_ht_png(
    function() do.call(visCluster, ht_args),
    paste0(ht_base, ".png"),
    width = 3000,
    height = 3600,
    res = 300
  )

  cat(sprintf("\n  ✓ %s 分析完成，结果保存至: %s\n", site, out_dir))
}

cat(
  "\n\n════════════════════════════════════\n  全部分析完成！\n════════════════════════════════════\n"
)
cat("\n注意事项：\n")
cat("  1. 请查看各部位的 *_01_elbow 图，若聚类数明显不合适，\n")
cat("     请修改 K_CLUSTERS 参数后重新运行。\n")
cat("  2. 带 GO 注释的图需要网络请求 org.Mm.eg.db，\n")
cat("     若失败可适当增大 pvalueCutoff 或减小 ENRICH_TOPN。\n")
cat("  3. 所有结果文件均保存在各部位的 ClusterGVis/ 子目录中。\n")
