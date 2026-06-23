# ================================================================
# 转录组下游分析脚本
# 项目：Semaglutide 干预糖尿病小鼠多部位转录组分析
# 分析部位：Lung、BO（肠道活检，对应16S的Gut）、Oral
# 4组设计：WT_Control / WT_Treated / DB_Control / DB_Treated
# 差异分析：DESeq2 Wald（两两对比）+ LRT（多组整体检验）
# 富集分析：clusterProfiler GO BP / KEGG + fgsea GSEA
# 联合分析接口：输出VST矩阵、metadata、LRT基因集供16S整合
# ================================================================

# ────────────────────────────────────────────────────────────
# 0. 包检查与加载
# ────────────────────────────────────────────────────────────

required_pkgs <- c(
  "DESeq2",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "RColorBrewer",
  "dplyr",
  "tidyr",
  "stringr",
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "fgsea",
  "msigdbr"
)

not_installed <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]

if (length(not_installed) > 0) {
  cat("── 请先安装以下包后重新运行 ──────────────────────\n")
  cat('if (!requireNamespace("BiocManager")) install.packages("BiocManager")\n')
  cat(sprintf(
    'BiocManager::install(c("%s"))\n',
    paste(not_installed, collapse = '", "')
  ))
  stop("缺少必要包，请安装后重新运行", call. = FALSE)
}

suppressPackageStartupMessages(
  for (pkg in required_pkgs) {
    library(pkg, character.only = TRUE)
  }
)

set.seed(42)
options(warn = 1)

# ────────────────────────────────────────────────────────────
# 1. 路径设置
# ────────────────────────────────────────────────────────────

DATA_DIR <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/data/RNA/1.count"
OUT_DIR <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/analysis_results/RNA"

# 创建输出目录
for (d in c(
  file.path(OUT_DIR, c("Lung", "BO", "Oral")),
  file.path(OUT_DIR, c("Lung", "BO", "Oral"), "GSEA"),
  file.path(OUT_DIR, "Integration")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("输出根目录:", OUT_DIR, "\n\n")

# ────────────────────────────────────────────────────────────
# 2. 分析参数（集中管理，便于调整）
# ────────────────────────────────────────────────────────────

P <- list(
  min_count = 10, # 低表达过滤：单个样本中的最低count
  min_samples = 3, # 至少有N个样本通过min_count
  padj_cutoff = 0.05, # 差异显著性阈值
  lfc_cutoff = 1.0, # Volcano图上下调着色阈值（log2FC）
  heatmap_n = 60, # 热图最多显示N个基因
  gsea_padj = 0.25, # GSEA显著性阈值（常用宽松阈值）
  gsea_min_size = 15,
  gsea_max_size = 500,
  gsea_nperm = 10000
)

# ────────────────────────────────────────────────────────────
# 3. 读取数据
# ────────────────────────────────────────────────────────────

message("读取 gene_count.xls ...")

count_raw <- read.table(
  file.path(DATA_DIR, "gene_count.xls"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1 # 第一列 gene_id 作为行名
)

# 识别样本列（后缀为 _Lung / _BO / _Oral）和注释列
is_sample_col <- grepl("_(Lung|BO|Oral)$", colnames(count_raw))
count_matrix <- count_raw[, is_sample_col, drop = FALSE]
gene_anno_raw <- count_raw[, !is_sample_col, drop = FALSE]

# 为注释表补回 gene_id 列（方便后续 left_join）
gene_anno <- gene_anno_raw
gene_anno$gene_id <- rownames(gene_anno)
gene_anno$gene_name <- as.character(gene_anno$gene_name)

message(sprintf(
  "  基因总数: %d  | 样本总数: %d",
  nrow(count_matrix),
  ncol(count_matrix)
))
message(sprintf(
  "  注释列: %s",
  paste(colnames(gene_anno_raw), collapse = ", ")
))

# ────────────────────────────────────────────────────────────
# 4. 构建 RNA 样本 metadata
# ────────────────────────────────────────────────────────────
# 命名规则：{db/WT}{N/S}{rep}_{Site}
#   db  = diabetic (DB)；  WT = wild-type
#   N   = Normal/No treatment (Control)；  S = Semaglutide (Treated)

rna_meta <- data.frame(
  rna_sample = colnames(count_matrix),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Site = str_extract(rna_sample, "Lung|BO|Oral"),
    prefix = str_remove(rna_sample, "_(Lung|BO|Oral)$"),
    Condition = case_when(
      str_starts(prefix, "WT") ~ "WT",
      str_starts(prefix, "db") ~ "DB",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      str_detect(prefix, "N\\d+$") ~ "Control",
      str_detect(prefix, "S\\d+$") ~ "Treated",
      TRUE ~ NA_character_
    ),
    Replicate = str_extract(prefix, "\\d+$"),
    Group = paste(Condition, Treatment, sep = "_"),
    # 对应 16S metadata 中的 Site 命名（BO → Gut）
    Site_16S = case_when(Site == "BO" ~ "Gut", TRUE ~ Site)
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("WT_Control", "WT_Treated", "DB_Control", "DB_Treated")
    )
  )

cat("\n样本分组统计：\n")
print(table(rna_meta$Site, rna_meta$Group))

write.csv(
  rna_meta,
  file.path(OUT_DIR, "Integration", "RNA_sample_metadata.csv"),
  row.names = FALSE
)

# ────────────────────────────────────────────────────────────
# 5. 全局配色方案
# ────────────────────────────────────────────────────────────

GROUP_COLORS <- c(
  "WT_Control" = "#4DBBD5",
  "WT_Treated" = "#00A087",
  "DB_Control" = "#E64B35",
  "DB_Treated" = "#F39B7F"
)

# ────────────────────────────────────────────────────────────
# 6. 加载 MSigDB 基因集（小鼠，用于 GSEA）
# ────────────────────────────────────────────────────────────

message("\n加载 MSigDB 基因集（Mus musculus）...")

.load_msig <- function(cat, subcat = NULL) {
  df <- if (is.null(subcat)) {
    msigdbr(species = "Mus musculus", category = cat)
  } else {
    msigdbr(species = "Mus musculus", category = cat, subcategory = subcat)
  }
  split(df$gene_symbol, df$gs_name)
}

msig_hallmark <- .load_msig("H")
msig_go_bp <- .load_msig("C5", "GO:BP")
msig_kegg_gsea <- tryCatch(
  .load_msig("C2", "CP:KEGG_LEGACY"),
  error = function(e) {
    message("  CP:KEGG_LEGACY 不存在，尝试 CP:KEGG ...")
    .load_msig("C2", "CP:KEGG")
  }
)

message(sprintf(
  "  Hallmark: %d | GO BP: %d | KEGG: %d",
  length(msig_hallmark),
  length(msig_go_bp),
  length(msig_kegg_gsea)
))

# 重点关注的 pathway 关键词（大写，匹配 MSigDB 命名）
# 涵盖：炎症反应、上皮分化/角化、屏障功能、细胞因子信号、氧化应激/创伤愈合
FOCUS_RE <- paste(
  c(
    "INFLAMMATORY",
    "INTERLEUKIN",
    "IL[0-9]",
    "TNF",
    "NFKB",
    "NF_KB",
    "CYTOKINE",
    "INTERFERON",
    "CHEMOKINE",
    "JAK_STAT",
    "KERATINIZATION",
    "KERATINOCYTE",
    "CORNIFICATION",
    "EPITHELIAL.*DIFF",
    "SKIN_DEVELOPMENT",
    "BARRIER",
    "TIGHT_JUNCTION",
    "OXIDATIVE_STRESS",
    "REACTIVE_OXYGEN",
    "WOUND"
  ),
  collapse = "|"
)

# ────────────────────────────────────────────────────────────
# 7. 辅助函数
# ────────────────────────────────────────────────────────────

# ── 7.1 Volcano 图 ──────────────────────────────────────────
.make_volcano <- function(res_df, title) {
  # 需要 gene_name 列；若无则用 gene_id
  if (!"gene_name" %in% colnames(res_df)) {
    res_df$gene_name <- rownames(res_df)
  }

  df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      sig = case_when(
        padj < P$padj_cutoff & log2FoldChange > P$lfc_cutoff ~ "Up",
        padj < P$padj_cutoff & log2FoldChange < -P$lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(
        sig != "NS" & abs(log2FoldChange) > 2 & padj < 0.01,
        gene_name,
        NA_character_
      )
    )

  n_up <- sum(df$sig == "Up", na.rm = TRUE)
  n_down <- sum(df$sig == "Down", na.rm = TRUE)

  ggplot(df, aes(log2FoldChange, -log10(padj + 1e-300), color = sig)) +
    geom_point(alpha = 0.45, size = 1.2) +
    geom_text_repel(
      aes(label = label),
      size = 2.8,
      max.overlaps = 20,
      na.rm = TRUE,
      segment.color = "grey60"
    ) +
    scale_color_manual(
      values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey75"),
      labels = c(
        Up = paste0("Up (n=", n_up, ")"),
        Down = paste0("Down (n=", n_down, ")"),
        NS = "NS"
      )
    ) +
    geom_vline(
      xintercept = c(-P$lfc_cutoff, P$lfc_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(P$padj_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey50"
    ) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(padj)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(size = 11))
}

# ── 7.2 单次 fgsea + 可视化 ─────────────────────────────────
.run_fgsea <- function(
  ranked_genes,
  gset_list,
  gset_label,
  contrast_label,
  site_label,
  out_gsea_dir
) {
  message(sprintf(
    "    [GSEA %s | %s]: %d genes",
    gset_label,
    contrast_label,
    length(ranked_genes)
  ))

  res <- tryCatch(
    fgsea(
      pathways = gset_list,
      stats = ranked_genes,
      minSize = P$gsea_min_size,
      maxSize = P$gsea_max_size,
      nPermSimple = P$gsea_nperm,
      eps = 0
    ),
    error = function(e) {
      message("      fgsea 出错: ", e$message)
      NULL
    }
  )

  if (is.null(res) || nrow(res) == 0) {
    return(invisible(NULL))
  }

  res <- res %>%
    arrange(pval) %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "/"))

  # 保存全部结果
  fname <- file.path(
    out_gsea_dir,
    paste0("GSEA_", gset_label, "_", contrast_label)
  )
  write.csv(res, paste0(fname, ".csv"), row.names = FALSE)

  # 显著结果
  res_sig <- res %>%
    filter(!is.na(padj), padj < P$gsea_padj) %>%
    mutate(is_focus = grepl(FOCUS_RE, pathway, ignore.case = TRUE))

  message(sprintf("      显著通路 (padj<%.2f): %d", P$gsea_padj, nrow(res_sig)))
  if (nrow(res_sig) == 0) {
    return(invisible(res))
  }

  # NES 条形图（top 30 by |NES|）
  top_n <- min(30, nrow(res_sig))
  plot_df <- res_sig %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n) %>%
    arrange(NES) %>%
    mutate(
      direction = factor(
        ifelse(NES > 0, "Activated", "Suppressed"),
        levels = c("Activated", "Suppressed")
      ),
      # 去掉固定前缀，重点通路末尾加 ★
      clean_name = str_remove(
        pathway,
        paste0("^(HALLMARK|KEGG_LEGACY|KEGG|GOBP|GOCC|GOMF|GO)_")
      ),
      pathway_label = ifelse(
        is_focus,
        paste0(str_trunc(clean_name, 52), " \u2605"),
        str_trunc(clean_name, 55)
      ),
      # 用于透明度映射
      focus_chr = ifelse(is_focus, "yes", "no")
    )

  p_bar <- ggplot(
    plot_df,
    aes(NES, reorder(pathway_label, NES), fill = direction, alpha = focus_chr)
  ) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c(Activated = "#E64B35", Suppressed = "#4DBBD5")
    ) +
    scale_alpha_manual(values = c("yes" = 1.0, "no" = 0.55), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    labs(
      title = paste0(site_label, "  |  GSEA \u2013 ", gset_label),
      subtitle = paste0(
        contrast_label,
        "   \u2605 = 重点关注通路   top ",
        top_n,
        " by |NES|"
      ),
      x = "NES (Normalized Enrichment Score)",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  ph <- max(5, top_n * 0.32 + 2.5)
  ggsave(paste0(fname, "_barplot.pdf"), p_bar, width = 11, height = ph)
  ggsave(
    paste0(fname, "_barplot.png"),
    p_bar,
    width = 11,
    height = ph,
    dpi = 180
  )

  # 重点 pathway 富集曲线（最多 6 条）
  focus_paths <- res_sig %>%
    filter(is_focus) %>%
    arrange(padj) %>%
    slice_head(n = 6) %>%
    pull(pathway)

  if (length(focus_paths) > 0) {
    pdf(paste0(fname, "_focus_enrichcurve.pdf"), width = 8, height = 5)
    for (fp in focus_paths) {
      tryCatch(
        {
          row_fp <- res[res$pathway == fp, ]
          p_c <- plotEnrichment(gset_list[[fp]], ranked_genes) +
            labs(
              title = str_trunc(fp, 65),
              subtitle = sprintf(
                "NES=%.3f  padj=%.4f  size=%d",
                row_fp$NES,
                row_fp$padj,
                row_fp$size
              )
            ) +
            theme_bw(base_size = 11)
          print(p_c)
        },
        error = function(e) NULL
      )
    }
    dev.off()
  }

  return(invisible(res))
}

# ── 7.3 EntrezID 转换辅助 ───────────────────────────────────
.sym2ez <- function(symbols) {
  tryCatch(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    ),
    error = function(e) {
      message("      bitr 转换出错: ", e$message)
      data.frame(SYMBOL = character(), ENTREZID = character())
    }
  )
}

# ────────────────────────────────────────────────────────────
# 8. 主分析函数（按部位调用）
# ────────────────────────────────────────────────────────────

analyze_site <- function(site_name) {
  cat(sprintf(
    "\n%s\n  \u25b6 分析部位: %s\n%s\n",
    strrep("\u2500", 55),
    site_name,
    strrep("\u2500", 55)
  ))

  out_s <- file.path(OUT_DIR, site_name)
  gsea_dir <- file.path(out_s, "GSEA")
  int_dir <- file.path(OUT_DIR, "Integration")

  # ── 8.1 数据子集 ────────────────────────────────────────
  meta_s <- rna_meta %>%
    filter(Site == site_name) %>%
    arrange(Group)
  rownames(meta_s) <- meta_s$rna_sample

  cnt_s <- count_matrix[, meta_s$rna_sample, drop = FALSE]

  # 仅保留 protein_coding 基因
  if ("gene_biotype" %in% colnames(gene_anno)) {
    pc_ids <- gene_anno$gene_id[
      !is.na(gene_anno$gene_biotype) &
        gene_anno$gene_biotype == "protein_coding"
    ]
    cnt_s <- cnt_s[rownames(cnt_s) %in% pc_ids, , drop = FALSE]
  }

  # 低表达过滤：至少 min_samples 个样本 count ≥ min_count
  cnt_s <- cnt_s[rowSums(cnt_s >= P$min_count) >= P$min_samples, , drop = FALSE]
  # data.frame → matrix，再转 integer（DESeq2 要求整数矩阵）
  cnt_s <- round(as.matrix(cnt_s))
  storage.mode(cnt_s) <- "integer"

  message(sprintf("  蛋白编码基因（过滤后）: %d", nrow(cnt_s)))

  # ── 8.2 DESeq2 对象 & VST 归一化 ────────────────────────
  dds <- DESeqDataSetFromMatrix(
    countData = cnt_s,
    colData = meta_s,
    design = ~Group
  )
  dds$Group <- relevel(dds$Group, ref = "WT_Control")

  # VST（blind=TRUE 用于无监督探索，不受分组影响）
  vsd <- vst(dds, blind = TRUE)
  vst_mat <- assay(vsd)

  # ── 8.3 PCA ─────────────────────────────────────────────
  message("  \u2192 PCA")
  pca_df <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
  pct_var <- round(100 * attr(pca_df, "percentVar"), 1)

  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4.5, alpha = 0.85) +
    # geom_text_repel(
    #   aes(label = name),
    #   size = 3,
    #   max.overlaps = 20,
    #   segment.color = "grey60"
    # ) +
    scale_color_manual(values = GROUP_COLORS) +
    labs(
      title = paste0(site_name, "  \u2014  PCA（VST 归一化）"),
      x = paste0("PC1 (", pct_var[1], "%)"),
      y = paste0("PC2 (", pct_var[2], "%)")
    ) +
    theme_bw(base_size = 13) +
    theme(legend.title = element_blank())

  ggsave(file.path(out_s, "PCA.pdf"), p_pca, width = 8, height = 6)
  ggsave(file.path(out_s, "PCA.png"), p_pca, width = 8, height = 6, dpi = 200)

  # ── 8.4 样本距离热图 ─────────────────────────────────────
  message("  \u2192 样本距离热图")
  dist_mat <- as.matrix(dist(t(vst_mat)))
  anno_col <- data.frame(Group = meta_s$Group, row.names = meta_s$rna_sample)

  pdf(file.path(out_s, "Sample_Distance_Heatmap.pdf"), width = 9, height = 8)
  pheatmap(
    dist_mat,
    annotation_col = anno_col,
    annotation_colors = list(Group = GROUP_COLORS),
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
    main = paste0(site_name, "  \u2014  样本间欧氏距离（VST）"),
    fontsize = 10
  )
  dev.off()

  # ── 8.5 DESeq2 Wald test（两两对比，获取 FC 和 stat） ───
  message("  \u2192 DESeq2 Wald test（两两对比）")
  dds_wald <- DESeq(dds, quiet = TRUE)

  # 四个关键对比
  #   DB_vs_WT        : 糖尿病效应（不治疗时 DB vs WT）
  #   DBTreated_vs_DB : 治疗效应（在 DB 小鼠中，药物的作用）
  #   DBTreated_vs_WT : 治疗后与正常对照的差距是否消除
  #   WTTreated_vs_WT : 药物本身对正常小鼠的效应（安全性对照）
  contrasts_def <- list(
    DB_vs_WT = c("Group", "DB_Control", "WT_Control"),
    DBTreated_vs_DB = c("Group", "DB_Treated", "DB_Control"),
    DBTreated_vs_WT = c("Group", "DB_Treated", "WT_Control"),
    WTTreated_vs_WT = c("Group", "WT_Treated", "WT_Control")
  )

  contrast_res <- lapply(names(contrasts_def), function(cn) {
    res <- results(
      dds_wald,
      contrast = contrasts_def[[cn]],
      alpha = P$padj_cutoff
    ) %>%
      as.data.frame() %>%
      mutate(gene_id = rownames(.)) %>%
      left_join(
        gene_anno %>% dplyr::select(gene_id, gene_name) %>% distinct(),
        by = "gene_id"
      )
    write.csv(
      res,
      file.path(out_s, paste0("DEG_", cn, ".csv")),
      row.names = FALSE
    )
    res
  })
  names(contrast_res) <- names(contrasts_def)

  # 打印各对比 DEG 数量汇总
  for (cn in names(contrast_res)) {
    n_up <- sum(
      contrast_res[[cn]]$padj < P$padj_cutoff &
        contrast_res[[cn]]$log2FoldChange > P$lfc_cutoff,
      na.rm = TRUE
    )
    n_down <- sum(
      contrast_res[[cn]]$padj < P$padj_cutoff &
        contrast_res[[cn]]$log2FoldChange < -P$lfc_cutoff,
      na.rm = TRUE
    )
    message(sprintf("    %-22s  Up=%d  Down=%d", cn, n_up, n_down))
  }

  # ── 8.6 DESeq2 LRT（多组整体检验）────────────────────────
  # LRT 检验"分组变量整体上是否对基因表达有影响"
  # 不预设哪两组比，是多组设计下识别组间差异基因最严谨的方式
  message("  \u2192 DESeq2 LRT（多组整体检验）")
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1, quiet = TRUE)

  lrt_all <- results(dds_lrt) %>%
    as.data.frame() %>%
    mutate(gene_id = rownames(.)) %>%
    left_join(
      gene_anno %>% dplyr::select(gene_id, gene_name) %>% distinct(),
      by = "gene_id"
    ) %>%
    arrange(padj)

  write.csv(lrt_all, file.path(out_s, "LRT_results_all.csv"), row.names = FALSE)

  lrt_sig <- lrt_all %>% filter(!is.na(padj), padj < P$padj_cutoff)
  write.csv(lrt_sig, file.path(out_s, "LRT_results_sig.csv"), row.names = FALSE)
  message(sprintf(
    "  LRT 显著基因 (padj<%.2f): %d",
    P$padj_cutoff,
    nrow(lrt_sig)
  ))

  # ── 8.7 Volcano 图 ───────────────────────────────────────
  message("  \u2192 Volcano 图")
  for (cn in c("DB_vs_WT", "DBTreated_vs_DB")) {
    pv <- .make_volcano(contrast_res[[cn]], paste0(site_name, "  |  ", cn))
    ggsave(
      file.path(out_s, paste0("Volcano_", cn, ".pdf")),
      pv,
      width = 7,
      height = 6
    )
    ggsave(
      file.path(out_s, paste0("Volcano_", cn, ".png")),
      pv,
      width = 7,
      height = 6,
      dpi = 200
    )
  }

  # ── 8.8 HVG 热图（LRT 显著基因 Z-score 热图）────────────
  message("  \u2192 HVG 热图")

  top_ids <- lrt_sig %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    head(P$heatmap_n) %>%
    pull(gene_id)

  # 若 LRT 显著基因 < 20，改用方差最高基因补充
  if (length(top_ids) < 20) {
    message("    LRT 显著基因较少，改用方差最高基因")
    gv <- sort(apply(vst_mat, 1, var), decreasing = TRUE)
    top_ids <- names(gv)[seq_len(min(P$heatmap_n, length(gv)))]
  }

  mat_z <- t(scale(t(vst_mat[top_ids, , drop = FALSE])))
  gn <- gene_anno$gene_name[match(rownames(mat_z), gene_anno$gene_id)]
  rownames(mat_z) <- ifelse(is.na(gn) | gn == "", rownames(mat_z), gn)

  # 按组排序样本列
  s_order <- meta_s %>% arrange(Group) %>% pull(rna_sample)
  anno_heat <- data.frame(Group = meta_s[s_order, "Group"], row.names = s_order)

  pdf(file.path(out_s, "HVG_Heatmap_top60.pdf"), width = 12, height = 14)
  pheatmap(
    mat_z[, s_order],
    annotation_col = anno_heat,
    annotation_colors = list(Group = GROUP_COLORS),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 8,
    fontsize_col = 9,
    main = paste0(
      site_name,
      "  \u2014  Top LRT 基因 Z-score 热图（n=",
      length(top_ids),
      "）"
    )
  )
  dev.off()

  # ── 8.9 GO BP 富集分析 ───────────────────────────────────
  message("  \u2192 GO BP 富集分析")

  sig_sym <- lrt_sig %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    pull(gene_name) %>%
    unique()
  bg_sym <- gene_anno %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    pull(gene_name) %>%
    unique()

  sig_ez <- .sym2ez(sig_sym)
  bg_ez <- .sym2ez(bg_sym)

  if (nrow(sig_ez) > 0) {
    go_res <- enrichGO(
      gene = sig_ez$ENTREZID,
      universe = if (nrow(bg_ez) > 0) bg_ez$ENTREZID else NULL,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      write.csv(
        go_res@result,
        file.path(out_s, "GO_BP_results.csv"),
        row.names = FALSE
      )

      p_go <- dotplot(
        go_res,
        showCategory = 20,
        title = paste0(site_name, "  \u2014  GO BP（LRT 显著基因）")
      )
      ggsave(
        file.path(out_s, "GO_BP_dotplot.pdf"),
        p_go,
        width = 10,
        height = 8
      )
      ggsave(
        file.path(out_s, "GO_BP_dotplot.png"),
        p_go,
        width = 10,
        height = 8,
        dpi = 200
      )

      # GO term 间相似性网络图
      tryCatch(
        {
          p_emap <- emapplot(
            pairwise_termsim(go_res),
            showCategory = 30,
            layout = "nicely"
          )
          ggsave(
            file.path(out_s, "GO_BP_emap.pdf"),
            p_emap,
            width = 12,
            height = 10
          )
        },
        error = function(e) message("    GO emap 跳过: ", e$message)
      )

      message(sprintf(
        "    GO BP 显著 terms: %d",
        sum(go_res@result$p.adjust < 0.05)
      ))
    } else {
      message("    GO BP: 无显著结果")
    }
  }

  # ── 8.10 KEGG 富集分析（需要网络连接）─────────────────────
  message("  \u2192 KEGG 富集分析（需网络）")

  if (nrow(sig_ez) > 0) {
    kegg_res <- tryCatch(
      enrichKEGG(
        gene = sig_ez$ENTREZID,
        organism = "mmu",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        universe = if (nrow(bg_ez) > 0) bg_ez$ENTREZID else NULL
      ),
      error = function(e) {
        message("    KEGG 连接失败（可检查网络或使用离线缓存）: ", e$message)
        NULL
      }
    )

    if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
      kegg_res <- setReadable(
        kegg_res,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID"
      )
      write.csv(
        kegg_res@result,
        file.path(out_s, "KEGG_results.csv"),
        row.names = FALSE
      )

      p_kegg <- dotplot(
        kegg_res,
        showCategory = 20,
        title = paste0(site_name, "  \u2014  KEGG（LRT 显著基因）")
      )
      ggsave(
        file.path(out_s, "KEGG_dotplot.pdf"),
        p_kegg,
        width = 10,
        height = 8
      )
      ggsave(
        file.path(out_s, "KEGG_dotplot.png"),
        p_kegg,
        width = 10,
        height = 8,
        dpi = 200
      )

      message(sprintf(
        "    KEGG 显著通路: %d",
        sum(kegg_res@result$p.adjust < 0.05)
      ))
    } else {
      message("    KEGG: 无显著结果")
    }
  }

  # ── 8.11 GSEA 分析 ───────────────────────────────────────
  message("  \u2192 GSEA 分析（Hallmark + GO BP + KEGG）")

  # 对 3 个关键对比运行 GSEA
  # ranked_genes 基于 Wald stat 排序（含正负方向信息）
  for (cn in c("DB_vs_WT", "DBTreated_vs_DB", "WTTreated_vs_WT")) {
    ranked <- contrast_res[[cn]] %>%
      filter(!is.na(stat), !is.na(gene_name), gene_name != "") %>%
      group_by(gene_name) %>%
      slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(stat)) %>%
      {
        setNames(.$stat, .$gene_name)
      }

    .run_fgsea(ranked, msig_hallmark, "Hallmark", cn, site_name, gsea_dir)
    .run_fgsea(ranked, msig_go_bp, "GOBP", cn, site_name, gsea_dir)
    .run_fgsea(ranked, msig_kegg_gsea, "KEGG", cn, site_name, gsea_dir)
  }

  # ── 8.12 导出 16S 联合分析接口数据 ──────────────────────
  message("  \u2192 导出整合接口数据")

  # VST 表达矩阵（含基因注释）
  vst_df <- as.data.frame(vst_mat) %>%
    mutate(gene_id = rownames(.)) %>%
    left_join(
      gene_anno %>% dplyr::select(gene_id, gene_name) %>% distinct(),
      by = "gene_id"
    ) %>%
    dplyr::select(gene_id, gene_name, everything())

  write.csv(
    vst_df,
    file.path(int_dir, paste0(site_name, "_VST_matrix.csv")),
    row.names = FALSE
  )
  write.csv(
    meta_s,
    file.path(int_dir, paste0(site_name, "_metadata.csv")),
    row.names = FALSE
  )
  write.csv(
    lrt_sig,
    file.path(int_dir, paste0(site_name, "_LRT_sig.csv")),
    row.names = FALSE
  )

  # DESeq2 对象（方便后续精细分析复用）
  saveRDS(dds_wald, file.path(int_dir, paste0(site_name, "_dds.rds")))
  saveRDS(vsd, file.path(int_dir, paste0(site_name, "_vsd.rds")))

  message(sprintf("  \u2713 %s 分析完成\n", site_name))

  # 返回关键对象供后续使用
  list(
    meta_s = meta_s,
    dds_wald = dds_wald,
    vsd = vsd,
    lrt_all = lrt_all,
    lrt_sig = lrt_sig,
    contrast_res = contrast_res
  )
}

# ────────────────────────────────────────────────────────────
# 9. 逐部位执行分析
# ────────────────────────────────────────────────────────────

all_results <- list()
for (site in c("Lung", "BO", "Oral")) {
  all_results[[site]] <- tryCatch(
    analyze_site(site),
    error = function(e) {
      message(sprintf("!!! 部位 %s 分析失败: %s", site, e$message))
      NULL
    }
  )
}

# ────────────────────────────────────────────────────────────
# 10. 汇总报告
# ────────────────────────────────────────────────────────────

cat(sprintf("\n%s\n   分析完成 — 汇总\n%s\n", strrep("=", 45), strrep("=", 45)))

summary_df <- data.frame(
  Site = names(all_results),
  LRT_sig = sapply(all_results, function(x) {
    if (!is.null(x)) nrow(x$lrt_sig) else NA
  }),
  stringsAsFactors = FALSE
)
print(summary_df)
write.csv(
  summary_df,
  file.path(OUT_DIR, "Summary_LRT_sig.csv"),
  row.names = FALSE
)

cat("\n输出目录:", OUT_DIR, "\n")

# ────────────────────────────────────────────────────────────
# 11. 写出 16S 联合分析说明
# ────────────────────────────────────────────────────────────

guide_text <- '
# ============================================================
# 转录组 x 16S 扩增子 联合分析接口说明
# ============================================================

## 部位对应关系

  RNA 部位   <->  16S 部位
  --------------------
  Lung       <->  Lung  (样本后缀 _L)
  BO         <->  Gut   (样本后缀 _B)
  Oral       <->  Oral  (样本后缀 _O)
  （RNA 无 Fecal 数据，16S 的 Fecal 无法与 RNA 直接匹配）

## 样本 ID 对应规则（RNA <-> 16S）

  RNA 样本名    16S 样本名
  dbN1_Lung  <-> dbN_1_L
  dbS2_BO    <-> dbS_2_B
  WTN3_Oral  <-> WTN_3_O
  WTS1_Lung  <-> WTS_1_L
  规则：RNA 去掉部位后缀 → 在数字前插入下划线 → 加部位缩写后缀

## Integration/ 目录文件说明

  *_VST_matrix.csv      VST 表达矩阵（行=基因，列=样本，含 gene_name）
  *_metadata.csv        样本信息（含 Site_16S 字段对应 16S 命名）
  *_LRT_sig.csv         LRT 多组显著基因列表（推荐作为整合分析 gene universe）
  *_dds.rds             DESeq2 对象（可提取 normalized counts 等）
  *_vsd.rds             VSD 对象（可直接提取 VST 矩阵）
  RNA_sample_metadata.csv  全部 RNA 样本元信息汇总

## 推荐联合分析方案

### 方案 A：直接 Spearman 相关（探索性，最简单）
  # 取 LRT 显著基因的 VST 矩阵 与 16S OTU 丰度表做相关
  # 注意先对齐样本顺序

  library(Hmisc)
  vst_lung <- read.csv("Integration/Lung_VST_matrix.csv")
  # expr_mat：行=基因，列=样本（取 LRT 显著基因）
  # otu_mat ：行=OTU，列=样本（来自 16S Gut 分析）
  cor_res <- rcorr(t(expr_mat), t(otu_mat), type = "spearman")

### 方案 B：MaAsLin2 多变量关联（校正混杂因素）
  library(Maaslin2)
  # 将 RNA 表达（或模块特征向量）作为固定效应协变量
  # 对 16S OTU/ASV 丰度做线性混合模型关联

### 方案 C：WGCNA 模块 x 16S 相关（最推荐，系统且可解释）

  library(WGCNA)

  # Step 1：构建 RNA 共表达网络（以 Lung 为例）
  vst_df   <- read.csv("Integration/Lung_VST_matrix.csv")
  expr_mat <- as.matrix(vst_df[, -(1:2)])   # 去掉 gene_id/gene_name
  rownames(expr_mat) <- vst_df$gene_name
  expr_t   <- t(expr_mat)                   # 转置：样本 x 基因

  # 选择软阈值
  sft <- pickSoftThreshold(expr_t, powerVector = 1:20,
                           networkType = "unsigned")
  power_use <- sft$powerEstimate

  # 构建网络
  net <- blockwiseModules(
    expr_t,
    power          = power_use,
    TOMType        = "unsigned",
    minModuleSize  = 30,
    mergeCutHeight = 0.25,
    verbose        = 3
  )

  # Step 2：提取模块特征向量（eigengene）
  MEs <- net$MEs  # 样本 x 模块（每个模块的第一主成分）

  # Step 3：与 16S OTU 丰度做相关
  # otu_gut <- read.csv("path/to/16S_Gut_OTU.csv")
  # 样本对齐后：
  # cor_mat <- cor(MEs, otu_gut, method = "spearman")
  # pheatmap(cor_mat) 可视化模块-菌群关联

## 注意事项
  - 联合分析建议在同一部位内进行（Lung-Lung，BO-Gut，Oral-Oral）
  - 样本量较小（每组 n=3），相关性分析结论需谨慎解读
  - 建议用 FDR 校正多重比较（p.adjust 或 qvalue）
'

writeLines(
  guide_text,
  file.path(OUT_DIR, "Integration", "00_Integration_Guide.txt")
)

cat("整合说明: Integration/00_Integration_Guide.txt\n")
cat("\n全部完成。\n")
