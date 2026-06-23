# ================================================================
# 转录组 × 16S 联合分析脚本
# 假设：糖尿病状态通过影响宿主基因表达进而影响微生物群落
#
# 核心流程：
#   1. 从挽救簇基因计算 Eigengene（PC1）—— 2个方向各1个
#   2. 识别 16S 菌群中的挽救菌属（DB→DBTreated方向逆转）
#   3. Eigengene × 菌属丰度 Spearman 相关热图
#   4. 整合趋势折线图（宿主 + 菌群共同恢复可视化）
#
# 前提：已运行 RNA_transcriptome_analysis.R 和 RNA_rescue_analysis.R
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(DESeq2)
  library(Hmisc) # rcorr（Spearman相关+p值）
  library(stringr)
  library(patchwork) # 图形拼接
})

set.seed(42)

# ────────────────────────────────────────────────────────────
# 0. 路径
# ────────────────────────────────────────────────────────────

BASE <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
RNA_DIR <- file.path(BASE, "analysis_results/RNA")
S16_DIR <- file.path(BASE, "analysis_results/tables")
INT_DIR <- file.path(RNA_DIR, "Integration")
OUT_DIR <- file.path(BASE, "analysis_results/Joint")

for (site in c("Lung", "BO", "Oral")) {
  dir.create(file.path(OUT_DIR, site), recursive = TRUE, showWarnings = FALSE)
}

# ────────────────────────────────────────────────────────────
# 1. 参数
# ────────────────────────────────────────────────────────────

P <- list(
  genus_min_prev = 0.3, # 菌属最低出现频率（30%样本中存在）
  genus_min_abund = 0.001, # 菌属最低平均相对丰度（0.1%）
  cor_padj_cutoff = 0.2, # 相关分析FDR阈值（小样本放宽）
  rescue_padj = 0.05, # 挽救菌属筛选阈值
  rescue_lfc = 0.5 # 挽救菌属最低|LFC|
)

GROUP_COLORS <- c(
  "WT_Control" = "#4DBBD5",
  "WT_Treated" = "#00A087",
  "DB_Control" = "#E64B35",
  "DB_Treated" = "#F39B7F"
)

# 部位后缀映射（RNA部位名 → 16S样本列后缀）
SITE_SUFFIX <- c(Lung = "L", BO = "B", Oral = "O")

# RNA部位名 → 16S metadata中的Site名（用于读取DESeq2结果文件）
SITE_16S_NAME <- c(Lung = "Lung", BO = "Gut", Oral = "Oral")

# ────────────────────────────────────────────────────────────
# 2. 读取全局数据
# ────────────────────────────────────────────────────────────

message("读取属水平丰度表...")
genus_all <- read.csv(
  file.path(S16_DIR, "g_abundance_all.csv"),
  row.names = 1,
  check.names = FALSE
)
# 去除行名（菌属名）中可能残留的引号
rownames(genus_all) <- gsub('"', '', rownames(genus_all))
# 列名形如 dbN.1.B / WTS.2.L，转换为统一格式备用
message(sprintf(
  "  菌属总数: %d  |  16S样本总数: %d",
  nrow(genus_all),
  ncol(genus_all)
))

# 读取16S的metadata（从RNA的metadata推导，因为样本命名规则一致）
rna_meta_all <- read.csv(
  file.path(INT_DIR, "RNA_sample_metadata.csv"),
  stringsAsFactors = FALSE
)

# ────────────────────────────────────────────────────────────
# 2b. 补充运行 16S DESeq2 差异分析（若结果文件不存在）
# microbiome_analysis.R 中因 TreatGroup 命名问题未能生成这些文件
# ────────────────────────────────────────────────────────────

message("检查 16S DESeq2 结果文件...")

# 读取16S metadata（与RNA metadata相同）
meta_16s_all <- read.table(
  file.path(BASE, "metadata.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1
)
meta_16s_all$SampleID_converted <- gsub("_", ".", rownames(meta_16s_all))
meta_16s_all$TreatGroup <- paste(
  meta_16s_all$Condition,
  meta_16s_all$Treatment,
  sep = "_"
)

# 读取 class 水平绝对丰度（与microbiome_analysis.R一致）
# 但这里用属水平绝对丰度做差异分析
genus_abs_path <- file.path(
  BASE,
  "data/02.ASVanalysis/Taxa_abundance/Evenabs/asv_table.g.absolute.xls"
)

# 若属水平绝对丰度不存在，尝试class水平
if (!file.exists(genus_abs_path)) {
  genus_abs_path <- file.path(
    BASE,
    "data/02.ASVanalysis/Taxa_abundance/Evenabs/asv_table.c.absolute.xls"
  )
  message("  属水平绝对丰度不存在，使用纲水平代替")
}

if (file.exists(genus_abs_path)) {
  genus_abs <- read.table(
    genus_abs_path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )
  # 找共同样本
  common_s <- intersect(meta_16s_all$SampleID_converted, colnames(genus_abs))
  genus_abs <- genus_abs[, common_s, drop = FALSE]
  meta_16s_all <- meta_16s_all[
    meta_16s_all$SampleID_converted %in% common_s,
    ,
    drop = FALSE
  ]

  for (site_16s in c("Lung", "Gut", "Oral")) {
    f_db_wt <- file.path(S16_DIR, paste0("DESeq2_", site_16s, "_DB_vs_WT.csv"))
    f_trt_db <- file.path(
      S16_DIR,
      paste0("DESeq2_", site_16s, "_Treatment_vs_Control.csv")
    )

    if (file.exists(f_db_wt) && file.exists(f_trt_db)) {
      message(sprintf("  %s DESeq2文件已存在，跳过", site_16s))
      next
    }

    site_samples <- meta_16s_all %>%
      filter(Site == site_16s) %>%
      pull(SampleID_converted)

    if (length(site_samples) < 6) {
      next
    }

    meta_site <- meta_16s_all[
      meta_16s_all$SampleID_converted %in% site_samples,
      ,
      drop = FALSE
    ]
    rownames(meta_site) <- meta_site$SampleID_converted
    meta_site$TreatGroup <- factor(meta_site$TreatGroup)

    cnt_site <- round(as.matrix(genus_abs[, site_samples, drop = FALSE]))
    cnt_site <- cnt_site[rowSums(cnt_site) >= 10, , drop = FALSE]
    storage.mode(cnt_site) <- "integer"

    message(sprintf(
      "  运行 %s 16S DESeq2 (%d taxa)...",
      site_16s,
      nrow(cnt_site)
    ))

    dds_16s <- tryCatch(
      {
        d <- DESeqDataSetFromMatrix(
          countData = cnt_site,
          colData = meta_site,
          design = ~TreatGroup
        )
        d$TreatGroup <- relevel(d$TreatGroup, ref = "WT_Control")
        DESeq(d, quiet = TRUE)
      },
      error = function(e) {
        message(sprintf("    DESeq2 失败: %s", e$message))
        NULL
      }
    )

    if (is.null(dds_16s)) {
      next
    }

    # DB_vs_WT
    if ("DB_Control" %in% levels(dds_16s$TreatGroup)) {
      res <- results(
        dds_16s,
        contrast = c("TreatGroup", "DB_Control", "WT_Control")
      ) %>%
        as.data.frame() %>%
        mutate(Taxon = rownames(.)) %>%
        arrange(padj)
      write.csv(res, f_db_wt, row.names = FALSE)
      message(sprintf("    %s DB_vs_WT: 写出 %s", site_16s, f_db_wt))
    }

    # DBTreated_vs_DB
    if (
      "DB_Treated" %in%
        levels(dds_16s$TreatGroup) &&
        "DB_Control" %in% levels(dds_16s$TreatGroup)
    ) {
      res2 <- results(
        dds_16s,
        contrast = c("TreatGroup", "DB_Treated", "DB_Control")
      ) %>%
        as.data.frame() %>%
        mutate(Taxon = rownames(.)) %>%
        arrange(padj)
      write.csv(res2, f_trt_db, row.names = FALSE)
      message(sprintf(
        "    %s Treatment_vs_Control: 写出 %s",
        site_16s,
        f_trt_db
      ))
    }
  }
} else {
  message("  未找到16S绝对丰度文件，跳过补充DESeq2分析")
  message(sprintf("  预期路径: %s", genus_abs_path))
}

# ────────────────────────────────────────────────────────────
# 3. 辅助函数
# ────────────────────────────────────────────────────────────

# ── 3.1 RNA样本名 → 16S样本名转换 ────────────────────────
# RNA: dbN1_Lung → 16S: dbN.1.L
rna_to_16s <- function(rna_name, site_suffix) {
  # 去掉部位后缀
  prefix <- str_remove(rna_name, paste0("_(Lung|BO|Oral)$"))
  # dbN1 → dbN.1，WTS3 → WTS.3（在数字前插入点）
  prefix_dot <- str_replace(prefix, "([A-Za-z]+)(\\d+)$", "\\1.\\2")
  paste0(prefix_dot, ".", site_suffix)
}

# ── 3.2 计算 Eigengene（PC1）──────────────────────────────
compute_eigengene <- function(expr_mat, gene_ids, label) {
  # expr_mat: 基因 × 样本
  ids_use <- intersect(gene_ids, rownames(expr_mat))
  if (length(ids_use) < 3) {
    message(sprintf("    [%s] 可用基因 < 3，跳过", label))
    return(NULL)
  }
  sub_mat <- expr_mat[ids_use, , drop = FALSE]
  # 去除方差为0的基因
  sub_mat <- sub_mat[apply(sub_mat, 1, var) > 0, , drop = FALSE]
  if (nrow(sub_mat) < 3) {
    return(NULL)
  }

  pca <- prcomp(t(sub_mat), scale. = TRUE, center = TRUE)
  pct <- round(100 * pca$sdev[1]^2 / sum(pca$sdev^2), 1)
  eg <- pca$x[, 1]

  message(sprintf(
    "    [%s] %d基因 → Eigengene PC1解释%.1f%%方差",
    label,
    nrow(sub_mat),
    pct
  ))
  list(eigengene = eg, pct_var = pct, n_genes = nrow(sub_mat))
}

# ── 3.3 CLR转换（组成数据相关前标准化）─────────────────────
clr_transform <- function(mat) {
  # mat: 菌属 × 样本，值为相对丰度（0~1）
  # 加伪计数后做CLR
  mat_pseudo <- mat + 1e-6
  log_mat <- log(mat_pseudo)
  # 按列减去几何均值
  sweep(log_mat, 2, colMeans(log_mat), "-")
}

# ── 3.4 相关热图 ──────────────────────────────────────────
plot_cor_heatmap <- function(
  cor_mat,
  pval_mat,
  title,
  sig_mark = TRUE,
  out_path
) {
  # 标记显著性
  sig_mat <- matrix(
    "",
    nrow = nrow(pval_mat),
    ncol = ncol(pval_mat),
    dimnames = dimnames(pval_mat)
  )
  if (sig_mark) {
    # FDR校正（BH）
    padj_vec <- p.adjust(as.vector(pval_mat), method = "BH")
    padj_mat <- matrix(
      padj_vec,
      nrow = nrow(pval_mat),
      dimnames = dimnames(pval_mat)
    )
    sig_mat[padj_mat < 0.2] <- "."
    sig_mat[padj_mat < 0.1] <- "*"
    sig_mat[padj_mat < 0.05] <- "**"
  }

  # 清理菌属名（去掉引号和前缀）
  clean_genus <- function(x) {
    x <- str_remove_all(x, '"')
    x <- str_replace(x, "^[kpcofgs]__", "")
    x
  }
  rownames(cor_mat) <- clean_genus(rownames(cor_mat))
  rownames(sig_mat) <- rownames(cor_mat)

  # 截断行名避免过长
  rownames(cor_mat) <- str_trunc(rownames(cor_mat), 35)
  rownames(sig_mat) <- rownames(cor_mat)

  ph <- max(4, nrow(cor_mat) * 0.28 + 2)
  pw <- max(4, ncol(cor_mat) * 1.2 + 3)

  pdf(paste0(out_path, ".pdf"), width = pw, height = ph)
  pheatmap(
    cor_mat,
    display_numbers = sig_mat,
    number_color = "black",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 9,
    fontsize_col = 11,
    fontsize_number = 10,
    main = title,
    border_color = "grey90"
  )
  dev.off()

  # 同时保存PNG
  png(paste0(out_path, ".png"), width = pw * 100, height = ph * 100, res = 120)
  pheatmap(
    cor_mat,
    display_numbers = sig_mat,
    number_color = "black",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 9,
    fontsize_col = 11,
    fontsize_number = 10,
    main = title,
    border_color = "grey90"
  )
  dev.off()
}

# ── 3.5 整合趋势折线图 ────────────────────────────────────
# 左：Eigengene趋势；右：菌属丰度趋势（仅挽救菌属）
plot_trend_panel <- function(eg_df, genus_trend_df, site_label, out_path) {
  # eg_df: Group, Eigengene_name, value（已是各样本的eigengene值）
  # genus_trend_df: Group, Genus, clr_value

  grp_order <- c("WT_Control", "DB_Control", "DB_Treated")

  # ── 左图：Eigengene按组均值±SE ──
  eg_summary <- eg_df %>%
    filter(Group %in% grp_order) %>%
    mutate(Group = factor(Group, levels = grp_order)) %>%
    group_by(Eigengene, Group) %>%
    summarise(
      mean_val = mean(value),
      se_val = sd(value) / sqrt(n()),
      .groups = "drop"
    )

  p_left <- ggplot(
    eg_summary,
    aes(Group, mean_val, color = Eigengene, group = Eigengene)
  ) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 4) +
    geom_errorbar(
      aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
      width = 0.15,
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.5
    ) +
    scale_color_manual(
      values = c(
        "EG_DB_up" = "#E64B35",
        "EG_DB_down" = "#4DBBD5"
      ),
      labels = c(
        "EG_DB_up" = "DB中上调\n（治疗后回落）",
        "EG_DB_down" = "DB中下调\n（治疗后回升）"
      )
    ) +
    labs(
      title = paste0(site_label, "  宿主 Eigengene"),
      subtitle = "PC1 均值 ± SE",
      x = NULL,
      y = "Eigengene（PC1得分）",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold")
    )

  # ── 右图：挽救菌属 CLR丰度按组均值 ──
  if (!is.null(genus_trend_df) && nrow(genus_trend_df) > 0) {
    genus_summary <- genus_trend_df %>%
      filter(Group %in% grp_order) %>%
      mutate(
        Group = factor(Group, levels = grp_order),
        Genus = str_trunc(str_remove_all(Genus, '"'), 30)
      ) %>%
      group_by(Genus, Group) %>%
      summarise(mean_clr = mean(clr_value), .groups = "drop")

    p_right <- ggplot(
      genus_summary,
      aes(Group, mean_clr, color = Genus, group = Genus)
    ) +
      geom_line(linewidth = 1.0, alpha = 0.85) +
      geom_point(size = 3.0, alpha = 0.85) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey60",
        linewidth = 0.5
      ) +
      labs(
        title = paste0(site_label, "  挽救菌属（属水平）"),
        subtitle = "CLR均值  |  仅显示挽救趋势菌属",
        x = NULL,
        y = "CLR 丰度",
        color = "菌属"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 11, face = "bold")
      )
  } else {
    p_right <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "无显著挽救菌属",
        size = 5,
        color = "grey50"
      ) +
      theme_void()
  }

  p_combined <- p_left +
    p_right +
    plot_layout(widths = c(1, 1.5)) +
    plot_annotation(
      title = paste0(site_label, "  ——  宿主-菌群共同恢复趋势"),
      theme = theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    )

  ggsave(paste0(out_path, ".pdf"), p_combined, width = 14, height = 5.5)
  ggsave(
    paste0(out_path, ".png"),
    p_combined,
    width = 14,
    height = 5.5,
    dpi = 200
  )
}

# ────────────────────────────────────────────────────────────
# 4. 主分析函数（按部位）
# ────────────────────────────────────────────────────────────

analyze_joint <- function(site_name) {
  cat(sprintf(
    "\n%s\n  \u25b6 联合分析: %s\n%s\n",
    strrep("\u2500", 55),
    site_name,
    strrep("\u2500", 55)
  ))

  out_s <- file.path(OUT_DIR, site_name)
  suffix_16s <- SITE_SUFFIX[site_name]

  # ── 4.1 加载 RNA 数据 ─────────────────────────────────────
  vsd <- readRDS(file.path(INT_DIR, paste0(site_name, "_vsd.rds")))
  vst_mat <- assay(vsd)

  meta_rna <- read.csv(
    file.path(INT_DIR, paste0(site_name, "_metadata.csv")),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Group = factor(
        Group,
        levels = c("WT_Control", "WT_Treated", "DB_Control", "DB_Treated")
      )
    )

  cluster_summary <- read.csv(
    file.path(RNA_DIR, "Rescue", site_name, "C_kmeans_cluster_summary.csv"),
    stringsAsFactors = FALSE
  )

  rescue_genes_all <- read.csv(
    file.path(RNA_DIR, "Rescue", site_name, "C_rescue_cluster_genes.csv"),
    stringsAsFactors = FALSE
  )

  # 识别挽救簇（is_rescue=TRUE 且 rescue_mag=TRUE）
  rescue_clusters <- cluster_summary %>%
    filter(is_rescue == TRUE, rescue_mag == TRUE) %>%
    pull(cluster)

  if (length(rescue_clusters) == 0) {
    message("  无挽救簇，跳过")
    return(invisible(NULL))
  }

  message(sprintf("  挽救簇编号: %s", paste(rescue_clusters, collapse = ", ")))

  # 按 disease_dir 合并为2个方向的基因集
  # disease_dir=1  → DB中上调，治疗后下调（"炎症/疾病激活型"）
  # disease_dir=-1 → DB中下调，治疗后上调（"保护性/功能性下调型"）
  dir_info <- cluster_summary %>%
    filter(cluster %in% rescue_clusters) %>%
    dplyr::select(cluster, disease_dir)

  genes_db_up <- rescue_genes_all %>%
    inner_join(dir_info, by = "cluster") %>%
    filter(disease_dir == 1) %>%
    pull(gene_id) %>%
    unique()

  genes_db_down <- rescue_genes_all %>%
    inner_join(dir_info, by = "cluster") %>%
    filter(disease_dir == -1) %>%
    pull(gene_id) %>%
    unique()

  message(sprintf(
    "  DB上调挽救基因: %d  |  DB下调挽救基因: %d",
    length(genes_db_up),
    length(genes_db_down)
  ))

  # ── 4.2 计算 Eigengene ────────────────────────────────────
  message("  \u2192 计算 Eigengene（PC1）")

  eg_up <- compute_eigengene(vst_mat, genes_db_up, "EG_DB_up")
  eg_down <- compute_eigengene(vst_mat, genes_db_down, "EG_DB_down")

  # 构建 Eigengene 数据框
  eg_list <- list()
  if (!is.null(eg_up)) {
    eg_list[["EG_DB_up"]] <- eg_up$eigengene
  }
  if (!is.null(eg_down)) {
    eg_list[["EG_DB_down"]] <- eg_down$eigengene
  }

  if (length(eg_list) == 0) {
    message("  Eigengene 计算失败，跳过")
    return(invisible(NULL))
  }

  # 确保 Eigengene 的方向语义正确：
  # EG_DB_up 应该在 DB_Control 样本中 > WT_Control → 若反了则取负
  for (eg_name in names(eg_list)) {
    eg_vals <- eg_list[[eg_name]]
    db_mean <- mean(eg_vals[meta_rna$rna_sample[
      meta_rna$Group == "DB_Control"
    ]])
    wt_mean <- mean(eg_vals[meta_rna$rna_sample[
      meta_rna$Group == "WT_Control"
    ]])
    expected_up <- eg_name == "EG_DB_up"
    if (expected_up && db_mean < wt_mean) {
      eg_list[[eg_name]] <- -eg_vals
      message(sprintf("    [%s] PC1方向已翻转以符合生物学方向", eg_name))
    }
    if (!expected_up && db_mean > wt_mean) {
      eg_list[[eg_name]] <- -eg_vals
      message(sprintf("    [%s] PC1方向已翻转以符合生物学方向", eg_name))
    }
  }

  # 保存 Eigengene 值
  eg_df_wide <- as.data.frame(eg_list) %>%
    mutate(rna_sample = rownames(.)) %>%
    left_join(
      meta_rna %>% dplyr::select(rna_sample, Group, Replicate),
      by = "rna_sample"
    )

  write.csv(
    eg_df_wide,
    file.path(out_s, "eigengene_values.csv"),
    row.names = FALSE
  )

  # 转长格式，供趋势图使用
  eg_df_long <- eg_df_wide %>%
    pivot_longer(
      cols = starts_with("EG_"),
      names_to = "Eigengene",
      values_to = "value"
    )

  # ── 4.3 处理 16S 属水平数据 ──────────────────────────────
  message("  \u2192 处理 16S 属水平数据")

  # 提取该部位的16S样本列
  # 16S列名格式: dbN.1.B / WTS.2.L 等，按后缀筛选
  cols_16s <- colnames(genus_all)[
    str_ends(colnames(genus_all), paste0("\\.", suffix_16s))
  ]
  # 只取4组（排除Fecal，16S中B/L/O各12样本）
  genus_site <- genus_all[, cols_16s, drop = FALSE]

  # 构建 16S metadata（从列名推导）
  meta_16s <- data.frame(
    s16_sample = cols_16s,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      prefix = str_remove(s16_sample, paste0("\\.", suffix_16s, "$")),
      Condition = case_when(
        str_starts(prefix, "WT") ~ "WT",
        str_starts(prefix, "db") ~ "DB",
        TRUE ~ NA_character_
      ),
      Treatment = case_when(
        str_detect(prefix, "N\\.\\d+$") ~ "Control",
        str_detect(prefix, "S\\.\\d+$") ~ "Treated",
        TRUE ~ NA_character_
      ),
      Group = paste(Condition, Treatment, sep = "_")
    ) %>%
    filter(!is.na(Condition))

  # 按 RNA metadata 的样本顺序建立 16S → RNA 对应关系
  meta_rna <- meta_rna %>%
    mutate(
      s16_sample = rna_to_16s(rna_sample, suffix_16s)
    )

  # 检查匹配
  matched <- intersect(meta_rna$s16_sample, meta_16s$s16_sample)
  message(sprintf(
    "  RNA-16S 样本匹配: %d / %d",
    length(matched),
    nrow(meta_rna)
  ))

  if (length(matched) < 6) {
    message("  匹配样本过少，跳过相关分析")
    return(invisible(NULL))
  }

  # 对齐顺序
  meta_matched <- meta_rna %>%
    filter(s16_sample %in% matched) %>%
    arrange(match(s16_sample, colnames(genus_site)))

  genus_matched <- genus_site[, meta_matched$s16_sample, drop = FALSE]

  # 过滤低丰度 / 低频率菌属
  # 1. 至少在 P$genus_min_prev 比例样本中出现
  prev_ok <- rowSums(genus_matched > 0) / ncol(genus_matched) >=
    P$genus_min_prev
  # 2. 平均丰度 >= 阈值
  abund_ok <- rowMeans(genus_matched) >= P$genus_min_abund
  genus_filt <- genus_matched[prev_ok & abund_ok, , drop = FALSE]

  message(sprintf(
    "  菌属过滤后: %d / %d",
    nrow(genus_filt),
    nrow(genus_matched)
  ))

  if (nrow(genus_filt) < 3) {
    message("  有效菌属 < 3，降低过滤阈值重试")
    genus_filt <- genus_matched[abund_ok, , drop = FALSE]
    message(sprintf("  放宽后菌属: %d", nrow(genus_filt)))
  }

  # CLR 转换
  genus_clr <- clr_transform(genus_filt)
  colnames(genus_clr) <- meta_matched$s16_sample

  # ── 4.4 识别 16S 挽救菌属（Rescue Score，基于组均值方向，不依赖统计阈值）──
  # 与 RNA 挽救分析方法B一致：
  #   RescueScore = -(mean_CLR[DB_Treated] - mean_CLR[DB_Control])
  #                / (mean_CLR[DB_Control] - mean_CLR[WT_Control])
  # Score > 0：治疗向 WT 方向回调（挽救方向正确）
  # Score > 1：完全挽救或过矫正
  # 不设 p 值阈值，仅要求：
  #   1. 疾病效应幅度 |DB - WT| >= min_disease_delta（CLR单位）
  #   2. RescueScore >= min_rescue_score
  message("  \u2192 识别挽救菌属（Rescue Score，基于CLR组均值）")

  # 获取各组样本名
  samp_wt <- meta_matched$s16_sample[meta_matched$Group == "WT_Control"]
  samp_db <- meta_matched$s16_sample[meta_matched$Group == "DB_Control"]
  samp_trt <- meta_matched$s16_sample[meta_matched$Group == "DB_Treated"]

  # 计算各菌属各组 CLR 均值
  genus_score_df <- data.frame(
    Genus = rownames(genus_clr),
    mean_wt = rowMeans(genus_clr[, samp_wt, drop = FALSE]),
    mean_db = rowMeans(genus_clr[, samp_db, drop = FALSE]),
    mean_trt = rowMeans(genus_clr[, samp_trt, drop = FALSE]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      delta_disease = mean_db - mean_wt, # 疾病效应（DB - WT）
      delta_treat = mean_trt - mean_db, # 治疗效应（Treated - DB）
      # Rescue Score：负号使"方向相反"为正值
      rescue_score = ifelse(
        abs(delta_disease) > 1e-6,
        -delta_treat / delta_disease,
        NA_real_
      ),
      # 挽救方向是否正确（治疗效应与疾病效应方向相反）
      is_rescue_dir = sign(delta_treat) != sign(delta_disease)
    ) %>%
    filter(!is.na(rescue_score)) %>%
    arrange(desc(rescue_score))

  write.csv(
    genus_score_df,
    file.path(out_s, "genus_rescue_score_all.csv"),
    row.names = FALSE
  )

  # 筛选挽救菌属：
  #   疾病效应幅度 >= 0.3 CLR（排除本底噪声）
  #   RescueScore >= 0.3（至少挽救30%）
  rescue_genus_df <- genus_score_df %>%
    filter(
      abs(delta_disease) >= 0.3,
      rescue_score >= 0.3,
      is_rescue_dir == TRUE
    )

  write.csv(
    rescue_genus_df,
    file.path(out_s, "rescue_genus_list.csv"),
    row.names = FALSE
  )

  rescue_genus <- rescue_genus_df$Genus
  message(sprintf(
    "  挽救菌属 (RescueScore≥0.3, |ΔDisease|≥0.3): %d 个",
    length(rescue_genus)
  ))

  # 挽救菌属趋势图（象限散点图：x=疾病效应，y=治疗效应，与RNA挽救图逻辑一致）
  p_rs <- ggplot(
    genus_score_df %>% filter(abs(delta_disease) >= 0.1),
    aes(delta_disease, delta_treat)
  ) +
    geom_point(
      data = ~ filter(.x, Genus %in% rescue_genus),
      aes(color = rescue_score),
      size = 3,
      alpha = 0.9
    ) +
    geom_point(
      data = ~ filter(.x, !Genus %in% rescue_genus),
      color = "grey75",
      size = 1.5,
      alpha = 0.4
    ) +
    geom_text_repel(
      data = ~ filter(.x, Genus %in% rescue_genus) %>%
        arrange(desc(rescue_score)) %>%
        head(20),
      aes(label = str_trunc(Genus, 25)),
      size = 2.8,
      max.overlaps = 20,
      segment.color = "grey50"
    ) +
    scale_color_gradientn(
      colors = c("#4DBBD5", "#F39B7F", "#E64B35"),
      name = "Rescue\nScore"
    ) +
    annotate(
      "rect",
      xmin = 0.3,
      xmax = Inf,
      ymin = -Inf,
      ymax = -0.1,
      alpha = 0.05,
      fill = "#E64B35"
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = -0.3,
      ymin = 0.1,
      ymax = Inf,
      alpha = 0.05,
      fill = "#4DBBD5"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    labs(
      title = paste0(site_name, "  菌属挽救趋势象限图"),
      subtitle = "x: DB−WT均值差（CLR）  y: Treated−DB均值差（CLR）\n着色点为挽救菌属（RescueScore≥0.3，|ΔDisease|≥0.3）",
      x = "疾病效应 (mean CLR: DB_Control − WT_Control)",
      y = "治疗效应 (mean CLR: DB_Treated − DB_Control)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8.5))

  ggsave(
    file.path(out_s, "genus_rescue_scatter.pdf"),
    p_rs,
    width = 8,
    height = 7
  )
  ggsave(
    file.path(out_s, "genus_rescue_scatter.png"),
    p_rs,
    width = 8,
    height = 7,
    dpi = 200
  )

  # ── 4.5 Eigengene × 菌属 Spearman 相关 ──────────────────
  message("  \u2192 Spearman 相关分析")

  # Eigengene 矩阵（样本顺序与 genus_clr 一致）
  eg_mat <- do.call(
    cbind,
    lapply(names(eg_list), function(nm) {
      eg_list[[nm]][meta_matched$rna_sample]
    })
  )
  rownames(eg_mat) <- meta_matched$rna_sample
  colnames(eg_mat) <- names(eg_list)

  # genus_clr: 菌属 × 样本 → 转置为 样本 × 菌属
  genus_t <- t(genus_clr)

  # rcorr（Hmisc）：计算 Spearman r 和 p
  combined_mat <- cbind(eg_mat, genus_t)
  cor_full <- rcorr(combined_mat, type = "spearman")

  n_eg <- ncol(eg_mat)
  cor_eg_genus <- cor_full$r[
    1:n_eg,
    (n_eg + 1):ncol(combined_mat),
    drop = FALSE
  ]
  pval_eg_genus <- cor_full$P[
    1:n_eg,
    (n_eg + 1):ncol(combined_mat),
    drop = FALSE
  ]

  # 转置为：菌属 × Eigengene（更符合热图习惯）
  cor_mat_plot <- t(cor_eg_genus)
  pval_mat_plot <- t(pval_eg_genus)

  write.csv(
    cbind(
      as.data.frame(cor_mat_plot),
      as.data.frame(pval_mat_plot) %>%
        setNames(paste0("pval_", colnames(.)))
    ),
    file.path(out_s, "cor_eigengene_genus.csv")
  )

  # ── 热图1：全部过滤后的菌属 ──
  plot_cor_heatmap(
    cor_mat = cor_mat_plot,
    pval_mat = pval_mat_plot,
    title = paste0(
      site_name,
      "  Eigengene × 菌属 Spearman 相关\n",
      "（. p<0.2  * p<0.1  ** p<0.05，BH校正）"
    ),
    out_path = file.path(out_s, "cor_heatmap_all_genus")
  )

  # ── 热图2：仅挽救菌属（若有）──
  rescue_in_mat <- intersect(rescue_genus, rownames(cor_mat_plot))
  if (length(rescue_in_mat) >= 2) {
    plot_cor_heatmap(
      cor_mat = cor_mat_plot[rescue_in_mat, , drop = FALSE],
      pval_mat = pval_mat_plot[rescue_in_mat, , drop = FALSE],
      title = paste0(
        site_name,
        "  Eigengene × 挽救菌属 Spearman 相关\n",
        "（. p<0.2  * p<0.1  ** p<0.05，BH校正）"
      ),
      out_path = file.path(out_s, "cor_heatmap_rescue_genus")
    )
  }

  # ── 4.6 整合趋势折线图 ────────────────────────────────────
  message("  \u2192 整合趋势折线图")

  # 挽救菌属的 CLR 趋势数据
  genus_trend_df <- NULL

  # 优先用显著挽救菌属；若无则用与EG相关最强的前10个菌属
  top_genus <- rescue_in_mat
  if (length(top_genus) == 0) {
    message("  无挽救菌属，改用EG相关最强的前10菌属")
    max_abs_cor <- apply(abs(cor_mat_plot), 1, max)
    top_genus <- names(sort(max_abs_cor, decreasing = TRUE))[
      1:min(10, length(max_abs_cor))
    ]
  }

  if (length(top_genus) > 0) {
    genus_trend_df <- genus_clr[top_genus, , drop = FALSE] %>%
      as.data.frame() %>%
      mutate(Genus = rownames(.)) %>%
      pivot_longer(
        cols = -Genus,
        names_to = "s16_sample",
        values_to = "clr_value"
      ) %>%
      left_join(
        meta_matched %>%
          dplyr::select(s16_sample, Group),
        by = "s16_sample"
      )
  }

  plot_trend_panel(
    eg_df = eg_df_long,
    genus_trend_df = genus_trend_df,
    site_label = site_name,
    out_path = file.path(out_s, "trend_host_microbiome")
  )

  # ── 4.7 显著相关对汇总表 ────────────────────────────────
  sig_pairs <- data.frame(
    Site = site_name,
    Genus = rep(rownames(cor_mat_plot), ncol(cor_mat_plot)),
    Eigengene = rep(colnames(cor_mat_plot), each = nrow(cor_mat_plot)),
    r = as.vector(cor_mat_plot),
    p_raw = as.vector(pval_mat_plot),
    p_adj = p.adjust(as.vector(pval_mat_plot), method = "BH"),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(r)) %>%
    mutate(
      is_rescue_genus = Genus %in% rescue_genus
    ) %>%
    arrange(p_adj, desc(abs(r)))

  write.csv(
    sig_pairs,
    file.path(out_s, "cor_significant_pairs.csv"),
    row.names = FALSE
  )

  n_sig <- sum(sig_pairs$p_adj < P$cor_padj_cutoff, na.rm = TRUE)
  message(sprintf("  显著相关对 (FDR<%.1f): %d", P$cor_padj_cutoff, n_sig))
  message(sprintf("  \u2713 %s 联合分析完成\n", site_name))

  list(
    eg_list = eg_list,
    eg_df_wide = eg_df_wide,
    genus_clr = genus_clr,
    rescue_genus = rescue_genus,
    cor_mat = cor_mat_plot,
    pval_mat = pval_mat_plot,
    sig_pairs = sig_pairs
  )
}

# ────────────────────────────────────────────────────────────
# 5. 逐部位执行
# ────────────────────────────────────────────────────────────

all_joint <- list()
for (site in c("Lung", "BO", "Oral")) {
  all_joint[[site]] <- tryCatch(
    analyze_joint(site),
    error = function(e) {
      message(sprintf("!!! %s 联合分析失败: %s", site, e$message))
      NULL
    }
  )
}

# ────────────────────────────────────────────────────────────
# 6. 跨部位汇总：显著相关对整合表
# ────────────────────────────────────────────────────────────

message("\n汇总跨部位结果...")

all_sig_raw <- bind_rows(lapply(all_joint, function(x) {
  if (!is.null(x) && !is.null(x$sig_pairs)) x$sig_pairs else NULL
}))

all_sig <- if (nrow(all_sig_raw) > 0 && "p_adj" %in% colnames(all_sig_raw)) {
  all_sig_raw %>%
    filter(!is.na(p_adj), p_adj < P$cor_padj_cutoff) %>%
    arrange(Site, p_adj)
} else {
  all_sig_raw
}

write.csv(
  all_sig,
  file.path(OUT_DIR, "All_sites_significant_cor_pairs.csv"),
  row.names = FALSE
)

cat(sprintf("\n%s\n   联合分析完成\n", strrep("=", 50)))
cat(sprintf("   显著相关对 (FDR<%.1f) 汇总:\n", P$cor_padj_cutoff))
if (nrow(all_sig) > 0) {
  print(
    all_sig %>%
      dplyr::select(Site, Genus, Eigengene, r, p_adj, is_rescue_genus) %>%
      head(20)
  )
} else {
  cat(
    "   (暂无达到FDR阈值的显著相关对，可查看各部位 cor_significant_pairs.csv 的原始p值)\n"
  )
}
cat(sprintf("\n输出目录: %s\n", OUT_DIR))
