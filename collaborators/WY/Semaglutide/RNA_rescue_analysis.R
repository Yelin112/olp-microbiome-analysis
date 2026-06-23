# ================================================================
# 挽救/恢复趋势基因分析脚本
# 前提：已成功运行 RNA_transcriptome_analysis.R
# 目标：识别 WT → DB → DB-Treated 呈恢复趋势的基因和通路
# 方法：
#   A. 交集法（方向相反的显著DEG）
#   B. Rescue Score 定量评分 + GSEA
#   C. LRT显著基因的模式聚类（k-means）
# ================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

set.seed(42)

# ────────────────────────────────────────────────────────────
# 0. 路径（与主脚本保持一致）
# ────────────────────────────────────────────────────────────

OUT_DIR <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide/analysis_results/RNA"
INT_DIR <- file.path(OUT_DIR, "Integration")

# 挽救分析输出目录
RSC_DIR <- file.path(OUT_DIR, "Rescue")
for (site in c("Lung", "BO", "Oral")) {
  dir.create(file.path(RSC_DIR, site), recursive = TRUE, showWarnings = FALSE)
}

# ────────────────────────────────────────────────────────────
# 1. 参数
# ────────────────────────────────────────────────────────────

P <- list(
  padj_cutoff    = 0.05,
  lfc_cutoff     = 0.5,    # 挽救基因筛选时的最低 |LFC|（主脚本中为1.0，这里放宽）
  rescue_min_lfc = 0.5,    # DB_vs_WT 的最低 |LFC|（疾病效应需够强才有挽救意义）
  heatmap_n      = 60,
  kmeans_k       = 6,      # 模式聚类的簇数（可调整）
  gsea_padj      = 0.25,
  gsea_nperm     = 10000,
  gsea_min_size  = 15,
  gsea_max_size  = 500
)

GROUP_COLORS <- c(
  "WT_Control"  = "#4DBBD5",
  "WT_Treated"  = "#00A087",
  "DB_Control"  = "#E64B35",
  "DB_Treated"  = "#F39B7F"
)

# 重点通路关键词
FOCUS_RE <- paste(c(
  "INFLAMMATORY", "INTERLEUKIN", "IL[0-9]", "TNF", "NFKB",
  "CYTOKINE", "INTERFERON", "CHEMOKINE", "JAK_STAT",
  "KERATINIZATION", "KERATINOCYTE", "CORNIFICATION",
  "EPITHELIAL", "SKIN_DEVELOPMENT",
  "BARRIER", "TIGHT_JUNCTION",
  "OXIDATIVE_STRESS", "REACTIVE_OXYGEN", "WOUND"
), collapse = "|")

# ────────────────────────────────────────────────────────────
# 2. 加载 MSigDB 基因集
# ────────────────────────────────────────────────────────────

message("加载 MSigDB 基因集...")

.load_msig <- function(cat, subcat = NULL) {
  df <- if (is.null(subcat))
    msigdbr(species = "Mus musculus", category = cat)
  else
    msigdbr(species = "Mus musculus", category = cat, subcategory = subcat)
  split(df$gene_symbol, df$gs_name)
}

msig_hallmark  <- .load_msig("H")
msig_go_bp     <- .load_msig("C5", "GO:BP")
msig_kegg_gsea <- tryCatch(
  .load_msig("C2", "CP:KEGG_LEGACY"),
  error = function(e) .load_msig("C2", "CP:KEGG")
)

# ────────────────────────────────────────────────────────────
# 3. 辅助函数
# ────────────────────────────────────────────────────────────

# ── 3.1 GSEA（专为 Rescue Score 排序设计）──────────────────
.run_rescue_gsea <- function(ranked_genes, gset_list, gset_label,
                             contrast_label, site_label, out_dir) {

  message(sprintf("  [GSEA %s | %s]: %d genes",
                  gset_label, contrast_label, length(ranked_genes)))

  res <- tryCatch(
    fgsea(
      pathways    = gset_list,
      stats       = ranked_genes,
      minSize     = P$gsea_min_size,
      maxSize     = P$gsea_max_size,
      nPermSimple = P$gsea_nperm,
      eps         = 0
    ),
    error = function(e) { message("fgsea 出错: ", e$message); NULL }
  )
  if (is.null(res) || nrow(res) == 0) return(invisible(NULL))

  res <- res %>%
    arrange(pval) %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "/"))

  fname <- file.path(out_dir,
                     paste0("GSEA_Rescue_", gset_label, "_", contrast_label))
  write.csv(res, paste0(fname, ".csv"), row.names = FALSE)

  res_sig <- res %>%
    filter(!is.na(padj), padj < P$gsea_padj) %>%
    mutate(is_focus = grepl(FOCUS_RE, pathway, ignore.case = TRUE))

  if (nrow(res_sig) == 0) {
    message(sprintf("  GSEA %s 无显著通路", gset_label))
    return(invisible(res))
  }

  message(sprintf("  显著通路 (padj<%.2f): %d", P$gsea_padj, nrow(res_sig)))

  # NES 条形图
  top_n   <- min(30, nrow(res_sig))
  plot_df <- res_sig %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n) %>%
    arrange(NES) %>%
    mutate(
      direction = factor(ifelse(NES > 0, "Rescued↑", "Rescued↓"),
                         levels = c("Rescued↑", "Rescued↓")),
      clean_name = stringr::str_remove(
        pathway,
        "^(HALLMARK|KEGG_LEGACY|KEGG|GOBP|GOCC|GOMF|GO)_"
      ),
      pathway_label = ifelse(
        is_focus,
        paste0(stringr::str_trunc(clean_name, 52), " \u2605"),
        stringr::str_trunc(clean_name, 55)
      ),
      alpha_val = ifelse(is_focus, "yes", "no")
    )

  p_bar <- ggplot(plot_df,
                  aes(NES, reorder(pathway_label, NES),
                      fill = direction, alpha = alpha_val)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Rescued↑" = "#E64B35",
                                 "Rescued↓" = "#4DBBD5")) +
    scale_alpha_manual(values = c("yes" = 1.0, "no" = 0.5), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    labs(
      title    = paste0(site_label, "  |  挽救趋势 GSEA \u2013 ", gset_label),
      subtitle = paste0("NES > 0 = DB中上调且DB-Treated中被抑制的通路",
                        "   \u2605 = 重点通路   top ", top_n),
      x = "NES（基于 Rescue Score 排序）",
      y = NULL, fill = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom",
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8.5))

  ph <- max(5, top_n * 0.33 + 2.5)
  ggsave(paste0(fname, "_barplot.pdf"), p_bar, width = 11, height = ph)
  ggsave(paste0(fname, "_barplot.png"), p_bar, width = 11, height = ph, dpi = 180)

  return(invisible(res))
}

# ── 3.2 GO BP 富集 ──────────────────────────────────────────
.enrich_go <- function(sig_sym, bg_sym, title, out_prefix) {
  ez_sig <- tryCatch(
    bitr(sig_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  ez_bg <- tryCatch(
    bitr(bg_sym,  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (nrow(ez_sig) == 0) { message("  GO: 无法转换 EntrezID"); return(NULL) }

  go_res <- tryCatch(
    enrichGO(
      gene          = ez_sig$ENTREZID,
      universe      = if (nrow(ez_bg) > 0) ez_bg$ENTREZID else NULL,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    ),
    error = function(e) { message("  GO enrichment 出错: ", e$message); NULL }
  )
  if (is.null(go_res) || nrow(go_res@result) == 0) return(NULL)

  write.csv(go_res@result, paste0(out_prefix, "_GO_BP.csv"), row.names = FALSE)

  p_dot <- dotplot(go_res, showCategory = 20, title = title)
  ggsave(paste0(out_prefix, "_GO_BP_dotplot.pdf"), p_dot, width = 10, height = 8)
  ggsave(paste0(out_prefix, "_GO_BP_dotplot.png"), p_dot, width = 10, height = 8, dpi = 200)

  go_res
}

# ────────────────────────────────────────────────────────────
# 4. 主分析函数（按部位）
# ────────────────────────────────────────────────────────────

analyze_rescue <- function(site_name) {

  cat(sprintf("\n%s\n  \u25b6 挽救分析: %s\n%s\n",
              strrep("\u2500", 55), site_name, strrep("\u2500", 55)))

  out_s <- file.path(RSC_DIR, site_name)

  # ── 4.1 读取主脚本生成的数据 ────────────────────────────
  deg_db_wt  <- read.csv(file.path(OUT_DIR, site_name, "DEG_DB_vs_WT.csv"))
  deg_trt_db <- read.csv(file.path(OUT_DIR, site_name, "DEG_DBTreated_vs_DB.csv"))
  deg_trt_wt <- read.csv(file.path(OUT_DIR, site_name, "DEG_DBTreated_vs_WT.csv"))
  lrt_sig    <- read.csv(file.path(OUT_DIR, site_name, "LRT_results_sig.csv"))
  meta_s     <- read.csv(file.path(INT_DIR, paste0(site_name, "_metadata.csv")))

  # 加载 VSD 对象用于提取表达值
  vsd      <- readRDS(file.path(INT_DIR, paste0(site_name, "_vsd.rds")))
  vst_mat  <- assay(vsd)

  bg_sym   <- unique(na.omit(deg_db_wt$gene_name))
  bg_sym   <- bg_sym[bg_sym != ""]

  # ── 4.2 方法A：交集法 ───────────────────────────────────
  message("  \u2192 方法A：交集法")

  # 合并两个对比的 LFC 和显著性
  df_merge <- deg_db_wt %>%
    dplyr::select(gene_id, gene_name,
                  lfc_db_wt  = log2FoldChange,
                  padj_db_wt = padj) %>%
    inner_join(
      deg_trt_db %>%
        dplyr::select(gene_id,
                      lfc_trt_db  = log2FoldChange,
                      padj_trt_db = padj),
      by = "gene_id"
    ) %>%
    inner_join(
      deg_trt_wt %>%
        dplyr::select(gene_id,
                      lfc_trt_wt  = log2FoldChange,
                      padj_trt_wt = padj),
      by = "gene_id"
    ) %>%
    filter(!is.na(lfc_db_wt), !is.na(lfc_trt_db),
           !is.na(gene_name), gene_name != "")

  # 挽救基因定义（交集法）：
  #   1. 疾病效应显著（DB_vs_WT 达到阈值）
  #   2. 治疗效应显著且与疾病效应方向相反（DBTreated_vs_DB）
  rescue_df <- df_merge %>%
    mutate(
      # 挽救程度：治疗挽救了多少比例的疾病效应（1=完全挽救，>1=过矫正）
      rescue_pct = -lfc_trt_db / lfc_db_wt,

      # 分类
      rescue_type = case_when(
        # 完全挽救：疾病效应显著 + 治疗效应反向显著 + DBTreated_vs_WT 不显著
        padj_db_wt  < P$padj_cutoff &
          abs(lfc_db_wt) >= P$rescue_min_lfc &
          padj_trt_db < P$padj_cutoff &
          sign(lfc_trt_db) != sign(lfc_db_wt) &
          (is.na(padj_trt_wt) | padj_trt_wt >= P$padj_cutoff) ~ "Complete Rescue",

        # 部分挽救：疾病效应显著 + 治疗效应反向显著 + DBTreated_vs_WT 仍显著但变小
        padj_db_wt  < P$padj_cutoff &
          abs(lfc_db_wt) >= P$rescue_min_lfc &
          padj_trt_db < P$padj_cutoff &
          sign(lfc_trt_db) != sign(lfc_db_wt) &
          !is.na(padj_trt_wt) & padj_trt_wt < P$padj_cutoff &
          abs(lfc_trt_wt) < abs(lfc_db_wt) ~ "Partial Rescue",

        TRUE ~ "Other"
      )
    )

  rescue_genes <- rescue_df %>%
    filter(rescue_type %in% c("Complete Rescue", "Partial Rescue")) %>%
    arrange(rescue_type, desc(abs(lfc_db_wt)))

  write.csv(rescue_df,   file.path(out_s, "A_rescue_all_genes.csv"),   row.names = FALSE)
  write.csv(rescue_genes, file.path(out_s, "A_rescue_sig_genes.csv"),   row.names = FALSE)

  n_complete <- sum(rescue_genes$rescue_type == "Complete Rescue")
  n_partial  <- sum(rescue_genes$rescue_type == "Partial Rescue")
  message(sprintf("  完全挽救: %d 基因  | 部分挽救: %d 基因", n_complete, n_partial))

  # ── 4.3 方法B：Rescue Score 可视化 & GSEA ───────────────
  message("  \u2192 方法B：Rescue Score")

  # 仅对疾病效应足够强的基因计算 Rescue Score
  df_score <- df_merge %>%
    filter(abs(lfc_db_wt) >= P$rescue_min_lfc) %>%
    mutate(
      rescue_score = -lfc_trt_db / lfc_db_wt,
      # 过滤极端outlier（通常是 lfc_db_wt 趋于0时的噪声）
      rescue_score = ifelse(abs(rescue_score) > 5, NA, rescue_score)
    ) %>%
    filter(!is.na(rescue_score)) %>%
    arrange(desc(rescue_score))

  write.csv(df_score, file.path(out_s, "B_rescue_score_all.csv"), row.names = FALSE)

  # ── LFC 散点图：DB_vs_WT (x) vs DBTreated_vs_DB (y) ────
  # 理想的挽救基因在第2、4象限（方向相反）
  plot_df <- df_merge %>%
    mutate(
      category = case_when(
        gene_id %in% rescue_genes$gene_id[rescue_genes$rescue_type == "Complete Rescue"] ~ "Complete Rescue",
        gene_id %in% rescue_genes$gene_id[rescue_genes$rescue_type == "Partial Rescue"]  ~ "Partial Rescue",
        TRUE ~ "Other"
      ),
      label = ifelse(
        category != "Other" & abs(lfc_db_wt) > 1.5,
        gene_name, NA_character_
      )
    )

  p_scatter <- ggplot(plot_df,
                      aes(lfc_db_wt, lfc_trt_db, color = category)) +
    geom_point(data = filter(plot_df, category == "Other"),
               alpha = 0.25, size = 1.0) +
    geom_point(data = filter(plot_df, category != "Other"),
               alpha = 0.80, size = 1.8) +
    geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 20,
                    na.rm = TRUE, segment.color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
    # 标记挽救象限
    annotate("rect", xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = -0.5,
             alpha = 0.06, fill = "#E64B35") +
    annotate("rect", xmin = -Inf, xmax = -0.5, ymin = 0.5, ymax = Inf,
             alpha = 0.06, fill = "#4DBBD5") +
    annotate("text", x =  2.5, y = -3, label = "DB\u2191 \u2192 Treated\u2193",
             color = "#E64B35", size = 3.5, fontface = "bold") +
    annotate("text", x = -2.5, y =  3, label = "DB\u2193 \u2192 Treated\u2191",
             color = "#4DBBD5", size = 3.5, fontface = "bold") +
    scale_color_manual(values = c(
      "Complete Rescue" = "#E64B35",
      "Partial Rescue"  = "#F39B7F",
      "Other"           = "grey75"
    )) +
    labs(
      title = paste0(site_name, "  \u2014  挽救趋势象限图"),
      subtitle = paste0("x轴: DB_vs_WT (log2FC)   y轴: DBTreated_vs_DB (log2FC)\n",
                        "着色基因同时在两个对比中显著且方向相反"),
      x = "log2FC  (DB_Control vs WT_Control)",
      y = "log2FC  (DB_Treated vs DB_Control)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9))

  ggsave(file.path(out_s, "B_LFC_scatter_rescue.pdf"), p_scatter, width = 8, height = 7)
  ggsave(file.path(out_s, "B_LFC_scatter_rescue.png"), p_scatter, width = 8, height = 7, dpi = 200)

  # ── Rescue Score GSEA ────────────────────────────────────
  ranked_rescue <- df_score %>%
    arrange(desc(rescue_score)) %>%
    group_by(gene_name) %>%
    slice_max(abs(rescue_score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    { setNames(.$rescue_score, .$gene_name) }

  gsea_dir <- file.path(out_s, "GSEA")
  dir.create(gsea_dir, showWarnings = FALSE)

  .run_rescue_gsea(ranked_rescue, msig_hallmark,  "Hallmark", site_name, site_name, gsea_dir)
  .run_rescue_gsea(ranked_rescue, msig_go_bp,     "GOBP",     site_name, site_name, gsea_dir)
  .run_rescue_gsea(ranked_rescue, msig_kegg_gsea, "KEGG",     site_name, site_name, gsea_dir)

  # ── 4.4 方法C：LRT 显著基因模式聚类 ─────────────────────
  message("  \u2192 方法C：模式聚类（k-means）")

  # 计算三组均值（只用 WT_Control / DB_Control / DB_Treated）
  meta_3grp <- meta_s %>%
    filter(Group %in% c("WT_Control", "DB_Control", "DB_Treated"))
  rownames(meta_3grp) <- meta_3grp$rna_sample

  lrt_ids <- lrt_sig$gene_id[lrt_sig$gene_id %in% rownames(vst_mat)]
  if (length(lrt_ids) == 0) {
    message("  LRT 显著基因无法在 VST 矩阵中找到，跳过聚类")
  } else {
    vst_3grp  <- vst_mat[lrt_ids, meta_3grp$rna_sample, drop = FALSE]

    # 按组计算均值
    grp_means <- sapply(c("WT_Control", "DB_Control", "DB_Treated"), function(g) {
      cols <- meta_3grp$rna_sample[meta_3grp$Group == g]
      rowMeans(vst_3grp[, cols, drop = FALSE])
    })
    colnames(grp_means) <- c("WT_Control", "DB_Control", "DB_Treated")

    # Z-score 归一化（按行）
    grp_z <- t(scale(t(grp_means)))

    # 去除含NA的行（方差为0）
    grp_z <- grp_z[complete.cases(grp_z), ]

    # k-means 聚类
    k_use  <- min(P$kmeans_k, nrow(grp_z) - 1)
    km_res <- kmeans(grp_z, centers = k_use, nstart = 50, iter.max = 100)

    # 识别"挽救"模式的簇：
    #   DB_Control 均值与 WT_Control 同侧差距最大，
    #   DB_Treated 均值回到 WT_Control 方向
    cluster_summary <- as.data.frame(km_res$centers) %>%
      mutate(
        cluster    = seq_len(nrow(.)),
        n_genes    = as.integer(table(km_res$cluster)),
        # 疾病变化方向
        disease_dir = sign(DB_Control - WT_Control),
        # 治疗方向与疾病方向相反则为挽救趋势
        is_rescue   = sign(DB_Treated - DB_Control) != sign(DB_Control - WT_Control),
        # 挽救程度：DB_Treated 相对 WT_Control 的偏差比 DB_Control 小
        rescue_mag  = abs(DB_Treated - WT_Control) < abs(DB_Control - WT_Control)
      ) %>%
      arrange(desc(is_rescue), desc(rescue_mag))

    write.csv(cluster_summary,
              file.path(out_s, "C_kmeans_cluster_summary.csv"), row.names = FALSE)

    message("  聚类结果概况：")
    print(cluster_summary %>%
            dplyr::select(cluster, n_genes, WT_Control, DB_Control, DB_Treated,
                          is_rescue, rescue_mag))

    rescue_clusters <- cluster_summary$cluster[cluster_summary$is_rescue &
                                                 cluster_summary$rescue_mag]
    message(sprintf("  \u2605 挽救模式簇: %s",
                    paste(rescue_clusters, collapse = ", ")))

    # ── 聚类趋势折线图（每簇一行）──────────────────────────
    plot_cluster <- grp_z %>%
      as.data.frame() %>%
      mutate(gene_id  = rownames(.),
             cluster  = km_res$cluster) %>%
      pivot_longer(cols = c(WT_Control, DB_Control, DB_Treated),
                   names_to = "Group", values_to = "z_score") %>%
      mutate(
        Group = factor(Group,
                       levels = c("WT_Control", "DB_Control", "DB_Treated")),
        cluster_label = paste0(
          "Cluster ", cluster, " (n=",
          cluster_summary$n_genes[match(cluster, cluster_summary$cluster)], ")",
          ifelse(cluster %in% rescue_clusters, " \u2605", "")
        )
      )

    # 各簇平均趋势线
    cluster_mean <- plot_cluster %>%
      group_by(cluster, cluster_label, Group) %>%
      summarise(mean_z = mean(z_score), .groups = "drop")

    p_cluster <- ggplot() +
      # 背景个体基因线（透明）
      geom_line(data = plot_cluster,
                aes(Group, z_score, group = gene_id),
                color = "grey80", alpha = 0.15, linewidth = 0.3) +
      # 簇均值趋势线
      geom_line(data = cluster_mean,
                aes(Group, mean_z, group = cluster,
                    color = factor(cluster %in% rescue_clusters)),
                linewidth = 1.4) +
      geom_point(data = cluster_mean,
                 aes(Group, mean_z,
                     color = factor(cluster %in% rescue_clusters)),
                 size = 3.5) +
      scale_color_manual(
        values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
        labels = c("TRUE" = "挽救趋势", "FALSE" = "其他模式"),
        name   = NULL
      ) +
      facet_wrap(~ cluster_label, ncol = 3) +
      labs(
        title    = paste0(site_name, "  \u2014  LRT 显著基因表达模式聚类（k=", k_use, "）"),
        subtitle = "\u2605 = 具有 WT\u2192DB\u2192DB-Treated 挽救趋势的簇",
        x = NULL, y = "Z-score（VST 组均值）"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 20, hjust = 1, size = 9),
        strip.text   = element_text(size = 9, face = "bold"),
        legend.position = "bottom"
      )

    ggsave(file.path(out_s, "C_kmeans_pattern_plot.pdf"),
           p_cluster, width = 12, height = 8)
    ggsave(file.path(out_s, "C_kmeans_pattern_plot.png"),
           p_cluster, width = 12, height = 8, dpi = 200)

    # ── 挽救簇基因列表 + 热图 ─────────────────────────────
    if (length(rescue_clusters) > 0) {

      rescue_cluster_genes <- lrt_sig %>%
        filter(gene_id %in% names(km_res$cluster[km_res$cluster %in% rescue_clusters])) %>%
        mutate(cluster = km_res$cluster[match(gene_id,
                                              names(km_res$cluster))]) %>%
        dplyr::select(gene_id, gene_name, cluster, everything())

      write.csv(rescue_cluster_genes,
                file.path(out_s, "C_rescue_cluster_genes.csv"), row.names = FALSE)

      message(sprintf("  挽救簇基因总数: %d", nrow(rescue_cluster_genes)))

      # 热图（top P$heatmap_n 基因）
      top_ids_heat <- rescue_cluster_genes %>%
        arrange(padj) %>%
        head(P$heatmap_n) %>%
        pull(gene_id)

      s_order  <- meta_3grp %>% arrange(Group) %>% pull(rna_sample)
      mat_heat <- vst_mat[top_ids_heat, s_order, drop = FALSE]
      mat_heat <- t(scale(t(mat_heat)))

      gn_heat  <- lrt_sig$gene_name[match(rownames(mat_heat), lrt_sig$gene_id)]
      rownames(mat_heat) <- ifelse(is.na(gn_heat) | gn_heat == "",
                                   rownames(mat_heat), gn_heat)

      anno_col <- data.frame(
        Group   = meta_3grp$Group[match(s_order, meta_3grp$rna_sample)],
        row.names = s_order
      )
      grp_cols <- GROUP_COLORS[unique(as.character(anno_col$Group))]

      pdf(file.path(out_s, "C_rescue_cluster_heatmap.pdf"), width = 11, height = 13)
      pheatmap(
        mat_heat,
        annotation_col    = anno_col,
        annotation_colors = list(Group = grp_cols),
        color        = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        fontsize_row = 8,
        fontsize_col = 9,
        main = paste0(site_name,
                      "  \u2014  挽救模式基因热图（top ", P$heatmap_n, "）")
      )
      dev.off()

      # ── 挽救簇基因的 GO 富集 ───────────────────────────
      message("  \u2192 挽救簇 GO BP 富集")
      rg_syms <- unique(na.omit(rescue_cluster_genes$gene_name))
      rg_syms <- rg_syms[rg_syms != ""]
      if (length(rg_syms) >= 10) {
        .enrich_go(
          sig_sym    = rg_syms,
          bg_sym     = bg_sym,
          title      = paste0(site_name, "  |  挽救基因簇 GO BP"),
          out_prefix = file.path(out_s, "C_rescue_cluster")
        )
      }
    }
  }

  message(sprintf("  \u2713 %s 挽救分析完成\n", site_name))
}

# ────────────────────────────────────────────────────────────
# 5. 逐部位执行
# ────────────────────────────────────────────────────────────

for (site in c("Lung", "BO", "Oral")) {
  tryCatch(
    analyze_rescue(site),
    error = function(e) message(sprintf("!!! %s 挽救分析失败: %s", site, e$message))
  )
}

cat(sprintf("\n%s\n   挽救分析全部完成\n   输出: %s\n%s\n",
            strrep("=", 50), RSC_DIR, strrep("=", 50)))
