# 20260512趋化因子组间表达量比较 ----------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(limma)

# 1. 读取数据 ------------------------------------------------------------------
dat_raw <- read_excel(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/YZX/BvsA_deg_all.xlsx"
)

# 2. 筛选趋化因子相关基因 -------------------------------------------------------
chemokine_pattern <- "^(Ccl|Cxcl|Cx3cl|Xcl|Ccr|Cxcr|Cx3cr|Xcr)\\d"

dat_chemo <- dat_raw %>%
  filter(grepl(chemokine_pattern, gene_name, ignore.case = TRUE))

cat("筛选到趋化因子相关基因数量:", nrow(dat_chemo), "\n")
print(dat_chemo$gene_name)

# 3. 宽转长格式（用于t/Wilcoxon/置换检验） -------------------------------------
dat_long <- dat_chemo %>%
  pivot_longer(cols = -gene_name, names_to = "sample", values_to = "expr") %>%
  mutate(
    group = case_when(grepl("^A", sample) ~ "A", grepl("^B", sample) ~ "B"),
    expr = as.numeric(expr)
  )

# 4-A. t检验 + Wilcoxon --------------------------------------------------------
base_stats <- dat_long %>%
  group_by(gene_name) %>%
  summarise(
    mean_A = mean(expr[group == "A"], na.rm = TRUE),
    mean_B = mean(expr[group == "B"], na.rm = TRUE),
    sd_A = sd(expr[group == "A"], na.rm = TRUE),
    sd_B = sd(expr[group == "B"], na.rm = TRUE),
    t_pval = tryCatch(
      t.test(expr[group == "A"], expr[group == "B"])$p.value,
      error = function(e) NA_real_
    ),
    w_pval = tryCatch(
      wilcox.test(
        expr[group == "A"],
        expr[group == "B"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(log2FC = log2(mean_B / mean_A)) # B vs A

# 4-B. limma moderated t-test --------------------------------------------------
# limma 对 log2 表达量建模；加1避免 log(0)
expr_mat <- dat_chemo %>%
  tibble::column_to_rownames("gene_name") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
expr_mat <- log2(expr_mat + 1)

group_vec <- factor(c("B", "B", "B", "A", "A", "A", "A")) # 列顺序：B1-B3, A1-A4
design_mat <- model.matrix(~ 0 + group_vec)
colnames(design_mat) <- c("A", "B")
contrast_mat <- makeContrasts(BvsA = B - A, levels = design_mat)

fit <- lmFit(expr_mat, design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

limma_res <- topTable(fit2, coef = "BvsA", number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("gene_name") %>%
  select(
    gene_name,
    limma_logFC = logFC,
    limma_pval = P.Value,
    limma_padj = adj.P.Val
  )

# 4-C. 置换检验（Permutation test）--------------------------------------------
# 对每个基因，置换组标签 n_perm 次，计算 t 统计量的经验分布
set.seed(42)
n_perm <- 10000

perm_test <- function(vals_A, vals_B, n_perm) {
  obs_t <- tryCatch(
    t.test(vals_A, vals_B)$statistic,
    error = function(e) NA_real_
  )
  if (is.na(obs_t)) {
    return(NA_real_)
  }
  combined <- c(vals_A, vals_B)
  nA <- length(vals_A)
  perm_t <- replicate(n_perm, {
    idx <- sample(length(combined), nA)
    tryCatch(
      t.test(combined[idx], combined[-idx])$statistic,
      error = function(e) NA_real_
    )
  })
  mean(abs(perm_t) >= abs(obs_t), na.rm = TRUE)
}

perm_pvals <- dat_long %>%
  group_by(gene_name) %>%
  summarise(
    perm_pval = perm_test(expr[group == "A"], expr[group == "B"], n_perm),
    .groups = "drop"
  )

# 5. 汇总所有检验结果 ----------------------------------------------------------
stat_results <- base_stats %>%
  left_join(limma_res, by = "gene_name") %>%
  left_join(perm_pvals, by = "gene_name") %>%
  mutate(
    t_padj = p.adjust(t_pval, method = "BH"),
    w_padj = p.adjust(w_pval, method = "BH"),
    perm_padj = p.adjust(perm_pval, method = "BH")
  ) %>%
  select(
    gene_name,
    log2FC,
    mean_A,
    sd_A,
    mean_B,
    sd_B,
    t_pval,
    t_padj,
    w_pval,
    w_padj,
    limma_pval,
    limma_padj,
    perm_pval,
    perm_padj
  ) %>%
  arrange(limma_pval)

cat("\n===== 统计结果汇总（按 limma p值排序）=====\n")
print(stat_results, n = Inf)

# 6. 可视化：pure_scatter 风格（参考 compare_plot_optimized）------------------
# n=3 vs n=4，样本量 < 5，使用 pure_scatter 策略：均值 crossbar + SE errorbar + 散点
# 主题：theme_classic；配色：NPG；散点：shape 21 白色填充

NPG_colors <- c(A = "#E64B35", B = "#4DBBD5")

gene_order <- stat_results$gene_name

# p值标注数据框（limma 为主，附带其他方法）
pval_df <- stat_results %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_order),
    pval_label = sprintf(
      "t=%.3f  W=%.3f\nlimma=%.3f  perm=%.3f",
      t_pval,
      w_pval,
      limma_pval,
      perm_pval
    )
  ) %>%
  select(gene_name, pval_label, limma_pval)

# 计算每个 facet 的 y_max，用于放置 p 值文字
y_pos_df <- dat_long %>%
  mutate(gene_name = factor(gene_name, levels = gene_order)) %>%
  group_by(gene_name) %>%
  summarise(
    y_data_max = max(expr, na.rm = TRUE),
    y_range = diff(range(expr, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  # 文字放在数据最大值上方 5% 的位置
  mutate(y_text = y_data_max + y_range * 0.05)

pval_pos_df <- pval_df %>%
  left_join(y_pos_df, by = "gene_name")

dat_plot <- dat_long %>%
  mutate(gene_name = factor(gene_name, levels = gene_order))

p <- ggplot(dat_plot, aes(x = group, y = expr, fill = group, color = group)) +
  # 均值 crossbar：NPG 彩色填充 + 同色边线
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.3,
    linewidth = 1.2,
    fatten = 1,
    alpha = 0.55
  ) +
  # SE errorbar：与组色一致
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.18,
    linewidth = 0.9
  ) +
  # 散点：NPG 彩色填充，黑色细边框
  geom_jitter(
    shape = 21,
    size = 3,
    stroke = 0.6,
    color = "black",
    width = 0.15,
    alpha = 0.9
  ) +
  # 四种方法的 p 值文字标注（放在数据上方）
  geom_text(
    data = pval_pos_df,
    aes(x = 1.5, y = y_text, label = pval_label),
    inherit.aes = FALSE,
    size = 2.4,
    vjust = 0,
    hjust = 0.5,
    color = "grey30",
    lineheight = 1.2
  ) +
  scale_fill_manual(values = NPG_colors) +
  scale_color_manual(values = NPG_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.42))) +
  facet_wrap(~gene_name, scales = "free_y", ncol = 4) +
  labs(
    x = NULL,
    y = "Expression",
    title = "Chemokine-related genes: B vs A",
    subtitle = "Mean ± SE  |  p-values: t / Wilcoxon / limma / Permutation(10000×)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "italic"),
    strip.background = element_blank(),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.5),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/YZX/chemokine_BvsA.pdf",
  p,
  width = 4 * 4,
  height = ceiling(nrow(dat_chemo) / 4) * 5 + 1.5,
  device = "pdf"
)

cat("\n图已保存至 chemokine_BvsA.pdf\n")

# 7. 火山图：全转录组 BvsA 差异分析结果 ----------------------------------------
library(readxl)
library(ggrepel)

deg <- read.delim(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/YZX/all_compare.xls",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    log2FC = as.numeric(BvsA_log2FoldChange),
    pval = as.numeric(BvsA_pvalue),
    neglog10p = -log10(pval)
  ) %>%
  filter(!is.na(log2FC), !is.na(pval), pval > 0)

# 阈值设定：pvalue < 0.05，|log2FC| > 1
fc_cut <- 1
p_cut <- 0.05

deg <- deg %>%
  mutate(
    sig = case_when(
      pval < p_cut & log2FC > fc_cut ~ "Up",
      pval < p_cut & log2FC < -fc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up", "Down", "NS"))
  )

cat(sprintf(
  "\nUp: %d  |  Down: %d  |  NS: %d\n",
  sum(deg$sig == "Up"),
  sum(deg$sig == "Down"),
  sum(deg$sig == "NS")
))

# 需要标注的基因
label_genes <- c("Cxcl13", "Ccl17")
label_df <- deg %>% filter(gene_name %in% label_genes)

# 配色：NPG Up/Down + 灰色 NS
vol_colors <- c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey70")

# 点大小映射到 |log2FC|，基础大小加大，NS/显著点透明度分层
p_vol <- ggplot(
  deg,
  aes(x = log2FC, y = neglog10p, color = sig, size = abs(log2FC))
) +
  # NS 点：先画在底层，高透明度压缩视觉噪音
  geom_point(
    data = filter(deg, sig == "NS"),
    alpha = 0.2,
    shape = 16
  ) +
  # Up/Down 点：后画在上层，较低透明度突出显著点
  geom_point(
    data = filter(deg, sig != "NS"),
    alpha = 0.6,
    shape = 16
  ) +
  # 阈值参考线
  geom_vline(
    xintercept = c(-fc_cut, fc_cut),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(p_cut),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  # 标注基因：白底圆圈 + repel 斜体标签（size 固定，不随 FC 映射）
  geom_point(
    data = label_df,
    aes(size = abs(log2FC)),
    shape = 21,
    stroke = 0.9,
    fill = "white"
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = gene_name),
    size = 5,
    fontface = "italic",
    color = "black",
    box.padding = 0.6,
    point.padding = 0.6,
    segment.color = "grey40",
    segment.size = 0.6,
    max.overlaps = Inf
  ) +
  # 颜色 scale
  scale_color_manual(
    values = vol_colors,
    labels = c(
      Up = sprintf("Up (n=%d)", sum(deg$sig == "Up")),
      Down = sprintf("Down (n=%d)", sum(deg$sig == "Down")),
      NS = sprintf("NS (n=%d)", sum(deg$sig == "NS"))
    )
  ) +
  # 大小 scale：range 下限对应 FC≈0，上限对应最大 FC
  scale_size_continuous(
    name = expression("|" ~ log[2] ~ "FC|"),
    range = c(0.5, 6),
    breaks = c(1, 2, 4, 6),
    guide = guide_legend(
      override.aes = list(alpha = 0.7, color = "grey40"),
      order = 2
    )
  ) +
  labs(
    x = expression(log[2] ~ "Fold Change (B vs A)"),
    y = expression(-log[10] ~ italic(P)),
    title = "Volcano plot: B vs A",
    subtitle = sprintf("Thresholds: |log2FC| > %g  &  P < %g", fc_cut, p_cut),
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    legend.background = element_blank(),
    legend.spacing.y = unit(0.1, "cm"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.5)
  )

ggsave(
  "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/YZX/volcano_BvsA.pdf",
  p_vol,
  width = 7,
  height = 6,
  device = "pdf"
)

cat("火山图已保存至 volcano_BvsA.pdf\n")

# 8. 富集分析可视化 ------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)

# 加载 enrich_plot 框架
enrich_base <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/Scripts/r_functions/dev/enrich_plot/R"
source(file.path(enrich_base, "std_format.R"))
source(file.path(enrich_base, "registry.R"))
source(file.path(enrich_base, "styles/bubble.R"))
source(file.path(enrich_base, "styles/bar.R"))
source(file.path(enrich_base, "styles/lollipop.R"))
source(file.path(enrich_base, "styles/combined.R"))
source(file.path(enrich_base, "styles/point_bar.R"))

# 8-A. 提取差异基因 ------------------------------------------------------------
deg_sig <- deg %>% filter(pval < p_cut, abs(log2FC) > fc_cut)

genes_up <- deg_sig %>% filter(log2FC > fc_cut) %>% pull(gene_name)
genes_down <- deg_sig %>% filter(log2FC < -fc_cut) %>% pull(gene_name)
genes_all <- deg_sig %>% pull(gene_name)

cat(sprintf(
  "Up: %d  Down: %d  Total DEG: %d\n",
  length(genes_up),
  length(genes_down),
  length(genes_all)
))

# symbol → Entrez ID（clusterProfiler 需要 Entrez）
sym2entrez <- function(symbols) {
  bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db,
    drop = TRUE
  )$ENTREZID
}

entrez_up <- sym2entrez(genes_up)
entrez_down <- sym2entrez(genes_down)
entrez_all <- sym2entrez(genes_all)

# 8-B. GO 富集（BP）-----------------------------------------------------------
run_go <- function(entrez_ids, label = "") {
  if (length(entrez_ids) < 3) {
    message(label, ": 基因数 < 3，跳过 GO 富集")
    return(NULL)
  }
  enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE # 将 Entrez 映射回 symbol
  )
}

go_up <- run_go(entrez_up, "Up")
go_down <- run_go(entrez_down, "Down")
go_all <- run_go(entrez_all, "All DEG")

# 8-C. KEGG 富集 --------------------------------------------------------------
run_kegg <- function(entrez_ids, label = "") {
  if (length(entrez_ids) < 3) {
    message(label, ": 基因数 < 3，跳过 KEGG 富集")
    return(NULL)
  }
  enrichKEGG(
    gene = entrez_ids,
    organism = "mmu",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}

kegg_up <- run_kegg(entrez_up, "Up")
kegg_down <- run_kegg(entrez_down, "Down")
kegg_all <- run_kegg(entrez_all, "All DEG")

# 8-D. 可视化 -----------------------------------------------------------------
save_enrich_plot <- function(p, filename, w = 9, h = 7) {
  if (is.null(p)) {
    return(invisible(NULL))
  }
  ggsave(
    file.path(
      "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/YZX",
      filename
    ),
    p,
    width = w,
    height = h,
    device = "pdf"
  )
  cat("已保存:", filename, "\n")
}

# GO BP 气泡图：All DEG（Top 20）
if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {
  p_go_bubble <- plot_enrich(
    go_all,
    type = "combined",
    top_n = 20,
    x_var = "gene_ratio",
    color_low = "#d6eaf8",
    color_high = "#1a5276",
    title = "GO BP Enrichment – All DEGs (B vs A)"
  )
  save_enrich_plot(p_go_bubble, "enrich_GO_bubble_all.pdf", w = 16, h = 8)
}

# GO BP 双向条形图：Up vs Down
if (!is.null(go_up) && !is.null(go_down)) {
  up_df <- nrow(as.data.frame(go_up))
  down_df <- nrow(as.data.frame(go_down))
  if (up_df > 0 && down_df > 0) {
    p_go_dual <- plot_enrich(
      list(Up = go_up, Down = go_down),
      type = "bar",
      top_n = 10,
      group_colors = c("#E64B35", "#4DBBD5"),
      title = "GO BP – Up vs Down (B vs A)"
    )
    save_enrich_plot(p_go_dual, "enrich_GO_bar_dual.pdf", w = 10, h = 8)
  }
}

# KEGG 气泡图：All DEG
if (!is.null(kegg_all) && nrow(as.data.frame(kegg_all)) > 0) {
  # KEGG 结果 readable=FALSE，需手动将 Entrez 映射回 symbol 用于展示
  kegg_all_r <- setReadable(
    kegg_all,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID"
  )
  p_kegg_bubble <- plot_enrich(
    kegg_all_r,
    type = "bubble",
    top_n = 20,
    x_var = "gene_ratio",
    color_low = "#4DBBD5",
    color_high = "#E64B35",
    title = "KEGG Enrichment – All DEGs (B vs A)"
  )
  save_enrich_plot(p_kegg_bubble, "enrich_KEGG_bubble_all.pdf", w = 9, h = 8)
}

# KEGG 双向条形图：Up vs Down
if (!is.null(kegg_up) && !is.null(kegg_down)) {
  up_k <- nrow(as.data.frame(kegg_up))
  down_k <- nrow(as.data.frame(kegg_down))
  if (up_k > 0 && down_k > 0) {
    p_kegg_dual <- plot_enrich(
      list(
        Up = setReadable(kegg_up, org.Mm.eg.db, "ENTREZID"),
        Down = setReadable(kegg_down, org.Mm.eg.db, "ENTREZID")
      ),
      type = "bar",
      top_n = 10,
      group_colors = c("#E64B35", "#4DBBD5"),
      title = "KEGG – Up vs Down (B vs A)"
    )
    save_enrich_plot(p_kegg_dual, "enrich_KEGG_bar_dual.pdf", w = 16, h = 8)
  }
}

# 8-E. 棒棒糖图（lollipop）---------------------------------------------------
# GO BP All DEG：单组渐变棒棒糖，点大小映射 gene_count
if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {
  p_go_lollipop <- plot_enrich(
    go_all,
    type = "lollipop",
    top_n = 20,
    size_var = "gene_count",
    size_range = c(3, 10),
    group_colors = c("#d6eaf8", "#1a5276"),
    title = "GO BP – All DEGs (lollipop)"
  )
  save_enrich_plot(p_go_lollipop, "enrich_GO_lollipop_all.pdf", w = 9, h = 8)
}

# GO BP 双向棒棒糖：Up vs Down
if (
  !is.null(go_up) &&
    !is.null(go_down) &&
    nrow(as.data.frame(go_up)) > 0 &&
    nrow(as.data.frame(go_down)) > 0
) {
  p_go_lollipop_dual <- plot_enrich(
    list(Up = go_up, Down = go_down),
    type = "lollipop",
    top_n = 10,
    size_var = "gene_count",
    group_colors = c("#E64B35", "#4DBBD5"),
    title = "GO BP – Up vs Down (dual lollipop)"
  )
  save_enrich_plot(
    p_go_lollipop_dual,
    "enrich_GO_lollipop_dual.pdf",
    w = 10,
    h = 8
  )
}

# KEGG 双向棒棒糖
if (
  !is.null(kegg_up) &&
    !is.null(kegg_down) &&
    nrow(as.data.frame(kegg_up)) > 0 &&
    nrow(as.data.frame(kegg_down)) > 0
) {
  p_kegg_lollipop_dual <- plot_enrich(
    list(
      Up = setReadable(kegg_up, org.Mm.eg.db, "ENTREZID"),
      Down = setReadable(kegg_down, org.Mm.eg.db, "ENTREZID")
    ),
    type = "lollipop",
    top_n = 10,
    size_var = "gene_count",
    group_colors = c("#E64B35", "#4DBBD5"),
    title = "KEGG – Up vs Down (dual lollipop)"
  )
  save_enrich_plot(
    p_kegg_lollipop_dual,
    "enrich_KEGG_lollipop_dual.pdf",
    w = 10,
    h = 8
  )
}

# 8-F. point_bar：GO BP + KEGG 多分类双面板拼图 --------------------------------
# 需要 patchwork；GO 需同时运行 CC/MF，这里用 GO BP + KEGG 两个分类演示
library(patchwork)

pb_input <- list()
if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {
  pb_input[["GO BP"]] <- go_all
}
if (!is.null(kegg_all) && nrow(as.data.frame(kegg_all)) > 0) {
  pb_input[["KEGG"]] <- setReadable(kegg_all, org.Mm.eg.db, "ENTREZID")
}


if (length(pb_input) >= 1) {
  p_point_bar <- plot_enrich(
    pb_input,
    type = "point_bar",
    top_n = 8,
    category_colors = c("GO BP" = "#4DBBD5", "KEGG" = "#E64B35"),
    title = "Enrichment – GO BP & KEGG (point_bar)"
  )
  save_enrich_plot(p_point_bar, "enrich_point_bar.pdf", w = 12, h = 9)
}

# 8-G. combined：多分类综合可视化（左侧色块+气泡，右侧圆角条形）--------------
# 同时展示 GO BP 和 KEGG，每类 Top 5
combined_input <- list()
if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {
  combined_input[["GO BP"]] <- go_all
}
if (!is.null(kegg_all) && nrow(as.data.frame(kegg_all)) > 0) {
  combined_input[["KEGG"]] <- setReadable(kegg_all, org.Mm.eg.db, "ENTREZID")
}

if (length(combined_input) >= 1) {
  p_combined <- plot_enrich(
    combined_input,
    type = "combined",
    top_n = 5,
    category_colors = c("GO BP" = "#4DBBD5", "KEGG" = "#E64B35"),
    show_gene_ids = TRUE,
    max_genes_shown = 4,
    title = "Enrichment – GO BP & KEGG (combined)"
  )
  save_enrich_plot(p_combined, "enrich_combined.pdf", w = 13, h = 8)
}

cat("\n富集分析完成（框架样式）。\n")

# 9. 参考富集重画样式的富集可视化 -----------------------------------------------
library(patchwork)
library(stringr)

# ── 辅助：enrichResult → 整理好的 data.frame -----------------------------------
prep_df <- function(enrich_res, top_n = 15, group_label = NULL) {
  if (is.null(enrich_res)) {
    return(NULL)
  }
  df <- as.data.frame(enrich_res)
  if (nrow(df) == 0) {
    return(NULL)
  }
  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      }),
      neg_log10_padj = -log10(p.adjust),
      Description_wrap = str_wrap(Description, width = 40)
    )
  if (!is.null(group_label)) {
    df$Change <- group_label
  }
  df
}

# 9-A. RichFactor 渐变条形图（参考 Red模块富集 样式）--------------------------
# x = GeneRatio，fill = -log10(p.adjust)，通路名用白色斜体写在条形内
plot_richbar <- function(
  df,
  title = "",
  fill_low = "#4575b4",
  fill_high = "#d73027"
) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  df <- df %>%
    mutate(
      Description_wrap = factor(
        Description_wrap,
        levels = rev(unique(Description_wrap))
      )
    )

  ggplot(
    df,
    aes(x = GeneRatio_num, y = Description_wrap, fill = neg_log10_padj)
  ) +
    geom_col(alpha = 0.85, width = 0.75) +
    geom_text(
      aes(x = max(GeneRatio_num) * 0.015, label = Description_wrap),
      hjust = 0,
      fontface = "italic",
      size = 3.2,
      color = "white",
      lineheight = 0.88
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradient(
      low = fill_low,
      high = fill_high,
      name = expression(-log[10](p.adjust))
    ) +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

two_colors <- c(Up = "#E64B35", Down = "#4DBBD5")

# GO BP All DEG 渐变条形图
p9a_go <- plot_richbar(
  prep_df(go_all, 15),
  title = "GO BP – All DEGs (RichFactor bar)",
  fill_low = "#4575b4",
  fill_high = "#d73027"
)
# KEGG All DEG 渐变条形图
kegg_all_r2 <- if (!is.null(kegg_all)) {
  setReadable(kegg_all, org.Mm.eg.db, "ENTREZID")
} else {
  NULL
}
p9a_kegg <- plot_richbar(
  prep_df(kegg_all_r2, 15),
  title = "KEGG – All DEGs (RichFactor bar)",
  fill_low = "#4575b4",
  fill_high = "#d73027"
)

if (!is.null(p9a_go)) {
  save_enrich_plot(p9a_go, "enrich2_richbar_GO.pdf", w = 9, h = 7)
}
if (!is.null(p9a_kegg)) {
  save_enrich_plot(p9a_kegg, "enrich2_richbar_KEGG.pdf", w = 9, h = 7)
}

# 9-B. 左右拼图（参考 make_pair_plot 样式）-------------------------------------
# 合并 Up/Down 富集结果；左=点线图(-log10 p.adjust, x轴反转)，右=双向条形(GeneRatio)
make_pair_plot2 <- function(df_up, df_down, category_name = "") {
  if (is.null(df_up) && is.null(df_down)) {
    return(NULL)
  }

  # 合并，给每个 term 打 Up/Down 标签
  plot_data <- bind_rows(
    if (!is.null(df_up)) mutate(df_up, Change = "Up", score = GeneRatio_num),
    if (!is.null(df_down)) {
      mutate(df_down, Change = "Down", score = -GeneRatio_num)
    }
  ) %>%
    arrange(desc(abs(score))) %>%
    slice_head(n = 16) %>%
    mutate(
      Description_wrap = factor(
        Description_wrap,
        levels = rev(unique(Description_wrap))
      ),
      label_x = ifelse(
        score >= 0,
        -max(abs(score)) * 0.02,
        max(abs(score)) * 0.02
      ),
      hjust = ifelse(score >= 0, 1, 0)
    )

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  max_sig <- max(plot_data$neg_log10_padj, na.rm = TRUE)
  max_score <- max(abs(plot_data$score), na.rm = TRUE)

  # 左图：-log10(p.adjust) 点线图，x 轴反转
  p_left <- ggplot(
    plot_data,
    aes(x = neg_log10_padj, y = Description_wrap, group = 1)
  ) +
    geom_line(color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = Change), size = 3.5) +
    scale_color_manual(values = two_colors) +
    scale_x_reverse(limits = c(max_sig * 1.08, 0)) +
    labs(x = expression(-log[10](p.adjust)), y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = margin(5, 2, 5, 5)
    )

  # 右图：GeneRatio 双向条形图，通路名写在条形内
  p_right <- ggplot(
    plot_data,
    aes(x = score, y = Description_wrap, fill = Change)
  ) +
    geom_col(width = 0.72, alpha = 0.85) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_text(
      aes(x = label_x, label = Description_wrap, hjust = hjust),
      size = 2.8,
      lineheight = 0.88,
      color = "black"
    ) +
    scale_fill_manual(values = two_colors) +
    scale_x_continuous(limits = c(-max_score * 1.3, max_score * 1.3)) +
    labs(x = "Gene Ratio (Up +, Down −)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = margin(5, 5, 5, 2)
    )

  # 分类标题用 theme_void 注释
  title_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = category_name,
      fontface = "bold",
      size = 5
    ) +
    theme_void()

  title_plot /
    (p_left + p_right + plot_layout(widths = c(0.55, 1.45))) +
    plot_layout(heights = c(0.08, 1))
}

# GO BP 左右拼图
go_up_r <- if (!is.null(go_up)) go_up else NULL
go_down_r <- if (!is.null(go_down)) go_down else NULL
p9b_go <- make_pair_plot2(prep_df(go_up_r, 10), prep_df(go_down_r, 10), "GO BP")
if (!is.null(p9b_go)) {
  save_enrich_plot(p9b_go, "enrich2_pairplot_GO.pdf", w = 11, h = 8)
}

# KEGG 左右拼图
kegg_up_r2 <- if (!is.null(kegg_up)) {
  setReadable(kegg_up, org.Mm.eg.db, "ENTREZID")
} else {
  NULL
}
kegg_down_r2 <- if (!is.null(kegg_down)) {
  setReadable(kegg_down, org.Mm.eg.db, "ENTREZID")
} else {
  NULL
}
p9b_kegg <- make_pair_plot2(
  prep_df(kegg_up_r2, 10),
  prep_df(kegg_down_r2, 10),
  "KEGG"
)
if (!is.null(p9b_kegg)) {
  save_enrich_plot(p9b_kegg, "enrich2_pairplot_KEGG.pdf", w = 11, h = 8)
}

# 9-C. 分面条形图（参考 ONTOLOGY facet 样式）-----------------------------------
# 将 GO BP 和 KEGG 作为两个 ONTOLOGY，facet_wrap 展示，条形内写通路名
make_facet_bar <- function(enrich_list, top_n_each = 8, title = "") {
  df_all <- bind_rows(lapply(names(enrich_list), function(nm) {
    res <- enrich_list[[nm]]
    if (is.null(res)) {
      return(NULL)
    }
    df <- prep_df(res, top_n_each)
    if (is.null(df)) {
      return(NULL)
    }
    df$ONTOLOGY <- nm
    df
  }))
  if (is.null(df_all) || nrow(df_all) == 0) {
    return(NULL)
  }

  # 每个 ONTOLOGY 内部按 GeneRatio 从小到大排，固定 y 轴顺序
  df_all <- df_all %>%
    group_by(ONTOLOGY) %>%
    arrange(GeneRatio_num, .by_group = TRUE) %>%
    mutate(
      Description_wrap = factor(
        Description_wrap,
        levels = unique(Description_wrap)
      )
    ) %>%
    ungroup()

  ont_colors <- c(
    "GO BP" = "#4DBBD5",
    "KEGG" = "#E64B35"
  )

  ggplot(
    df_all,
    aes(x = GeneRatio_num, y = Description_wrap, fill = ONTOLOGY)
  ) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(
      aes(x = max(GeneRatio_num) * 0.01, label = Description_wrap),
      hjust = 0,
      size = 2.8,
      lineheight = 0.88,
      fontface = "italic"
    ) +
    facet_wrap(~ONTOLOGY, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = ont_colors) +
    expand_limits(x = max(df_all$GeneRatio_num) * 1.2) +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )
}

facet_input <- list()
if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {
  facet_input[["GO BP"]] <- go_all
}
if (!is.null(kegg_all_r2) && nrow(as.data.frame(kegg_all_r2)) > 0) {
  facet_input[["KEGG"]] <- kegg_all_r2
}

if (length(facet_input) > 0) {
  p9c <- make_facet_bar(
    facet_input,
    top_n_each = 10,
    title = "Enrichment – All DEGs (faceted bar)"
  )
  if (!is.null(p9c)) {
    save_enrich_plot(p9c, "enrich2_facetbar.pdf", w = 14, h = 9)
  }
}

# 9-D. 直接 ggplot2 气泡图（带外圈、网格背景）---------------------------------
plot_bubble_direct <- function(
  enrich_res,
  top_n = 15,
  title = "",
  color_low = "#4DBBD5",
  color_high = "#E64B35",
  fill_alpha = 0.8,
  ring_stroke = 1.8
) {
  if (is.null(enrich_res)) {
    return(NULL)
  }
  df <- as.data.frame(enrich_res)
  if (nrow(df) == 0) {
    return(NULL)
  }

  df <- df %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      }),
      neg_log10_padj = -log10(p.adjust)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Description = str_wrap(Description, width = 40),
      Description = factor(Description, levels = rev(unique(Description)))
    )

  size_range <- c(3, 12)

  ggplot(df, aes(x = GeneRatio_num, y = Description)) +
    # 单层 shape=21：fill 映射渐变（半透明），color 由 after_scale 派生为同色但强制不透明
    # stroke 控制外圈线宽
    geom_point(
      aes(
        size = Count,
        fill = neg_log10_padj,
        color = after_scale(scales::alpha(fill, 1)) # 外圈与 fill 同色，alpha 强制为 1
      ),
      shape = 21,
      alpha = fill_alpha, # 仅影响 fill；color 已通过 after_scale 锁定 alpha=1
      stroke = ring_stroke # 外圈宽度（单位：mm），调大可加粗外圈
    ) +
    scale_size_continuous(
      name = "Gene Count",
      range = size_range,
      breaks = pretty(df$Count, n = 4)
    ) +
    scale_fill_gradient(
      low = color_low,
      high = color_high,
      name = expression(-log[10](p.adjust))
    ) +
    guides(
      fill = guide_colorbar(barwidth = 0.8, barheight = 6, order = 1),
      size = guide_legend(
        override.aes = list(
          shape = 21,
          fill = "grey60",
          color = "grey40",
          stroke = 1.2
        ),
        order = 2
      )
    ) +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey60", linewidth = 0.6),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key = element_blank()
    )
}

# KEGG 气泡图（带外圈+网格）
if (exists("kegg_all_r") && !is.null(kegg_all_r)) {
  p9d_kegg <- plot_bubble_direct(
    kegg_all_r,
    top_n = 15,
    title = "KEGG Enrichment – All DEGs (B vs A)",
    color_low = "#4DBBD5",
    color_high = "#E64B35"
  )
  if (!is.null(p9d_kegg)) {
    save_enrich_plot(p9d_kegg, "enrich3_bubble_KEGG.pdf", w = 9, h = 7)
  }
}

# GO BP 气泡图（带外圈+网格）
if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {
  p9d_go <- plot_bubble_direct(
    go_all,
    top_n = 15,
    title = "GO BP Enrichment – All DEGs (B vs A)",
    color_low = "#74b9ff",
    color_high = "#d63031"
  )
  if (!is.null(p9d_go)) {
    save_enrich_plot(p9d_go, "enrich3_bubble_GO.pdf", w = 8, h = 7)
  }
}

cat("\n富集分析完成。\n")
