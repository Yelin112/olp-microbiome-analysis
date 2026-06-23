# ============================================================================
# 使用示例: plot_enrich_trinity
# 富集分析三位一体可视化函数
# ============================================================================

# 加载函数
source("function.R")

# === 准备工作 ===

# 如果缺少必需包,请先安装:
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
# install.packages(c("ggplot2", "dplyr", "stringr", "ggalluvial",
#                    "patchwork", "ggnewscale", "tidyr", "dittoSeq"))

library(clusterProfiler)
library(org.Hs.eg.db)

# === 示例1: 返回独立对象(推荐方式) ===

cat("\n=== 示例1: 返回三个独立的ggplot对象 ===\n")

# 1.1 准备更多的基因列表(25个炎症/免疫相关基因)
my_genes <- c(
  "NFKB1",
  "TNF",
  "IL1B",
  "IL6",
  "CXCL8",
  "PTGS2",
  "MMP9",
  "TLR4",
  "MAPK14",
  "ICAM1",
  "JUN",
  "STAT3",
  "MYD88",
  "CCL2",
  "IL10",
  "TGFB1",
  "VEGFA",
  "HIF1A",
  "NFKBIA",
  "RELA",
  "IKBKB",
  "CHUK",
  "IL1A",
  "CXCL10",
  "CCL5"
)

cat("输入基因数:", length(my_genes), "\n")
cat("基因列表:", paste(head(my_genes, 10), collapse = ", "), "...\n")

# 1.2 转换ID
gene_ids <- bitr(
  my_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat("成功转换:", nrow(gene_ids), "个基因\n")

# 1.3 KEGG富集分析
kegg_result <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

cat("富集到的通路数:", nrow(kegg_result), "\n")

# 1.4 准备差异表达数据(模拟更真实的数据)
set.seed(123)
diff_genes <- data.frame(
  gene = my_genes,
  log2FC = rnorm(length(my_genes), mean = 1, sd = 1.5), # 更真实的log2FC分布
  padj = 10^(-runif(length(my_genes), 2, 8)), # 更真实的p值分布
  stringsAsFactors = FALSE
)

cat("\n差异表达基因数:", nrow(diff_genes), "\n")
cat("log2FC范围:", round(range(diff_genes$log2FC), 2), "\n")

# 1.5 返回三个独立对象(默认combine=FALSE)
plots <- plot_enrich_trinity(
  enrich_result = kegg_result,
  top_n = 6,
  diff_data = diff_genes
)

cat("\n✓ 返回了三个独立对象:\n")
cat("  - plots$bar_plot: 基因条形图\n")
cat("  - plots$sankey_plot: 基因-通路桑基图\n")
cat("  - plots$bubble_plot: 通路气泡图\n")

# 1.6 查看单个图形
cat("\n查看桑基图:\n")
print(plots$sankey_plot)

# 1.7 自定义组合 - 桑基图+气泡图(最常用)
cat("\n自定义组合(桑基图+气泡图):\n")
combined_sankey_bubble <- plots$sankey_plot +
  plots$bubble_plot +
  plot_layout(widths = c(3, 1))

ggsave(
  "example1_sankey_bubble.png",
  combined_sankey_bubble,
  width = 10,
  height = 8
)

# 1.8 自定义组合 - 三个都要但调整比例
cat("\n自定义组合(三个图形,自定义比例):\n")
combined_custom <- plots$bar_plot +
  plots$sankey_plot +
  plots$bubble_plot +
  plot_layout(widths = c(0.8, 2, 0.8))

ggsave("example1_custom_layout.png", combined_custom, width = 14, height = 8)

cat("\n✓ 示例1完成!生成文件:\n")
cat("  - example1_sankey_bubble.png (桑基图+气泡图)\n")
cat("  - example1_custom_layout.png (三图自定义布局)\n")

# === 示例2: 自动组合并保存文件 ===

cat("\n=== 示例2: 使用output_file自动组合并保存 ===\n")

# 指定output_file会自动设置combine=TRUE
plot_enrich_trinity(
  enrich_result = kegg_result,
  top_n = 6,
  diff_data = diff_genes,
  output_file = "example2_auto_combined.png"
)

cat("✓ 示例2完成!图片保存为: example2_auto_combined.png\n")
cat("提示: 指定output_file会自动组合三个图形\n")

# === 示例3: 手动设置combine=TRUE ===

cat("\n=== 示例3: 手动组合(不保存文件) ===\n")

# 设置combine=TRUE返回组合对象
combined <- plot_enrich_trinity(
  enrich_result = kegg_result,
  top_n = 6,
  diff_data = diff_genes,
  combine = TRUE, # 手动设置组合
  layout_widths = c(0.8, 2, 0.8) # 调整比例
)

cat("✓ 返回组合后的patchwork对象\n")
print(combined)

# 可以进一步调整
library(patchwork)
combined_final <- combined +
  plot_annotation(
    title = "富集分析三位一体可视化",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave("example3_manual_combined.png", combined_final, width = 12, height = 8)

cat("✓ 示例3完成!图片保存为: example3_manual_combined.png\n")

# === 示例4: 最常用组合 - 桑基图+气泡图 ===

cat("\n=== 示例4: 桑基图+气泡图组合(推荐) ===\n")

# 获取独立对象
plots <- plot_enrich_trinity(
  enrich_result = kegg_result,
  top_n = 8, # 更多通路
  diff_data = diff_genes
)

# 只组合桑基图和气泡图(最常用)
library(patchwork)
sankey_bubble <- plots$sankey_plot +
  plots$bubble_plot +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "KEGG通路富集分析",
    subtitle = paste0("富集基因: ", length(my_genes), "个"),
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave(
  "example4_sankey_bubble_recommended.png",
  sankey_bubble,
  width = 12,
  height = 8
)

cat("✓ 示例4完成!图片保存为: example4_sankey_bubble_recommended.png\n")
cat("提示: 这是最常用的组合方式,简洁清晰\n")

# === 示例5: 指定特定通路 ===

cat("\n=== 示例5: 手动选择特定通路 ===\n")

# 指定感兴趣的通路
pathways_of_interest <- c(
  "NF-kappa B signaling pathway",
  "Toll-like receptor signaling pathway",
  "TNF signaling pathway",
  "Cytokine-cytokine receptor interaction"
)

plots <- plot_enrich_trinity(
  enrich_result = kegg_result,
  pathways = pathways_of_interest,
  diff_data = diff_genes
)

# 保存桑基图+气泡图
combined <- plots$sankey_plot +
  plots$bubble_plot +
  plot_layout(widths = c(3, 1))

ggsave("example5_selected_pathways.png", combined, width = 12, height = 8)

cat("✓ 示例5完成!图片保存为: example5_selected_pathways.png\n")

# === 示例6: GO富集分析 ===

cat("\n=== 示例6: GO富集分析可视化 ===\n")

# GO富集分析(生物学过程)
go_result <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

cat("富集到的GO terms:", nrow(go_result), "\n")

if (nrow(go_result) > 0) {
  plots_go <- plot_enrich_trinity(
    enrich_result = go_result,
    top_n = 6,
    diff_data = diff_genes
  )

  # 只保存桑基图+气泡图
  go_plot <- plots_go$sankey_plot +
    plots_go$bubble_plot +
    plot_layout(widths = c(3, 1))

  ggsave("example6_GO_BP.png", go_plot, width = 12, height = 8)
  cat("✓ 示例6完成!图片保存为: example6_GO_BP.png\n")
} else {
  cat("⚠️  未找到显著的GO terms\n")
}

# === 示例7: 自定义配色 ===

cat("\n=== 示例7: 自定义颜色方案 ===\n")

plots_custom <- plot_enrich_trinity(
  enrich_result = kegg_result,
  top_n = 6,
  diff_data = diff_genes,
  color_pvalue = c("blue", "orange"),
  color_fc = c("green", "purple"),
  bubble_size_range = c(3, 15)
)

custom_plot <- plots_custom$sankey_plot +
  plots_custom$bubble_plot +
  plot_layout(widths = c(3, 1))

ggsave("example7_custom_colors.png", custom_plot, width = 12, height = 8)

cat("✓ 示例7完成!图片保存为: example7_custom_colors.png\n")

# === 示例8: 批量处理GO三个本体 ===

cat("\n=== 示例8: 批量处理GO三个本体 ===\n")

go_onts <- c("BP", "CC", "MF")
go_names <- c("生物学过程", "细胞组分", "分子功能")

for (i in seq_along(go_onts)) {
  ont <- go_onts[i]
  ont_name <- go_names[i]

  cat("\n处理:", ont_name, "\n")

  go_res <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  if (nrow(go_res) > 5) {
    plots_ont <- plot_enrich_trinity(
      enrich_result = go_res,
      top_n = 6,
      diff_data = diff_genes
    )

    ont_plot <- plots_ont$sankey_plot +
      plots_ont$bubble_plot +
      plot_layout(widths = c(3, 1)) +
      plot_annotation(title = paste0("GO ", ont_name))

    ggsave(
      paste0("example8_GO_", ont, ".png"),
      ont_plot,
      width = 12,
      height = 8
    )

    cat("✓ 保存:", paste0("example8_GO_", ont, ".png\n"))
  } else {
    cat("⚠️  结果较少,跳过\n")
  }
}

cat("\n✓ 示例8完成!\n")

# === 完成 ===

cat("\n" %s% "=" %s% 50 %s% "\n")
cat("所有示例完成!\n")
cat("=" %s% 50 %s% "\n\n")

cat("生成的文件:\n")
cat("  - example1_sankey_bubble.png          (推荐: 桑基图+气泡图)\n")
cat("  - example1_custom_layout.png          (三图自定义布局)\n")
cat("  - example2_auto_combined.png          (自动组合)\n")
cat("  - example3_manual_combined.png        (手动组合)\n")
cat("  - example4_sankey_bubble_recommended.png  (最常用组合)\n")
cat("  - example5_selected_pathways.png      (指定通路)\n")
cat("  - example6_GO_BP.png                  (GO富集)\n")
cat("  - example7_custom_colors.png          (自定义配色)\n")
cat("  - example8_GO_*.png                   (GO三个本体)\n")

cat("\n核心功能说明:\n")
cat("1. 默认返回三个独立对象: $bar_plot, $sankey_plot, $bubble_plot\n")
cat("2. 推荐组合: sankey_plot + bubble_plot (最常用)\n")
cat("3. 自由组合: 使用patchwork自定义布局\n")
cat("4. 自动保存: 指定output_file会自动组合并保存\n")
cat("5. 使用更多基因(25个)获得更好的富集效果\n")

cat("\n使用建议:\n")
cat("✓ 优先使用独立对象方式,灵活性最高\n")
cat("✓ 桑基图+气泡图组合最清晰简洁\n")
cat("✓ 基因数建议20-50个,富集结果更有意义\n")
cat("✓ 可以只展示感兴趣的通路子集\n")
cat("✓ 使用patchwork自定义任意布局\n")
