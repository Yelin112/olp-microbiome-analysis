# ============================================================
# F/B比值分析：Firmicutes/Bacteroidota
# 三组比较：WT_Control vs DB_Control vs DB_Treated
# ============================================================

library(ggplot2)
library(dplyr)
library(ggpubr)

base_path <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
setwd(base_path)

if (!dir.exists("analysis_results/fb_ratio")) {
  dir.create("analysis_results/fb_ratio", recursive = TRUE)
}

# ============================================================
# 1. 读取数据
# ============================================================

# 读取门水平相对丰度
phylum_file <- "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.p.relative.xls"
phylum_data <- read.table(
  phylum_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# 移除最后的Tax_detail列
if ("Tax_detail" %in% colnames(phylum_data)) {
  phylum_data <- phylum_data[, !colnames(phylum_data) %in% "Tax_detail"]
}

# 读取metadata
metadata <- read.table(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1
)

# 统一样本名（下划线转点）
metadata$SampleID_converted <- gsub("_", ".", rownames(metadata))
metadata$TreatGroup <- paste(metadata$Condition, metadata$Treatment, sep = "_")

# 只保留三组
metadata_3group <- metadata %>%
  dplyr::filter(TreatGroup %in% c("WT_Control", "DB_Control", "DB_Treated"))

# 找共同样本
common_samples <- intersect(metadata_3group$SampleID_converted, colnames(phylum_data))
cat("共有样本数:", length(common_samples), "\n")

# ============================================================
# 2. 计算F/B比值
# ============================================================

# 提取Firmicutes和Bacteroidota（取交集样本）
phylum_sub <- phylum_data[, common_samples]

# 检查是否存在
cat("Firmicutes存在:", "Firmicutes" %in% rownames(phylum_sub), "\n")
cat("Bacteroidota存在:", "Bacteroidota" %in% rownames(phylum_sub), "\n")

firm <- as.numeric(phylum_sub["Firmicutes", ])
bact <- as.numeric(phylum_sub["Bacteroidota", ])

# 计算F/B比值（避免除以零）
fb_ratio <- ifelse(bact > 0, firm / bact, NA)

# 组合成数据框
fb_df <- data.frame(
  SampleID_converted = common_samples,
  Firmicutes = firm,
  Bacteroidota = bact,
  FB_ratio = fb_ratio,
  stringsAsFactors = FALSE
)

# 合并分组信息
fb_df <- fb_df %>%
  dplyr::left_join(
    metadata_3group %>%
      dplyr::select(SampleID_converted, Site, Condition, Treatment, TreatGroup),
    by = "SampleID_converted"
  )

# 设置分组顺序
fb_df$TreatGroup <- factor(
  fb_df$TreatGroup,
  levels = c("WT_Control", "DB_Control", "DB_Treated")
)

# 保存结果
write.csv(fb_df, "analysis_results/fb_ratio/fb_ratio_all_samples.csv", row.names = FALSE)
cat("F/B比值数据已保存\n")

# ============================================================
# 3. 统计检验
# ============================================================

sites <- unique(fb_df$Site)

stats_all <- data.frame()

for (site in sites) {
  site_data <- fb_df %>% dplyr::filter(Site == site, !is.na(FB_ratio))

  if (nrow(site_data) < 3) next

  groups <- levels(droplevels(site_data$TreatGroup))

  # Kruskal-Wallis整体检验
  kw <- kruskal.test(FB_ratio ~ TreatGroup, data = site_data)

  # Wilcoxon两两比较
  comparisons_pairs <- list(
    c("WT_Control", "DB_Control"),
    c("DB_Control", "DB_Treated"),
    c("WT_Control", "DB_Treated")
  )

  for (pair in comparisons_pairs) {
    if (all(pair %in% groups)) {
      pair_data <- site_data %>% dplyr::filter(TreatGroup %in% pair)
      wt <- wilcox.test(FB_ratio ~ TreatGroup, data = pair_data, exact = FALSE)

      # 计算中位数和效应量（rank-biserial correlation）
      n1 <- sum(pair_data$TreatGroup == pair[1])
      n2 <- sum(pair_data$TreatGroup == pair[2])
      r_effect <- 1 - (2 * wt$statistic) / (n1 * n2)

      stats_all <- rbind(stats_all, data.frame(
        Site = site,
        Comparison = paste(pair[1], "vs", pair[2]),
        KW_p = kw$p.value,
        Wilcox_p = wt$p.value,
        Effect_r = round(as.numeric(r_effect), 3),
        stringsAsFactors = FALSE
      ))
    }
  }
}

stats_all$Wilcox_p_adj <- p.adjust(stats_all$Wilcox_p, method = "BH")

write.csv(stats_all, "analysis_results/fb_ratio/fb_ratio_statistics.csv", row.names = FALSE)
cat("\n统计检验结果：\n")
print(stats_all)

# ============================================================
# 4. 可视化
# ============================================================

group_colors <- c(
  "WT_Control"  = "#8DD3C7",
  "DB_Control"  = "#BEBADA",
  "DB_Treated"  = "#FB8072"
)

comparisons_list <- list(
  c("WT_Control", "DB_Control"),
  c("DB_Control", "DB_Treated"),
  c("WT_Control", "DB_Treated")
)

theme_clean <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.8),
    axis.line.y = element_line(linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title  = element_text(size = 13, face = "bold"),
    plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(size = 11, face = "bold")
  )

# ===== 4a. 各Site单独出图 =====
plot_list <- list()

for (site in sites) {
  site_data <- fb_df %>%
    dplyr::filter(Site == site, !is.na(FB_ratio))

  if (nrow(site_data) < 3) next

  # 动态调整y轴上限（留空间给显著性标注）
  y_max <- max(site_data$FB_ratio, na.rm = TRUE) * 1.5

  p <- ggplot(site_data, aes(x = TreatGroup, y = FB_ratio, fill = TreatGroup)) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2,
      linewidth = 0.8,
      color = "black"
    ) +
    geom_point(
      shape = 21,
      size = 3,
      stroke = 0.8,
      color = "black",
      position = position_jitter(width = 0.12, seed = 42)
    ) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.signif",
      size = 4.5,
      tip.length = 0.02
    ) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = site,
      x = "",
      y = "Firmicutes/Bacteroidota Ratio"
    ) +
    theme_clean

  plot_list[[site]] <- p

  ggsave(
    paste0("analysis_results/fb_ratio/fb_ratio_", site, ".pdf"),
    p,
    width = 3.5,
    height = 4.5
  )
  cat("已保存:", site, "F/B比值图\n")
}

# ===== 4b. 四个Site合并出图 =====
if (length(plot_list) >= 2) {
  n_sites <- length(plot_list)
  ncols <- min(4, n_sites)

  p_combined <- ggpubr::ggarrange(
    plotlist = plot_list,
    ncol = ncols,
    nrow = ceiling(n_sites / ncols),
    common.legend = FALSE
  )

  ggsave(
    "analysis_results/fb_ratio/fb_ratio_all_sites_combined.pdf",
    p_combined,
    width = 3.5 * ncols,
    height = 4.5 * ceiling(n_sites / ncols)
  )
  cat("已保存合并图\n")
}

# ===== 4c. 按Site分面的总览图 =====
p_facet <- ggplot(
  fb_df %>% dplyr::filter(!is.na(FB_ratio)),
  aes(x = TreatGroup, y = FB_ratio, fill = TreatGroup)
) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    shape = 21,
    size = 2.5,
    stroke = 0.7,
    color = "black",
    position = position_jitter(width = 0.12, seed = 42)
  ) +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    label = "p.signif",
    size = 3.5,
    tip.length = 0.02
  ) +
  facet_wrap(~Site, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    x = "",
    y = "Firmicutes/Bacteroidota Ratio",
    title = "F/B Ratio across Sites"
  ) +
  theme_clean +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )

ggsave(
  "analysis_results/fb_ratio/fb_ratio_facet_all_sites.pdf",
  p_facet,
  width = 4 * length(sites),
  height = 5
)
cat("已保存分面总览图\n")

cat("\n所有F/B比值分析完成！\n")
cat("结果保存在: analysis_results/fb_ratio/\n")