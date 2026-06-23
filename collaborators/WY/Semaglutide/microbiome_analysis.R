# ============================================================
# 16S扩增子微生物组分析：Semaglutide对糖尿病小鼠微生态的影响
# 版本：2.1 (修正因子水平重复错误)
# ============================================================

# 安装和加载必需的包
packages <- c(
  "phyloseq",
  "ggplot2",
  "vegan",
  "dplyr",
  "tidyr",
  "DESeq2",
  "ape",
  "ComplexHeatmap",
  "RColorBrewer",
  "gridExtra",
  "reshape2",
  "scales",
  "microbiome",
  "ggalluvial"
)

# 检查并安装缺失的包
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("phyloseq", "DESeq2", "ComplexHeatmap", "microbiome")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 自定义函数：微生物丰度可视化（修正版）
# ============================================================

plot_microbial_abundance <- function(
  abundance_matrix,
  group_df = NULL,
  group_col = NULL,
  top_n = 10,
  plot_type = "relative",
  style = "bar",
  aggregate_by_group = FALSE,
  sort_by = "dominant",
  sort_decreasing = FALSE,
  palette = "Set3"
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggalluvial)

  # ================= 1. 数据预处理 =================

  if (plot_type == "relative") {
    mat_processed <- prop.table(as.matrix(abundance_matrix), 2)
    y_label <- "Relative Abundance"
    y_scale <- scale_y_continuous(labels = scales::percent, expand = c(0, 0))
  } else {
    mat_processed <- as.matrix(abundance_matrix)
    y_label <- "Absolute Abundance"
    y_scale <- scale_y_continuous(expand = c(0, 0))
  }

  target_taxon <- NULL
  if (!is.null(sort_by) && sort_by == "dominant") {
    target_taxon <- names(which.max(rowSums(mat_processed)))
  } else if (!is.null(sort_by) && sort_by %in% rownames(mat_processed)) {
    target_taxon <- sort_by
  }

  df <- as.data.frame(mat_processed)
  df$Taxonomy <- rownames(df)
  df_long <- tidyr::pivot_longer(
    df,
    cols = -Taxonomy,
    names_to = "SampleID",
    values_to = "Abundance"
  )

  if (!is.null(group_df) && !is.null(group_col)) {
    group_df$SampleID <- rownames(group_df)
    df_long <- dplyr::left_join(
      df_long,
      group_df[, c("SampleID", group_col)],
      by = "SampleID"
    )
    df_long <- df_long %>% dplyr::rename(Group = dplyr::all_of(group_col))
  }

  if (aggregate_by_group && !is.null(group_df) && !is.null(group_col)) {
    df_long <- df_long %>%
      dplyr::group_by(Group, Taxonomy) %>%
      dplyr::summarise(
        Abundance = mean(Abundance, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::rename(SampleID = Group)

    if (!is.null(target_taxon)) {
      order_df <- df_long %>%
        dplyr::filter(Taxonomy == target_taxon) %>%
        dplyr::arrange(
          if (sort_decreasing) dplyr::desc(Abundance) else Abundance
        )
      df_long$SampleID <- factor(df_long$SampleID, levels = order_df$SampleID)
    }
    need_facet <- FALSE
    x_axis_title <- "Groups"
  } else {
    if (!is.null(target_taxon)) {
      sort_values <- mat_processed[target_taxon, ]
      ordered_samples <- names(sort(sort_values, decreasing = sort_decreasing))
      df_long$SampleID <- factor(df_long$SampleID, levels = ordered_samples)
    }
    need_facet <- (!is.null(group_df) && !is.null(group_col))
    x_axis_title <- "Samples"
  }

  species_order <- df_long %>%
    dplyr::group_by(Taxonomy) %>%
    dplyr::summarise(Total = sum(Abundance)) %>%
    dplyr::arrange(dplyr::desc(Total)) %>%
    dplyr::pull(Taxonomy)
  target_species <- species_order[1:min(length(species_order), top_n)]

  df_final <- df_long %>%
    dplyr::mutate(
      Display_Taxonomy = ifelse(
        Taxonomy %in% target_species,
        Taxonomy,
        "Others"
      )
    ) %>%
    dplyr::group_by(SampleID, Display_Taxonomy)

  if (need_facet) {
    df_final <- df_final %>% dplyr::group_by(Group, .add = TRUE)
  }

  df_final <- df_final %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = 'drop')

  # ===== 修正：防止因子水平重复 =====
  # 如果target_species中已经包含"Others"，则不再添加
  if ("Others" %in% target_species) {
    # 移除target_species中的"Others"，然后在最后添加
    target_species_clean <- target_species[target_species != "Others"]
    unique_levels <- c(target_species_clean, "Others")
  } else {
    unique_levels <- c(target_species, "Others")
  }

  df_final$Display_Taxonomy <- factor(
    df_final$Display_Taxonomy,
    levels = unique_levels
  )

  # ================= 2. 主题定义 =================

  theme_for_alluvial <- theme_bw() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      panel.spacing.x = unit(0, units = "cm"),
      strip.background = element_rect(
        color = "black",
        fill = "white",
        linewidth = 0.8
      ),
      strip.placement = "outside",
      axis.line.y.left = element_line(color = "black", linewidth = 0.8),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),
      strip.text.x = element_text(size = 14, face = "bold"),
      axis.text = element_text(face = "bold", size = 12, color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      legend.title = element_text(face = "bold", size = 12, color = "black"),
      legend.text = element_text(face = "bold", size = 12, color = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 0.3)
    )

  theme_for_bar <- theme_bw() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

  # ================= 3. 绘图核心 =================

  if (style == "alluvial") {
    p <- ggplot(
      df_final,
      aes(
        x = SampleID,
        y = Abundance,
        fill = Display_Taxonomy,
        alluvium = Display_Taxonomy,
        stratum = Display_Taxonomy
      )
    ) +
      geom_flow(
        stat = "alluvium",
        lode.guidance = "frontback",
        width = 0.5,
        alpha = 0.5
      ) +
      geom_stratum(width = 0.5, alpha = 1, color = "black", linewidth = 0.1) +
      theme_for_alluvial
  } else {
    p <- ggplot(
      df_final,
      aes(x = SampleID, y = Abundance, fill = Display_Taxonomy)
    ) +
      geom_col(
        position = "stack",
        width = 0.8,
        color = "black",
        linewidth = 0.05
      ) +
      theme_for_bar
  }

  # ================= 4. 修饰与输出 =================
  p <- p +
    y_scale +
    labs(y = y_label, fill = "Taxonomy", x = x_axis_title)

  if (need_facet) {
    switch_param <- if (style == "alluvial") "x" else NULL
    p <- p +
      facet_grid(
        ~Group,
        scales = "free_x",
        space = "free_x",
        switch = switch_param
      )
  }

  if (length(unique(df_final$Display_Taxonomy)) <= 12) {
    p <- p + scale_fill_brewer(palette = palette)
  }

  return(p)
}

# ============================================================
# 1. 数据读取和预处理
# ============================================================

# 设置工作路径
base_path <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
setwd(base_path)

# 创建输出文件夹
if (!dir.exists("analysis_results")) {
  dir.create("analysis_results")
}
if (!dir.exists("analysis_results/figures")) {
  dir.create("analysis_results/figures")
}
if (!dir.exists("analysis_results/tables")) {
  dir.create("analysis_results/tables")
}

cat("开始读取数据...\n")

# 读取metadata
metadata <- read.table(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1
)

# 处理样本名命名不一致问题
convert_sample_name <- function(name) {
  gsub("_", ".", name)
}

metadata$SampleID_converted <- convert_sample_name(rownames(metadata))

# 读取ASV均一化count表（用于多样性分析）
asv_file <- "data/02.ASVanalysis/asv_even_taxon.txt"

if (file.exists(asv_file)) {
  lines <- readLines(asv_file)
  lines[1] <- sub("^#", "", lines[1])
  asv_count <- read.table(
    text = lines,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )

  if ("Taxonomy" %in% colnames(asv_count)) {
    taxonomy_info <- asv_count$Taxonomy
    asv_count <- asv_count[, !colnames(asv_count) %in% "Taxonomy"]
  }

  cat("ASV count表读取成功！维度:", dim(asv_count), "\n")
} else {
  stop("错误：未找到ASV count表文件！")
}

# 读取分类信息
taxonomy <- read.table(
  "data/02.ASVanalysis/Seq_taxonomy/seq_taxonomy_qza/taxonomy.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE,
  comment.char = ""
)

# 解析分类信息为各个层级
parse_taxonomy <- function(tax_string) {
  levels <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
  tax_split <- strsplit(tax_string, ";")[[1]]

  tax_clean <- sapply(tax_split, function(x) {
    gsub("^[kpcofgs]__", "", x)
  })

  if (length(tax_clean) < 7) {
    tax_clean <- c(tax_clean, rep("Unclassified", 7 - length(tax_clean)))
  }

  names(tax_clean) <- levels
  return(tax_clean)
}

tax_mat <- t(sapply(taxonomy$taxonomy, parse_taxonomy))
rownames(tax_mat) <- rownames(taxonomy)

# 确保样本名一致
metadata_samples <- metadata$SampleID_converted
asv_samples <- colnames(asv_count)

common_samples <- intersect(metadata_samples, asv_samples)
cat("共有样本数:", length(common_samples), "\n")
cat("Metadata中的样本数:", nrow(metadata), "\n")
cat("ASV count表中的样本数:", ncol(asv_count), "\n")

# 过滤并排序样本
metadata_filtered <- metadata[metadata$SampleID_converted %in% common_samples, ]
asv_count_filtered <- asv_count[, common_samples]

# 创建phyloseq对象（用于多样性分析）
OTU <- otu_table(as.matrix(asv_count_filtered), taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
SAMP <- sample_data(metadata_filtered)
rownames(SAMP) <- metadata_filtered$SampleID_converted

ps <- phyloseq(OTU, TAX, SAMP)

cat("\nPhyloseq对象创建成功！\n")
print(ps)

# ============================================================
# 2. 数据概览和质控
# ============================================================

cat("\n开始数据质控...\n")

# 计算每个样本的reads数
sample_sums_df <- data.frame(
  Sample = sample_names(ps),
  Reads = sample_sums(ps),
  Site = sample_data(ps)$Site,
  Condition = sample_data(ps)$Condition,
  Treatment = sample_data(ps)$Treatment
)

write.csv(
  sample_sums_df,
  "analysis_results/tables/sample_summary.csv",
  row.names = FALSE
)

# 绘制每个样本的reads数分布
p_reads <- ggplot(sample_sums_df, aes(x = Sample, y = Reads, fill = Site)) +
  geom_bar(stat = "identity") +
  facet_grid(Condition ~ Treatment, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "每个样本的Reads数", y = "Read Counts", x = "Sample")

ggsave(
  "analysis_results/figures/sample_reads_distribution.pdf",
  p_reads,
  width = 14,
  height = 8
)

# ============================================================
# 3. 丰度分析 - 各分类水平（按Site分开绘制）
# ============================================================

cat("\n开始丰度分析...\n")

# 读取各分类水平的相对丰度数据并绘图
levels <- c(
  "p" = "Phylum",
  "c" = "Class",
  "o" = "Order",
  "f" = "Family",
  "g" = "Genus",
  "s" = "Species"
)

# 获取所有Site
sites <- unique(metadata_filtered$Site)

for (level_code in names(levels)) {
  level_name <- levels[level_code]

  file_path <- paste0(
    "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.",
    level_code,
    ".relative.xls"
  )

  if (file.exists(file_path)) {
    # 读取数据（不设置行名）
    abundance_data <- read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE
    )

    # 第一列是分类名
    tax_col <- colnames(abundance_data)[1]

    # 如果有重复的分类名，按名称分组求和合并
    if (any(duplicated(abundance_data[, 1]))) {
      # 只对数值列求和，排除文本列（如Tax_detail）
      abundance_data <- abundance_data %>%
        dplyr::group_by(!!rlang::sym(tax_col)) %>%
        dplyr::summarise(
          dplyr::across(where(is.numeric), sum),
          .groups = 'drop'
        ) %>%
        as.data.frame()
    }

    # 设置行名并移除第一列
    rownames(abundance_data) <- abundance_data[, 1]
    abundance_data <- abundance_data[, -1]

    # 如果还存在Tax_detail等文本列，移除
    if ("Tax_detail" %in% colnames(abundance_data)) {
      abundance_data <- abundance_data[,
        !colnames(abundance_data) %in% "Tax_detail"
      ]
    }

    # 筛选共同样本
    abundance_filtered <- abundance_data[, common_samples]

    # 保存完整丰度表
    write.csv(
      abundance_filtered,
      paste0("analysis_results/tables/", level_code, "_abundance_all.csv")
    )

    # ===== 为每个Site单独绘图 =====
    for (site in sites) {
      # 筛选该Site的样本
      site_samples <- metadata_filtered %>%
        dplyr::filter(Site == site) %>%
        dplyr::pull(SampleID_converted)

      if (length(site_samples) > 0) {
        # 筛选该Site的丰度数据
        abundance_site <- abundance_filtered[, site_samples]

        # 准备该Site的分组信息
        group_info_site <- metadata_filtered[
          metadata_filtered$SampleID_converted %in% site_samples,
        ]
        rownames(group_info_site) <- group_info_site$SampleID_converted

        # 创建组合分组变量（只包含Condition和Treatment）
        group_info_site$CombinedGroup <- paste(
          group_info_site$Condition,
          group_info_site$Treatment,
          sep = "_"
        )

        # ===== 样式1：堆积柱状图（按样本） =====
        p_bar_sample <- plot_microbial_abundance(
          abundance_matrix = abundance_site,
          group_df = group_info_site,
          group_col = "CombinedGroup",
          top_n = 10,
          plot_type = "relative",
          style = "bar",
          aggregate_by_group = FALSE,
          sort_by = "dominant",
          palette = "Set3"
        ) +
          labs(title = paste(site, "-", level_name, "丰度分布（堆积柱状图）"))

        ggsave(
          paste0(
            "analysis_results/figures/",
            level_code,
            "_bar_",
            site,
            "_by_sample.pdf"
          ),
          p_bar_sample,
          width = 12,
          height = 8
        )

        # ===== 样式2：冲击图（按样本） =====
        p_alluvial_sample <- plot_microbial_abundance(
          abundance_matrix = abundance_site,
          group_df = group_info_site,
          group_col = "CombinedGroup",
          top_n = 10,
          plot_type = "relative",
          style = "alluvial",
          aggregate_by_group = FALSE,
          sort_by = "dominant",
          palette = "Set3"
        ) +
          labs(title = paste(site, "-", level_name, "丰度分布（冲击图）"))

        ggsave(
          paste0(
            "analysis_results/figures/",
            level_code,
            "_alluvial_",
            site,
            "_by_sample.pdf"
          ),
          p_alluvial_sample,
          width = 12,
          height = 8
        )

        # ===== 样式3：按分组聚合（堆积柱状图） =====
        p_bar_group <- plot_microbial_abundance(
          abundance_matrix = abundance_site,
          group_df = group_info_site,
          group_col = "CombinedGroup",
          top_n = 10,
          plot_type = "relative",
          style = "bar",
          aggregate_by_group = TRUE,
          sort_by = "dominant",
          palette = "Set3"
        ) +
          labs(title = paste(site, "-", level_name, "平均丰度（按分组）"))

        ggsave(
          paste0(
            "analysis_results/figures/",
            level_code,
            "_bar_",
            site,
            "_by_group.pdf"
          ),
          p_bar_group,
          width = 10,
          height = 7
        )

        cat(paste("  完成", site, "-", level_name, "水平分析\n"))
      }
    }

    cat(paste("完成", level_name, "水平所有Site的分析\n\n"))
  }
}

# ============================================================
# 4. Alpha多样性分析
# ============================================================

cat("\n开始Alpha多样性分析...\n")

# 计算多种alpha多样性指数
alpha_div <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon", "Simpson")
)

# 添加metadata信息
alpha_div$Sample <- rownames(alpha_div)
alpha_div <- alpha_div %>%
  dplyr::left_join(
    metadata_filtered %>%
      dplyr::select(SampleID_converted, Site, Condition, Treatment) %>%
      dplyr::rename(Sample = SampleID_converted),
    by = "Sample"
  )

write.csv(
  alpha_div,
  "analysis_results/tables/alpha_diversity.csv",
  row.names = FALSE
)

# 定义清新配色
fresh_colors <- c(
  "#8DD3C7",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#FCCDE5"
)

# 为每个Site绘制alpha多样性箱线图
sites <- unique(alpha_div$Site)

for (site in sites) {
  site_data <- alpha_div %>% dplyr::filter(Site == site)

  site_data$TreatGroup <- paste(
    site_data$Condition,
    site_data$Treatment,
    sep = "_"
  )

  # Shannon指数
  p_shannon <- ggplot(
    site_data,
    aes(x = TreatGroup, y = Shannon, fill = TreatGroup)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
    scale_fill_manual(values = fresh_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = paste("Shannon多样性 -", site),
      x = "Group",
      y = "Shannon Index"
    )

  # Observed ASVs
  p_observed <- ggplot(
    site_data,
    aes(x = TreatGroup, y = Observed, fill = TreatGroup)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
    scale_fill_manual(values = fresh_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = paste("观察到的ASV数 -", site),
      x = "Group",
      y = "Observed ASVs"
    )

  # 合并图
  p_combined <- grid.arrange(p_shannon, p_observed, ncol = 2)

  ggsave(
    paste0("analysis_results/figures/alpha_diversity_", site, ".pdf"),
    p_combined,
    width = 9,
    height = 5
  )
}

# 统计检验
alpha_stats <- data.frame()

for (site in sites) {
  site_data <- alpha_div %>% dplyr::filter(Site == site)
  site_data$TreatGroup <- paste(
    site_data$Condition,
    site_data$Treatment,
    sep = "_"
  )

  kw_shannon <- kruskal.test(Shannon ~ TreatGroup, data = site_data)
  kw_observed <- kruskal.test(Observed ~ TreatGroup, data = site_data)

  alpha_stats <- rbind(
    alpha_stats,
    data.frame(
      Site = site,
      Metric = "Shannon",
      p_value = kw_shannon$p.value
    )
  )

  alpha_stats <- rbind(
    alpha_stats,
    data.frame(
      Site = site,
      Metric = "Observed",
      p_value = kw_observed$p.value
    )
  )
}

write.csv(
  alpha_stats,
  "analysis_results/tables/alpha_diversity_statistics.csv",
  row.names = FALSE
)

cat("\nAlpha多样性分析完成！\n")

# ============================================================
# 5. Beta多样性分析
# ============================================================

cat("\n开始Beta多样性分析...\n")

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

for (site in sites) {
  ps_site <- subset_samples(ps_rel, Site == site)

  if (nsamples(ps_site) > 3) {
    ord_site <- ordinate(ps_site, method = "PCoA", distance = "bray")

    sample_data(ps_site)$TreatGroup <- paste(
      sample_data(ps_site)$Condition,
      sample_data(ps_site)$Treatment,
      sep = "_"
    )

    p_pcoa <- plot_ordination(
      ps_site,
      ord_site,
      color = "TreatGroup",
      shape = "Condition"
    ) +
      geom_point(size = 5) +
      theme_bw() +
      stat_ellipse(aes(group = TreatGroup), type = "t", linetype = 2) +
      scale_color_manual(values = fresh_colors) +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 11),
        panel.grid = element_blank()
      ) + # 去除网格线
      labs(title = paste("PCoA (Bray-Curtis) -", site))

    ggsave(
      paste0("analysis_results/figures/beta_diversity_PCoA_", site, ".pdf"),
      p_pcoa,
      width = 5,
      height = 3.5
    )

    metadata_site <- data.frame(sample_data(ps_site))
    dist_site <- phyloseq::distance(ps_site, method = "bray")

    perm_result <- adonis2(
      dist_site ~ TreatGroup,
      data = metadata_site,
      permutations = 999
    )

    write.csv(
      as.data.frame(perm_result),
      paste0("analysis_results/tables/PERMANOVA_", site, ".csv")
    )
  }
}

cat("\nBeta多样性分析完成！\n")

# ============================================================
# 6. 差异丰度分析 (DESeq2)
# ============================================================

cat("\n开始差异丰度分析...\n")

class_abs <- read.table(
  "data/02.ASVanalysis/Taxa_abundance/Evenabs/asv_table.c.absolute.xls",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

class_abs_filtered <- class_abs[, common_samples]

for (site in sites) {
  site_samples <- metadata_filtered %>%
    dplyr::filter(Site == site) %>%
    dplyr::pull(SampleID_converted)

  if (length(site_samples) > 3) {
    class_site <- class_abs_filtered[, site_samples]
    metadata_site <- metadata_filtered[
      metadata_filtered$SampleID_converted %in% site_samples,
    ]

    metadata_site$TreatGroup <- paste(
      metadata_site$Condition,
      metadata_site$Treatment,
      sep = "_"
    )

    class_site_filtered <- class_site[rowSums(class_site) >= 10, ]

    dds <- DESeqDataSetFromMatrix(
      countData = round(class_site_filtered),
      colData = metadata_site,
      design = ~TreatGroup
    )

    dds <- DESeq(dds)

    if (
      "DB_Control" %in%
        metadata_site$TreatGroup &&
        "WT_Control" %in% metadata_site$TreatGroup
    ) {
      res_DB_WT <- results(
        dds,
        contrast = c("TreatGroup", "DB_Control", "WT_Control")
      )
      res_DB_WT_df <- as.data.frame(res_DB_WT)
      res_DB_WT_df$Taxon <- rownames(res_DB_WT_df)
      res_DB_WT_df <- res_DB_WT_df[order(res_DB_WT_df$padj), ]

      write.csv(
        res_DB_WT_df,
        paste0("analysis_results/tables/DESeq2_", site, "_DB_vs_WT.csv"),
        row.names = FALSE
      )
    }

    if (
      "DB_Semaglutide" %in%
        metadata_site$TreatGroup &&
        "DB_Control" %in% metadata_site$TreatGroup
    ) {
      res_Treat_Control <- results(
        dds,
        contrast = c("TreatGroup", "DB_Semaglutide", "DB_Control")
      )
      res_Treat_Control_df <- as.data.frame(res_Treat_Control)
      res_Treat_Control_df$Taxon <- rownames(res_Treat_Control_df)
      res_Treat_Control_df <- res_Treat_Control_df[
        order(res_Treat_Control_df$padj),
      ]

      write.csv(
        res_Treat_Control_df,
        paste0(
          "analysis_results/tables/DESeq2_",
          site,
          "_Treatment_vs_Control.csv"
        ),
        row.names = FALSE
      )

      res_Treat_Control_df$significant <- ifelse(
        res_Treat_Control_df$padj < 0.05 &
          abs(res_Treat_Control_df$log2FoldChange) > 1,
        "Significant",
        "Not Significant"
      )

      p_volcano <- ggplot(
        res_Treat_Control_df,
        aes(x = log2FoldChange, y = -log10(padj), color = significant)
      ) +
        geom_point(alpha = 0.6) +
        scale_color_manual(
          values = c("Significant" = "red", "Not Significant" = "gray")
        ) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        theme_bw() +
        labs(
          title = paste("差异丰度分析 -", site, "(Treatment vs Control)"),
          x = "Log2 Fold Change",
          y = "-Log10(Adjusted P-value)"
        )

      ggsave(
        paste0(
          "analysis_results/figures/volcano_",
          site,
          "_Treatment_vs_Control.pdf"
        ),
        p_volcano,
        width = 10,
        height = 7
      )
    }
  }
}

cat("\n差异丰度分析完成！\n")

# ============================================================
# 7. 生成分析报告摘要
# ============================================================

summary_text <- paste0(
  "16S扩增子微生物组分析完成！\n\n",
  "分析样本数: ",
  nrow(metadata_filtered),
  "\n",
  "ASV数: ",
  ntaxa(ps),
  "\n",
  "分析部位: ",
  paste(unique(metadata_filtered$Site), collapse = ", "),
  "\n\n",
  "输出文件位置:\n",
  "1. 图表: analysis_results/figures/\n",
  "2. 数据表: analysis_results/tables/\n\n",
  "主要分析内容:\n",
  "- 样本质控和reads分布\n",
  "- 各分类水平丰度分布（3种可视化风格）\n",
  "- Alpha多样性分析（Shannon、Observed ASVs）\n",
  "- Beta多样性分析（PCoA、PERMANOVA）\n",
  "- 差异丰度分析（DESeq2）\n\n",
  "版本：2.1 - 修正因子水平重复错误\n"
)

writeLines(summary_text, "analysis_results/ANALYSIS_SUMMARY.txt")
cat(summary_text)

cat("\n所有分析完成！请查看 analysis_results 文件夹中的结果。\n")


# ============================================================
# 8. 菌群失调评分分析 (Dysbiosis Score Analysis)
# ============================================================

cat("\n开始菌群失调评分分析...\n")

# 安装和加载dysbiosisR包
if (!require("dysbiosisR", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("microsud/dysbiosisR")
  library(dysbiosisR)
}

library(pROC)
library(ggpubr)

# 创建存储结果的文件夹
if (!dir.exists("analysis_results/dysbiosis")) {
  dir.create("analysis_results/dysbiosis")
}

# 定义清新配色
fresh_colors <- c(
  "#8DD3C7",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#FCCDE5"
)

# 定义比较组合
comparisons_list <- list(
  c("WT_Control", "DB_Control"),
  c("DB_Control", "DB_Treated"),
  c("WT_Control", "DB_Treated")
)

# 为每个Site计算菌群失调评分
for (site in sites) {
  cat(paste("\n========== 处理", site, "==========\n"))

  # 筛选该Site的样本
  ps_site <- subset_samples(ps_rel, Site == site)

  # 检查样本数
  if (nsamples(ps_site) < 5) {
    cat(paste("  ", site, "样本数不足，跳过\n"))
    next
  }

  # 获取样本元数据
  meta_site <- data.frame(sample_data(ps_site))
  meta_site$TreatGroup <- paste(
    meta_site$Condition,
    meta_site$Treatment,
    sep = "_"
  )

  # 定义参考样本（WT_Control = 健康对照）
  ref_samples <- rownames(meta_site[meta_site$TreatGroup == "WT_Control", ])

  if (length(ref_samples) < 2) {
    cat(paste("  ", site, "参考样本不足，跳过\n"))
    next
  }

  # 计算Bray-Curtis距离矩阵
  dist_mat <- phyloseq::distance(ps_site, method = "bray")

  # ===== 方法1: dysbiosisMedianCLV =====
  cat("  计算Median CLV评分...\n")
  dysbiosis_clv <- dysbiosisMedianCLV(
    ps_site,
    dist_mat = dist_mat,
    reference_samples = ref_samples
  )

  # ===== 方法2: combinedShannonJSD =====
  cat("  计算Shannon-JSD评分...\n")
  dysbiosis_shannon <- combinedShannonJSD(
    ps_site,
    reference_samples = ref_samples
  )

  # ===== 方法3: cloudStatistic =====
  cat("  计算CLOUD评分...\n")
  dysbiosis_cloud <- cloudStatistic(
    ps_site,
    dist_mat = dist_mat,
    reference_samples = ref_samples,
    ndim = -1,
    k_num = 5
  )

  # ===== 合并所有评分 =====
  dysbiosis_all <- dysbiosis_clv %>%
    dplyr::select(SampleID_converted, score, Site, Condition, Treatment) %>%
    dplyr::rename(MedianCLV = score) %>%
    dplyr::left_join(
      dysbiosis_shannon %>% dplyr::select(SampleID_converted, ShannonJSDScore),
      by = "SampleID_converted"
    ) %>%
    dplyr::left_join(
      dysbiosis_cloud %>% dplyr::select(SampleID_converted, log2Stats),
      by = "SampleID_converted"
    ) %>%
    dplyr::rename(CloudScore = log2Stats)

  # 添加TreatGroup
  dysbiosis_all <- dysbiosis_all %>%
    dplyr::left_join(
      meta_site %>%
        dplyr::select(SampleID_converted, TreatGroup),
      by = "SampleID_converted"
    )

  # ===== 保存评分表（重要！）=====
  cat("  保存评分数据...\n")
  write.csv(
    dysbiosis_all,
    paste0("analysis_results/dysbiosis/dysbiosis_scores_", site, ".csv"),
    row.names = FALSE
  )

  # 设置分组顺序
  dysbiosis_all$TreatGroup <- factor(
    dysbiosis_all$TreatGroup,
    levels = c("WT_Control", "DB_Control", "DB_Treated", "WT_Treated")
  )

  # ===== 可视化 =====
  cat("  生成可视化图表...\n")

  # MedianCLV
  p_clv <- ggplot(
    dysbiosis_all,
    aes(x = TreatGroup, y = MedianCLV, fill = TreatGroup)
  ) +
    # stat_summary(
    #   fun = mean,
    #   geom = "bar",
    #   alpha = 0.8,
    #   linewidth = 1.2,
    #   width = 0.4
    # ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2,
      linewidth = 1
    ) +
    geom_point(
      shape = 21,
      size = 2,
      stroke = 1,
      # fill = "white",
      color = "black",
      position = position_jitter(width = 0.1, seed = 123)
    ) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif",
      size = 4
    ) +
    scale_fill_manual(values = fresh_colors) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.8),
      axis.line.y = element_line(linewidth = 0.8)
    ) +
    labs(
      title = paste(site, "- Median CLV Dysbiosis Score"),
      x = "",
      y = "Median CLV Score"
    )

  ggsave(
    paste0("analysis_results/dysbiosis/", site, "_MedianCLV_barplot.pdf"),
    p_clv,
    width = 3.5,
    height = 4
  )

  # ShannonJSD
  p_shannon <- ggplot(
    dysbiosis_all,
    aes(x = TreatGroup, y = ShannonJSDScore, fill = TreatGroup)
  ) +
    # stat_summary(
    #   fun = mean,
    #   geom = "bar",
    #   alpha = 0.8,
    #   linewidth = 1.2,
    #   width = 0.4
    # ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2,
      linewidth = 1
    ) +
    geom_point(
      shape = 21,
      size = 2,
      stroke = 1,
      # fill = "white",
      color = "black",
      position = position_jitter(width = 0.1, seed = 123)
    ) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif",
      size = 4
    ) +
    scale_fill_manual(values = fresh_colors) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.8),
      axis.line.y = element_line(linewidth = 0.8)
    ) +
    labs(
      title = paste(site, "- Shannon-JSD Dysbiosis Score"),
      x = "",
      y = "Shannon-JSD Score"
    )

  ggsave(
    paste0("analysis_results/dysbiosis/", site, "_ShannonJSD_barplot.pdf"),
    p_shannon,
    width = 3.5,
    height = 4
  )

  # CloudScore
  p_cloud <- ggplot(
    dysbiosis_all,
    aes(x = TreatGroup, y = CloudScore, fill = TreatGroup)
  ) +
    # stat_summary(
    #   fun = mean,
    #   geom = "bar",
    #   alpha = 0.8,
    #   linewidth = 1.2,
    #   width = 0.4
    # ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2,
      linewidth = 1
    ) +
    geom_point(
      shape = 21,
      size = 2,
      stroke = 1,
      # fill = "white",
      color = "black",
      position = position_jitter(width = 0.1, seed = 123)
    ) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif",
      size = 4
    ) +
    scale_fill_manual(values = fresh_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) + # Cloud可能有负值
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.8),
      axis.line.y = element_line(linewidth = 0.8)
    ) +
    labs(
      title = paste(site, "- CLOUD Dysbiosis Score"),
      x = "",
      y = "CLOUD Score (log2)"
    )

  ggsave(
    paste0("analysis_results/dysbiosis/", site, "_CloudScore_barplot.pdf"),
    p_cloud,
    width = 3.5,
    height = 4
  )

  # ===== 统计检验 =====
  cat("  进行统计检验...\n")

  stats_results <- data.frame()

  if (
    "DB_Control" %in%
      dysbiosis_all$TreatGroup &&
      "DB_Treated" %in% dysbiosis_all$TreatGroup
  ) {
    # MedianCLV
    test_clv <- wilcox.test(
      MedianCLV ~ TreatGroup,
      data = dysbiosis_all %>%
        dplyr::filter(TreatGroup %in% c("DB_Control", "DB_Treated"))
    )

    # Shannon-JSD
    test_shannon <- wilcox.test(
      ShannonJSDScore ~ TreatGroup,
      data = dysbiosis_all %>%
        dplyr::filter(TreatGroup %in% c("DB_Control", "DB_Treated"))
    )

    # CloudScore
    test_cloud <- wilcox.test(
      CloudScore ~ TreatGroup,
      data = dysbiosis_all %>%
        dplyr::filter(TreatGroup %in% c("DB_Control", "DB_Treated"))
    )

    stats_results <- data.frame(
      Site = site,
      Method = c("MedianCLV", "ShannonJSD", "CloudScore"),
      P_value = c(test_clv$p.value, test_shannon$p.value, test_cloud$p.value),
      Comparison = "DB_Treated_vs_DB_Control"
    )

    write.csv(
      stats_results,
      paste0(
        "analysis_results/dysbiosis/",
        site,
        "_treatment_effect_stats.csv"
      ),
      row.names = FALSE
    )
  }

  cat(paste("  ", site, "分析完成！\n"))
}

cat("\n菌群失调评分分析完成！\n")
cat("结果保存在: analysis_results/dysbiosis/\n")

# 列出生成的文件
cat("\n生成的文件列表：\n")
dysbiosis_files <- list.files(
  "analysis_results/dysbiosis/",
  pattern = "*.csv|*.pdf"
)
for (f in dysbiosis_files) {
  cat("  -", f, "\n")
}


# ============================================================
# 差异分析：属水平 + 三组比较（WT_Control, DB_Control, DB_Treated）
# ============================================================

library(dplyr)


# 创建输出文件夹
if (!dir.exists("analysis_results/differential")) {
  dir.create("analysis_results/differential")
}

# ===== Step 1：读取属水平丰度数据 =====
genus_abs <- read.table(
  "data/02.ASVanalysis/Taxa_abundance/Evenabs/asv_table.g.absolute.xls",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

genus_abs_filtered <- genus_abs[, common_samples]

# 创建TreatGroup
metadata_filtered$TreatGroup <- paste(
  metadata_filtered$Condition,
  metadata_filtered$Treatment,
  sep = "_"
)

# 按Site分析
for (site in sites) {
  cat(paste("\n========== 分析", site, "（属水平）==========\n"))

  # 筛选该Site的样本
  site_samples <- metadata_filtered %>%
    dplyr::filter(Site == site) %>%
    dplyr::pull(SampleID_converted)

  if (length(site_samples) < 6) {
    cat("  样本数不足，跳过\n")
    next
  }

  # 准备数据
  abund_site <- genus_abs_filtered[, site_samples]

  # 获取分组信息
  meta_site <- metadata_filtered %>%
    dplyr::filter(SampleID_converted %in% site_samples)

  # ===== 只保留三组：WT_Control, DB_Control, DB_Treated =====
  meta_site <- meta_site %>%
    dplyr::filter(TreatGroup %in% c("WT_Control", "DB_Control", "DB_Treated"))

  if (nrow(meta_site) < 6) {
    cat("  三组样本数不足，跳过\n")
    next
  }

  # 更新样本列表和丰度矩阵
  site_samples_filtered <- meta_site$SampleID_converted
  abund_site <- abund_site[, site_samples_filtered]

  # 创建分组因子
  group_site <- factor(
    meta_site$TreatGroup,
    levels = c("WT_Control", "DB_Control", "DB_Treated") # 明确顺序
  )
  names(group_site) <- meta_site$SampleID_converted

  cat("  样本数:", length(site_samples_filtered), "\n")
  cat("  分组分布:\n")
  print(table(group_site))

  # ===== Step 2：标准化 =====
  cat("  标准化数据...\n")

  norm_data <- norm_abundance(abund_site, method = "TSS")

  # ===== Step 3：过滤低丰度（属水平阈值更低）=====

  # 确保顺序一致
  group_site <- group_site[colnames(norm_data)]

  # 计算每组平均丰度
  group_means <- do.call(
    cbind,
    lapply(levels(group_site), function(g) {
      samples_in_group <- names(group_site)[group_site == g]
      rowMeans(norm_data[, samples_in_group, drop = FALSE])
    })
  )
  colnames(group_means) <- levels(group_site)

  # 属水平过滤阈值：0.05%（更宽松）
  keep <- apply(group_means, 1, max) > 0.0005
  norm_data_filtered <- norm_data[keep, ]

  cat("  过滤后保留", nrow(norm_data_filtered), "个属\n")

  if (nrow(norm_data_filtered) < 3) {
    cat("  过滤后属数过少，跳过\n")
    next
  }

  # ===== Step 4：Dunn两两比较 =====
  cat("  进行Dunn两两比较检验...\n")

  results_dunn <- tryCatch(
    {
      diff_abundance(
        abund_table = norm_data_filtered,
        group = group_site,
        method = "KW_dunn",
        p_adjust_method = "bh"
      )
    },
    error = function(e) {
      cat("  错误:", e$message, "\n")
      NULL
    }
  )

  if (is.null(results_dunn) || nrow(results_dunn) == 0) {
    cat("  Dunn检验失败，跳过\n")
    next
  }

  # 保存完整结果
  write.csv(
    results_dunn,
    paste0("analysis_results/differential/", site, "_Genus_Dunn_all.csv"),
    row.names = FALSE
  )

  cat("  检验完成，共", nrow(results_dunn), "个两两比较\n")

  # ===== Step 5：筛选显著差异（放宽标准）=====

  # P.adj < 0.1（更宽松）, 效应量 r > 0.05（更宽松）
  sig_results <- results_dunn %>%
    dplyr::filter(!is.na(P.adj), P.adj < 0.1, abs(EffectSize) > 0.05) %>%
    dplyr::arrange(Comparison, desc(abs(EffectSize)))

  cat("  显著差异比较数:", nrow(sig_results), "\n")

  if (nrow(sig_results) > 0) {
    # 保存显著结果
    write.csv(
      sig_results,
      paste0(
        "analysis_results/differential/",
        site,
        "_Genus_Dunn_significant.csv"
      ),
      row.names = FALSE
    )

    # 提取显著差异的属
    sig_features <- unique(sig_results$Feature)
    cat("  显著差异属数:", length(sig_features), "\n")

    # 保存候选属丰度
    candidate_abund <- norm_data_filtered[sig_features, ]

    # 保存完整数据包
    candidate_data <- list(
      abundance = candidate_abund,
      metadata = meta_site,
      group = group_site,
      site = site,
      level = "Genus"
    )

    saveRDS(
      candidate_data,
      paste0(
        "analysis_results/differential/",
        site,
        "_Genus_candidates_data.rds"
      )
    )

    # 保存CSV
    write.csv(
      candidate_abund,
      paste0(
        "analysis_results/differential/",
        site,
        "_Genus_candidates_abundance.csv"
      )
    )

    # ===== 生成摘要统计 =====

    comparison_summary <- sig_results %>%
      dplyr::group_by(Comparison) %>%
      dplyr::summarise(
        N_significant = dplyr::n(),
        Mean_EffectSize = mean(abs(EffectSize), na.rm = TRUE),
        Max_EffectSize = max(abs(EffectSize), na.rm = TRUE),
        .groups = 'drop'
      )

    cat("\n  各比较的显著属数:\n")
    print(comparison_summary)

    write.csv(
      comparison_summary,
      paste0(
        "analysis_results/differential/",
        site,
        "_Genus_comparison_summary.csv"
      ),
      row.names = FALSE
    )

    # ===== 治疗效应统计 =====

    # 特别关注 DB_Control vs DB_Treated
    treatment_effect <- sig_results %>%
      dplyr::filter(grepl(
        "DB_Control.*DB_Treated|DB_Treated.*DB_Control",
        Comparison
      ))

    if (nrow(treatment_effect) > 0) {
      cat("\n  治疗相关的显著差异属数:", nrow(treatment_effect), "\n")

      write.csv(
        treatment_effect,
        paste0(
          "analysis_results/differential/",
          site,
          "_Genus_treatment_effect.csv"
        ),
        row.names = FALSE
      )
    }
  } else {
    cat("  无显著差异\n")
  }

  cat("  完成！\n")
}

cat("\n\n差异分析完成！\n")
cat("结果保存在: analysis_results/differential/\n")


# ============================================================
# 差异分析结果可视化：展示候选属的变化模式
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 创建可视化输出文件夹
if (!dir.exists("analysis_results/differential/figures")) {
  dir.create("analysis_results/differential/figures", recursive = TRUE)
}

# 定义清新配色
fresh_colors <- c(
  "WT_Control" = "#8DD3C7",
  "DB_Control" = "#BEBADA",
  "DB_Treated" = "#FB8072"
)


# 为每个Site生成可视化
for (site in sites) {
  cat(paste("\n========== 可视化", site, "==========\n"))

  # 读取候选数据
  rds_file <- paste0(
    "analysis_results/differential/",
    site,
    "_Genus_candidates_data.rds"
  )

  if (!file.exists(rds_file)) {
    cat("  未找到候选数据，跳过\n")
    next
  }

  candidate_data <- readRDS(rds_file)

  # 提取数据
  abund <- candidate_data$abundance
  metadata <- candidate_data$metadata
  group <- candidate_data$group

  if (nrow(abund) == 0) {
    cat("  无候选属，跳过\n")
    next
  }

  cat("  候选属数:", nrow(abund), "\n")

  # ===== 1. 热图：展示候选属在三组中的丰度模式 =====
  cat("  生成热图...\n")

  # 计算z-score标准化
  abund_zscore <- t(scale(t(abund)))

  # 按组排序样本
  sample_order <- order(group)
  abund_zscore_ordered <- abund_zscore[, sample_order]
  group_ordered <- group[sample_order]

  # 计算每个属的平均轨迹（用于聚类）
  group_means <- do.call(
    cbind,
    lapply(levels(group), function(g) {
      samples_in_group <- names(group)[group == g]
      rowMeans(abund[, samples_in_group, drop = FALSE])
    })
  )
  colnames(group_means) <- levels(group)

  # 对属进行层次聚类
  if (nrow(abund_zscore_ordered) > 2) {
    row_order <- hclust(dist(abund_zscore_ordered))$order
    abund_zscore_ordered <- abund_zscore_ordered[row_order, ]
  }

  # 设置Spectral颜色（反转使红色表示高值）
  col_fun <- colorRamp2(
    c(-2, 0, 2),
    rev(brewer.pal(11, "Spectral")[c(2, 6, 10)]) # 蓝-白-红
  )

  # 列注释（分组）
  col_anno <- HeatmapAnnotation(
    Group = group_ordered,
    col = list(Group = fresh_colors),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Group = list(title = "Group", nrow = 1)
    )
  )

  # 计算cell尺寸（正方形）
  n_rows <- nrow(abund_zscore_ordered)
  n_cols <- ncol(abund_zscore_ordered)
  cell_size <- unit(5, "mm") # 每个cell 5mm

  # 绘制热图
  pdf(
    paste0("analysis_results/differential/figures/", site, "_heatmap.pdf"),
    width = max(8, n_cols * 0.2 + 3),
    height = max(6, n_rows * 0.2 + 2)
  )

  ht <- Heatmap(
    abund_zscore_ordered,
    name = "Z-score",
    col = col_fun,
    top_annotation = col_anno,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_title = paste(site, "- Candidate Genera (n =", nrow(abund), ")"),

    # 设置cell为正方形 + 浅边框
    width = ncol(abund_zscore_ordered) * cell_size,
    height = nrow(abund_zscore_ordered) * cell_size,
    rect_gp = gpar(col = "grey90", lwd = 0.5), # 浅灰色边框，细线

    heatmap_legend_param = list(
      title = "Abundance\nZ-score",
      direction = "horizontal"
    )
  )

  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()

  # ===== 2. 轨迹图：展示每个属从WT→DB→Treated的变化 =====
  cat("  生成轨迹图...\n")

  # 确保group和abund列名对应
  group <- group[colnames(abund)]

  # 准备轨迹数据
  trajectory_data <- lapply(rownames(abund), function(genus) {
    means <- sapply(levels(group), function(g) {
      # 找到该组的样本
      samples_in_group <- names(group)[group == g]
      # 提取该属在这些样本中的丰度，转为数值向量
      values <- as.numeric(abund[genus, samples_in_group])
      mean(values, na.rm = TRUE)
    })

    data.frame(
      Genus = genus,
      Group = factor(levels(group), levels = levels(group)),
      Mean = means,
      stringsAsFactors = FALSE
    )
  })

  trajectory_df <- do.call(rbind, trajectory_data)

  # 检查
  cat("  轨迹数据前几行:\n")
  print(head(trajectory_df))

  # 计算变化模式类型
  trajectory_types <- sapply(rownames(abund), function(genus) {
    # 提取三组的平均值
    wt_samples <- names(group)[group == "WT_Control"]
    db_samples <- names(group)[group == "DB_Control"]
    tr_samples <- names(group)[group == "DB_Treated"]

    wt <- mean(as.numeric(abund[genus, wt_samples]), na.rm = TRUE)
    db <- mean(as.numeric(abund[genus, db_samples]), na.rm = TRUE)
    tr <- mean(as.numeric(abund[genus, tr_samples]), na.rm = TRUE)

    # 判断模式
    if (is.na(wt) || is.na(db) || is.na(tr)) {
      return("No Data")
    }

    if (db > wt && tr < db && tr > wt) {
      "Recovery (U-shape)"
    } else if (db > wt && tr <= wt) {
      "Full Recovery"
    } else if (db < wt && tr > db && tr < wt) {
      "Recovery (Inverted-U)"
    } else if (db < wt && tr >= wt) {
      "Full Recovery"
    } else if (db > wt && tr > db) {
      "Worsening"
    } else if (db < wt && tr < db) {
      "Worsening"
    } else {
      "No Clear Pattern"
    }
  })

  trajectory_df$Pattern <- trajectory_types[trajectory_df$Genus]

  # 统计各模式数量
  pattern_summary <- table(trajectory_types)
  cat("\n  变化模式统计:\n")
  print(pattern_summary)

  # 按模式分面绘制轨迹图
  p_trajectory <- ggplot(
    trajectory_df,
    aes(x = Group, y = Mean, group = Genus)
  ) +
    geom_line(aes(color = Pattern), alpha = 0.6, linewidth = 0.8) +
    geom_point(aes(color = Pattern), size = 2, alpha = 0.7) +
    facet_wrap(~Pattern, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Set2") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(
      title = paste(site, "- Trajectory Patterns of Candidate Genera"),
      x = "",
      y = "Mean Relative Abundance",
      subtitle = paste("Total candidates:", nrow(abund))
    )

  ggsave(
    paste0(
      "analysis_results/differential/figures/",
      site,
      "_trajectory_patterns.pdf"
    ),
    p_trajectory,
    width = 12,
    height = max(6, ceiling(length(unique(trajectory_types)) / 2) * 3)
  )

  # ===== 3. 平均轨迹图（按模式分组）=====
  cat("  生成平均轨迹图...\n")

  trajectory_summary <- trajectory_df %>%
    group_by(Pattern, Group) %>%
    summarise(
      Mean = mean(Mean),
      SE = sd(Mean) / sqrt(n()),
      .groups = 'drop'
    )

  p_trajectory_mean <- ggplot(
    trajectory_summary,
    aes(x = Group, y = Mean, color = Pattern, group = Pattern)
  ) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    scale_color_brewer(palette = "Set2") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.8)
    ) +
    labs(
      title = paste(site, "- Average Trajectory by Pattern"),
      x = "",
      y = "Mean Relative Abundance",
      color = "Pattern"
    )

  ggsave(
    paste0(
      "analysis_results/differential/figures/",
      site,
      "_trajectory_mean.pdf"
    ),
    p_trajectory_mean,
    width = 8,
    height = 6
  )

  # ===== 4. Top候选属的箱线图 =====
  cat("  生成Top候选箱线图...\n")

  # 读取显著性结果，找出效应量最大的
  sig_file <- paste0(
    "analysis_results/differential/",
    site,
    "_Genus_Dunn_significant.csv"
  )

  if (file.exists(sig_file)) {
    sig_results <- read.csv(sig_file)

    # 提取治疗相关比较
    treatment_comparison <- sig_results %>%
      filter(grepl(
        "DB_Control.*DB_Treated|DB_Treated.*DB_Control",
        Comparison
      )) %>%
      arrange(desc(abs(EffectSize))) %>%
      head(12) # Top 12

    if (nrow(treatment_comparison) > 0) {
      top_genera <- unique(treatment_comparison$Feature)

      # 准备箱线图数据
      top_data <- abund[top_genera, , drop = FALSE]
      top_long <- top_data %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Genus") %>%
        pivot_longer(-Genus, names_to = "Sample", values_to = "Abundance") %>%
        left_join(
          metadata %>% select(SampleID_converted, TreatGroup),
          by = c("Sample" = "SampleID_converted")
        )

      top_long$TreatGroup <- factor(
        top_long$TreatGroup,
        levels = c("WT_Control", "DB_Control", "DB_Treated")
      )

      # 绘制
      p_boxplot <- ggplot(
        top_long,
        aes(x = TreatGroup, y = Abundance, fill = TreatGroup)
      ) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        facet_wrap(~Genus, scales = "free_y", ncol = 3) +
        scale_fill_manual(values = fresh_colors) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90")
        ) +
        labs(
          title = paste(
            site,
            "- Top Candidate Genera (by Treatment Effect Size)"
          ),
          x = "",
          y = "Relative Abundance",
          fill = "Group"
        )

      ggsave(
        paste0(
          "analysis_results/differential/figures/",
          site,
          "_top_candidates_boxplot.pdf"
        ),
        p_boxplot,
        width = 8,
        height = max(6, ceiling(length(top_genera) / 3) * 3)
      )
    }
  }

  # ===== 5. 变化模式统计柱状图 =====
  cat("  生成模式统计图...\n")

  pattern_df <- data.frame(
    Pattern = names(pattern_summary),
    Count = as.numeric(pattern_summary)
  )

  p_pattern_bar <- ggplot(
    pattern_df,
    aes(x = reorder(Pattern, Count), y = Count, fill = Pattern)
  ) +
    geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.5) +
    geom_text(aes(label = Count), hjust = -0.2, size = 4) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.8)
    ) +
    labs(
      title = paste(site, "- Distribution of Change Patterns"),
      x = "",
      y = "Number of Genera"
    )

  ggsave(
    paste0(
      "analysis_results/differential/figures/",
      site,
      "_pattern_distribution.pdf"
    ),
    p_pattern_bar,
    width = 8,
    height = 5
  )

  # ===== 保存模式分类结果 =====
  pattern_classification <- data.frame(
    Genus = rownames(abund),
    Pattern = trajectory_types,
    Mean_WT = group_means[, "WT_Control"],
    Mean_DB = group_means[, "DB_Control"],
    Mean_Treated = group_means[, "DB_Treated"],
    stringsAsFactors = FALSE
  )

  write.csv(
    pattern_classification,
    paste0(
      "analysis_results/differential/",
      site,
      "_pattern_classification.csv"
    ),
    row.names = FALSE
  )

  cat("  完成！\n")
}

cat("\n\n可视化完成！\n")
cat("图表保存在: analysis_results/differential/figures/\n")


# ============================================================
# 各部位候选属的网络分析
# ============================================================

library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)

# 创建输出文件夹
if (!dir.exists("analysis_results/networks")) {
  dir.create("analysis_results/networks", recursive = TRUE)
}

# 定义变化模式配色
pattern_colors <- c(
  "Full Recovery" = "#2ecc71",
  "Recovery (U-shape)" = "#3498db",
  "Recovery (Inverted-U)" = "#9b59b6",
  "Worsening" = "#e74c3c",
  "No Clear Pattern" = "#95a5a6"
)

# ===== 为每个Site构建网络 =====

for (site in sites) {
  cat(paste("\n========== 构建", site, "网络 ==========\n"))

  # 读取数据
  rds_file <- paste0(
    "analysis_results/differential/",
    site,
    "_Genus_candidates_data.rds"
  )
  pattern_file <- paste0(
    "analysis_results/differential/",
    site,
    "_pattern_classification.csv"
  )
  sig_file <- paste0(
    "analysis_results/differential/",
    site,
    "_Genus_Dunn_significant.csv"
  )

  if (!file.exists(rds_file) || !file.exists(pattern_file)) {
    cat("  数据文件缺失，跳过\n")
    next
  }

  candidate_data <- readRDS(rds_file)
  patterns <- read.csv(pattern_file)
  sig_results <- read.csv(sig_file)

  abund <- candidate_data$abundance
  group <- candidate_data$group

  if (nrow(abund) < 3) {
    cat("  候选属过少（<3），跳过\n")
    next
  }

  cat("  候选属数:", nrow(abund), "\n")

  # ===== 1. 计算属间相关性（Spearman）=====
  cat("  计算相关性矩阵...\n")

  # 转置：行=样本，列=属
  abund_t <- t(abund)

  # 计算Spearman相关
  cor_matrix <- cor(abund_t, method = "spearman", use = "pairwise.complete.obs")

  # 计算p值
  cor_test_results <- matrix(NA, nrow = ncol(abund_t), ncol = ncol(abund_t))
  rownames(cor_test_results) <- colnames(cor_test_results) <- colnames(abund_t)

  for (i in 1:(ncol(abund_t) - 1)) {
    for (j in (i + 1):ncol(abund_t)) {
      test <- cor.test(
        abund_t[, i],
        abund_t[, j],
        method = "spearman",
        exact = FALSE
      )
      cor_test_results[i, j] <- test$p.value
      cor_test_results[j, i] <- test$p.value
    }
  }

  # FDR校正
  p_values <- cor_test_results[upper.tri(cor_test_results)]
  p_adjusted <- p.adjust(p_values, method = "fdr")

  cor_test_adj <- cor_test_results
  cor_test_adj[upper.tri(cor_test_adj)] <- p_adjusted
  cor_test_adj[lower.tri(cor_test_adj)] <- t(cor_test_adj)[lower.tri(
    cor_test_adj
  )]

  # ===== 2. 构建边表（显著相关 + 相关性阈值）=====
  cat("  构建网络边...\n")

  edge_list <- data.frame()

  for (i in 1:(nrow(cor_matrix) - 1)) {
    for (j in (i + 1):nrow(cor_matrix)) {
      cor_val <- cor_matrix[i, j]
      p_val <- cor_test_adj[i, j]

      # 筛选条件：|r| > 0.5 且 P.adj < 0.05
      if (
        !is.na(cor_val) && !is.na(p_val) && abs(cor_val) > 0.5 && p_val < 0.05
      ) {
        edge_list <- rbind(
          edge_list,
          data.frame(
            from = rownames(cor_matrix)[i],
            to = rownames(cor_matrix)[j],
            weight = abs(cor_val),
            correlation = cor_val,
            p_value = p_val,
            edge_type = ifelse(cor_val > 0, "Positive", "Negative"),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if (nrow(edge_list) == 0) {
    cat("  无显著相关性，跳过网络绘制\n")
    next
  }

  cat("  显著相关边数:", nrow(edge_list), "\n")
  cat("    - 正相关:", sum(edge_list$edge_type == "Positive"), "\n")
  cat("    - 负相关:", sum(edge_list$edge_type == "Negative"), "\n")

  # 保存边表
  write.csv(
    edge_list,
    paste0("analysis_results/networks/", site, "_edges.csv"),
    row.names = FALSE
  )

  # ===== 3. 构建节点属性 =====
  cat("  构建节点属性...\n")

  # 提取网络中的节点
  nodes_in_network <- unique(c(edge_list$from, edge_list$to))

  node_data <- data.frame(
    Genus = nodes_in_network,
    stringsAsFactors = FALSE
  )

  # 添加变化模式
  node_data <- node_data %>%
    left_join(patterns %>% select(Genus, Pattern), by = "Genus")

  # 添加平均丰度（log10转换）
  node_data$Mean_abundance <- sapply(node_data$Genus, function(g) {
    mean(as.numeric(abund[g, ]), na.rm = TRUE)
  })

  node_data$Log_abundance <- log10(node_data$Mean_abundance + 1e-6)

  # 添加度中心性（先构建临时图）
  temp_graph <- graph_from_data_frame(
    edge_list,
    directed = FALSE,
    vertices = node_data
  )
  node_data$Degree <- degree(temp_graph)

  # 添加效应量（从显著性结果中提取治疗相关的）
  treatment_effects <- sig_results %>%
    filter(grepl(
      "DB_Control.*DB_Treated|DB_Treated.*DB_Control",
      Comparison
    )) %>%
    group_by(Feature) %>%
    summarise(
      Max_EffectSize = max(abs(EffectSize), na.rm = TRUE),
      .groups = 'drop'
    )

  node_data <- node_data %>%
    left_join(treatment_effects, by = c("Genus" = "Feature"))

  node_data$Max_EffectSize[is.na(node_data$Max_EffectSize)] <- 0

  # 保存节点表
  write.csv(
    node_data,
    paste0("analysis_results/networks/", site, "_nodes.csv"),
    row.names = FALSE
  )

  # ===== 4. 构建igraph对象 =====
  cat("  构建igraph对象...\n")

  g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = node_data)

  # ===== 5. 网络统计 =====
  cat("  计算网络统计...\n")

  network_stats <- data.frame(
    Site = site,
    N_nodes = vcount(g),
    N_edges = ecount(g),
    Density = edge_density(g),
    Avg_degree = mean(degree(g)),
    N_components = components(g)$no,
    Clustering_coef = transitivity(g, type = "global"),
    Avg_path_length = tryCatch(
      mean_distance(g, directed = FALSE),
      error = function(e) NA
    ),
    stringsAsFactors = FALSE
  )

  cat("\n  网络统计:\n")
  print(network_stats)

  write.csv(
    network_stats,
    paste0("analysis_results/networks/", site, "_network_stats.csv"),
    row.names = FALSE
  )

  # ===== 6. 可视化：ggraph网络图 =====
  cat("  绘制网络图...\n")

  # 转换为tidygraph对象
  tg <- as_tbl_graph(g)

  # 布局算法：Fruchterman-Reingold
  set.seed(123)

  # 网络图1：按变化模式着色
  p_network <- ggraph(tg, layout = "fr") +
    # 绘制边
    geom_edge_link(aes(
      edge_width = weight,
      color = edge_type,
      alpha = weight
    )) +
    scale_edge_width_continuous(range = c(0.3, 2), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    scale_edge_color_manual(
      values = c("Positive" = "#2ecc71", "Negative" = "#e74c3c"),
      name = "Correlation"
    ) +
    # 绘制节点
    geom_node_point(
      aes(size = Degree, fill = Pattern),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_size_continuous(range = c(3, 10), name = "Degree") +
    scale_fill_manual(values = pattern_colors, name = "Change Pattern") +
    # 添加标签
    geom_node_text(
      aes(label = name),
      size = 2.5,
      repel = TRUE,
      max.overlaps = 20,
      bg.color = "white",
      bg.r = 0.1
    ) +
    theme_graph() +
    labs(
      title = paste(site, "- Co-abundance Network"),
      subtitle = paste0(
        "Nodes: ",
        vcount(g),
        " | Edges: ",
        ecount(g),
        " | Density: ",
        round(edge_density(g), 3)
      )
    ) +
    theme(legend.position = "right")

  ggsave(
    paste0("analysis_results/networks/", site, "_network.pdf"),
    p_network,
    width = 8,
    height = 6,
    device = cairo_pdf
  )

  # 网络图2：节点大小反映效应量
  p_network_effect <- ggraph(tg, layout = "fr") +
    geom_edge_link(aes(
      edge_width = weight,
      color = edge_type,
      alpha = weight
    )) +
    scale_edge_width_continuous(range = c(0.3, 2), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    scale_edge_color_manual(
      values = c("Positive" = "#2ecc71", "Negative" = "#e74c3c"),
      name = "Correlation"
    ) +
    geom_node_point(
      aes(size = Max_EffectSize, fill = Pattern),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_size_continuous(range = c(3, 10), name = "Effect Size") +
    scale_fill_manual(values = pattern_colors, name = "Change Pattern") +
    geom_node_text(
      aes(label = name),
      size = 2.5,
      repel = TRUE,
      max.overlaps = 20,
      bg.color = "white",
      bg.r = 0.1
    ) +
    theme_graph() +
    labs(
      title = paste(site, "- Co-abundance Network (Effect Size)"),
      subtitle = paste0("Node size = Treatment effect size")
    ) +
    theme(legend.position = "right")

  ggsave(
    paste0("analysis_results/networks/", site, "_network_effectsize.pdf"),
    p_network_effect,
    width = 8,
    height = 6,
    device = cairo_pdf
  )

  # ===== 7. 识别关键节点（Hub）=====
  cat("  识别Hub节点...\n")

  # 计算中心性指标
  centrality <- data.frame(
    Genus = V(g)$name,
    Degree = degree(g),
    Betweenness = betweenness(g),
    Closeness = closeness(g),
    Eigenvector = eigen_centrality(g)$vector,
    stringsAsFactors = FALSE
  )

  # 标准化
  centrality$Degree_scaled <- scale(centrality$Degree)[, 1]
  centrality$Betweenness_scaled <- scale(centrality$Betweenness)[, 1]
  centrality$Closeness_scaled <- scale(centrality$Closeness)[, 1]
  centrality$Eigenvector_scaled <- scale(centrality$Eigenvector)[, 1]

  # Hub定义：度中心性 > 均值 + 1*SD
  degree_threshold <- mean(centrality$Degree) + sd(centrality$Degree)
  centrality$Is_hub <- centrality$Degree > degree_threshold

  hubs <- centrality %>%
    filter(Is_hub) %>%
    arrange(desc(Degree))

  cat("  Hub节点数:", nrow(hubs), "\n")
  if (nrow(hubs) > 0) {
    cat("  Top Hub:\n")
    print(head(hubs %>% select(Genus, Degree, Betweenness), 5))
  }

  write.csv(
    centrality,
    paste0("analysis_results/networks/", site, "_centrality.csv"),
    row.names = FALSE
  )

  # ===== 8. 社群检测 =====
  cat("  检测网络社群...\n")

  # Louvain算法
  communities <- cluster_louvain(g)

  community_info <- data.frame(
    Genus = V(g)$name,
    Community = membership(communities),
    stringsAsFactors = FALSE
  )

  # 统计每个社群的大小和模式
  community_summary <- community_info %>%
    left_join(node_data %>% select(Genus, Pattern), by = "Genus") %>%
    group_by(Community) %>%
    summarise(
      Size = n(),
      Dominant_pattern = names(which.max(table(Pattern))),
      .groups = 'drop'
    ) %>%
    arrange(desc(Size))

  cat("  检测到", length(unique(community_info$Community)), "个社群\n")

  write.csv(
    community_info,
    paste0("analysis_results/networks/", site, "_communities.csv"),
    row.names = FALSE
  )

  write.csv(
    community_summary,
    paste0("analysis_results/networks/", site, "_community_summary.csv"),
    row.names = FALSE
  )

  # 将社群信息添加到图对象
  V(g)$community <- community_info$Community

  # 重新创建tidygraph对象（包含community信息）
  tg_community <- as_tbl_graph(g)

  # 绘制社群网络
  p_network_community <- ggraph(tg_community, layout = "fr") +
    geom_edge_link(aes(edge_width = weight, alpha = weight), color = "grey70") +
    scale_edge_width_continuous(range = c(0.3, 2), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    geom_node_point(
      aes(size = Degree, fill = factor(community)), # 现在可以访问community了
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_size_continuous(range = c(3, 10), name = "Degree") +
    scale_fill_brewer(palette = "Set3", name = "Community") +
    geom_node_text(
      aes(label = name),
      size = 2.5,
      repel = TRUE,
      max.overlaps = 20
    ) +
    theme_graph() +
    labs(
      title = paste(site, "- Network Communities"),
      subtitle = paste0(
        "Detected ",
        length(unique(community_info$Community)),
        " communities"
      )
    ) +
    theme(legend.position = "right")
  ggsave(
    paste0("analysis_results/networks/", site, "_network_communities.pdf"),
    p_network_community,
    width = 8,
    height = 6,
    device = cairo_pdf
  )

  cat("  完成！\n")
}

# ===== 汇总所有Site的网络统计 =====
cat("\n===== 汇总网络统计 =====\n")

all_network_stats <- data.frame()

for (site in sites) {
  stats_file <- paste0("analysis_results/networks/", site, "_network_stats.csv")
  if (file.exists(stats_file)) {
    stats <- read.csv(stats_file)
    all_network_stats <- rbind(all_network_stats, stats)
  }
}

if (nrow(all_network_stats) > 0) {
  write.csv(
    all_network_stats,
    "analysis_results/networks/all_sites_network_stats.csv",
    row.names = FALSE
  )

  cat("\n所有Site的网络统计:\n")
  print(all_network_stats)
}

cat("\n网络分析完成！\n")
cat("结果保存在: analysis_results/networks/\n")


# ============================================================
# 改进版网络分析：放宽标准 + 分组网络
# ============================================================

library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)

if (!dir.exists("analysis_results/networks_improved")) {
  dir.create("analysis_results/networks_improved", recursive = TRUE)
}

pattern_colors <- c(
  "Full Recovery" = "#2ecc71",
  "Recovery (U-shape)" = "#3498db",
  "Recovery (Inverted-U)" = "#9b59b6",
  "Worsening" = "#e74c3c",
  "No Clear Pattern" = "#95a5a6"
)

# ===== 改进的网络构建（更多节点）=====

for (site in sites) {
  cat(paste("\n========== 分析", site, "==========\n"))

  # 读取该Site的所有属丰度（不限于显著差异）
  site_samples <- metadata_filtered %>%
    filter(
      Site == site,
      TreatGroup %in% c("WT_Control", "DB_Control", "DB_Treated")
    ) %>%
    pull(SampleID_converted)

  if (length(site_samples) < 6) {
    cat("  样本数不足\n")
    next
  }

  # 读取属水平数据
  abund_site <- genus_abs_filtered[, site_samples]

  # 标准化
  norm_data <- norm_abundance(abund_site, method = "TSS")

  # 获取分组
  meta_site <- metadata_filtered %>%
    filter(SampleID_converted %in% site_samples)

  group_site <- factor(
    meta_site$TreatGroup,
    levels = c("WT_Control", "DB_Control", "DB_Treated")
  )
  names(group_site) <- meta_site$SampleID_converted

  # ===== 放宽过滤标准：只保留平均丰度 > 0.01%（更宽松）=====
  mean_abund <- rowMeans(norm_data)
  keep <- mean_abund > 0.005 # 0.01%

  norm_data_filtered <- norm_data[keep, ]

  cat("  过滤后属数:", nrow(norm_data_filtered), "\n")

  if (nrow(norm_data_filtered) < 5) {
    cat("  属数过少\n")
    next
  }

  # ===== 为每个组单独构建网络 =====

  group_levels <- c("WT_Control", "DB_Control", "DB_Treated")

  for (grp in group_levels) {
    cat(paste("\n  --- 构建", grp, "网络 ---\n"))

    # 筛选该组的样本
    grp_samples <- names(group_site)[group_site == grp]

    if (length(grp_samples) < 3) {
      cat("    样本数不足(<3)\n")
      next
    }

    grp_abund <- norm_data_filtered[, grp_samples]

    # 进一步过滤：在该组中至少在60%样本中存在
    prevalence <- rowSums(grp_abund > 0) / ncol(grp_abund)
    keep_grp <- prevalence >= 0.6

    grp_abund_filtered <- grp_abund[keep_grp, ]

    cat("    该组属数:", nrow(grp_abund_filtered), "\n")

    if (nrow(grp_abund_filtered) < 5) {
      cat("    属数过少\n")
      next
    }

    # ===== 计算相关性（Spearman）=====
    abund_t <- t(grp_abund_filtered)

    cor_matrix <- cor(
      abund_t,
      method = "spearman",
      use = "pairwise.complete.obs"
    )

    # 计算p值
    n_genera <- ncol(abund_t)
    cor_test_results <- matrix(1, nrow = n_genera, ncol = n_genera)
    rownames(cor_test_results) <- colnames(cor_test_results) <- colnames(
      abund_t
    )

    for (i in 1:(n_genera - 1)) {
      for (j in (i + 1):n_genera) {
        if (sum(!is.na(abund_t[, i]) & !is.na(abund_t[, j])) >= 3) {
          test <- cor.test(
            abund_t[, i],
            abund_t[, j],
            method = "spearman",
            exact = FALSE
          )
          cor_test_results[i, j] <- test$p.value
          cor_test_results[j, i] <- test$p.value
        }
      }
    }

    # ===== 构建边表（放宽标准）=====
    edge_list <- data.frame()

    for (i in 1:(nrow(cor_matrix) - 1)) {
      for (j in (i + 1):nrow(cor_matrix)) {
        cor_val <- cor_matrix[i, j]
        p_val <- cor_test_results[i, j]

        # 放宽标准：|r| > 0.4 且 P < 0.05（不校正）
        if (
          !is.na(cor_val) && !is.na(p_val) && abs(cor_val) > 0.6 && p_val < 0.05
        ) {
          edge_list <- rbind(
            edge_list,
            data.frame(
              from = rownames(cor_matrix)[i],
              to = rownames(cor_matrix)[j],
              weight = abs(cor_val),
              correlation = cor_val,
              p_value = p_val,
              edge_type = ifelse(cor_val > 0, "Positive", "Negative"),
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }

    if (nrow(edge_list) == 0) {
      cat("    无显著相关\n")
      next
    }

    cat("    边数:", nrow(edge_list), "\n")

    # 保存边表
    write.csv(
      edge_list,
      paste0(
        "analysis_results/networks_improved/",
        site,
        "_",
        grp,
        "_edges.csv"
      ),
      row.names = FALSE
    )

    # ===== 构建节点属性 =====
    nodes_in_network <- unique(c(edge_list$from, edge_list$to))

    node_data <- data.frame(
      Genus = nodes_in_network,
      Mean_abundance = sapply(nodes_in_network, function(g) {
        mean(as.numeric(grp_abund_filtered[g, ]), na.rm = TRUE)
      }),
      stringsAsFactors = FALSE
    )

    node_data$Log_abundance <- log10(node_data$Mean_abundance + 1e-6)

    # 先构建图
    g <- graph_from_data_frame(
      edge_list,
      directed = FALSE,
      vertices = node_data
    )

    # 计算度并添加到node_data
    node_data$Degree <- degree(g)

    # 将Degree添加到图对象的顶点属性
    V(g)$Degree <- node_data$Degree
    V(g)$Mean_abundance <- node_data$Mean_abundance
    V(g)$Log_abundance <- node_data$Log_abundance

    # 网络统计
    network_stats <- data.frame(
      Site = site,
      Group = grp,
      N_nodes = vcount(g),
      N_edges = ecount(g),
      Density = edge_density(g),
      Avg_degree = mean(degree(g)),
      Clustering_coef = transitivity(g, type = "global"),
      stringsAsFactors = FALSE
    )

    cat(
      "    节点:",
      vcount(g),
      "| 边:",
      ecount(g),
      "| 密度:",
      round(edge_density(g), 3),
      "\n"
    )

    write.csv(
      network_stats,
      paste0(
        "analysis_results/networks_improved/",
        site,
        "_",
        grp,
        "_stats.csv"
      ),
      row.names = FALSE
    )

    # ===== 绘制网络 =====
    # 现在创建tidygraph对象，Degree已经在图中了
    tg <- as_tbl_graph(g)

    set.seed(123)

    p_network <- ggraph(tg, layout = "fr") +
      geom_edge_link(aes(
        edge_width = weight,
        color = edge_type,
        alpha = weight
      )) +
      scale_edge_width_continuous(range = c(0.3, 2), guide = "none") +
      scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
      scale_edge_color_manual(
        values = c("Positive" = "#2ecc71", "Negative" = "#e74c3c"),
        name = "Correlation"
      ) +
      geom_node_point(
        aes(size = Degree),
        fill = "#3498db", # 现在可以访问了
        shape = 21,
        color = "black",
        stroke = 0.5
      ) +
      scale_size_continuous(range = c(3, 10), name = "Degree") +
      geom_node_text(
        aes(label = name),
        size = 2.5,
        repel = TRUE,
        max.overlaps = 20
      ) +
      theme_graph() +
      labs(
        title = paste(site, "-", grp, "Network"),
        subtitle = paste0(
          "Nodes: ",
          vcount(g),
          " | Edges: ",
          ecount(g),
          " | Density: ",
          round(edge_density(g), 3)
        )
      ) +
      theme(legend.position = "right")

    ggsave(
      paste0(
        "analysis_results/networks_improved/",
        site,
        "_",
        grp,
        "_network.pdf"
      ),
      p_network,
      width = 8,
      height = 8,
      device = cairo_pdf,
    )
  }
}

# ===== 网络比较分析 =====
cat("\n===== 网络比较分析 =====\n")

all_stats <- data.frame()

for (site in sites) {
  for (grp in c("WT_Control", "DB_Control", "DB_Treated")) {
    stats_file <- paste0(
      "analysis_results/networks_improved/",
      site,
      "_",
      grp,
      "_stats.csv"
    )
    if (file.exists(stats_file)) {
      stats <- read.csv(stats_file)
      all_stats <- rbind(all_stats, stats)
    }
  }
}

if (nrow(all_stats) > 0) {
  write.csv(
    all_stats,
    "analysis_results/networks_improved/all_network_stats_comparison.csv",
    row.names = FALSE
  )

  cat("\n网络统计比较:\n")
  print(all_stats)

  # ===== 可视化：网络指标比较 =====

  library(tidyr)

  # 转为长格式
  stats_long <- all_stats %>%
    pivot_longer(
      cols = c(N_nodes, N_edges, Density, Avg_degree, Clustering_coef),
      names_to = "Metric",
      values_to = "Value"
    )

  stats_long$Group <- factor(
    stats_long$Group,
    levels = c("WT_Control", "DB_Control", "DB_Treated")
  )

  # 绘制比较图
  p_compare <- ggplot(stats_long, aes(x = Group, y = Value, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.5) +
    facet_grid(Metric ~ Site, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "WT_Control" = "#8DD3C7",
        "DB_Control" = "#FFFFB3",
        "DB_Treated" = "#BEBADA"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(
      title = "Network Topology Comparison Across Groups",
      x = "",
      y = "Value",
      fill = "Group"
    )

  ggsave(
    "analysis_results/networks_improved/network_comparison.pdf",
    p_compare,
    width = 12,
    height = 10,
    device = cairo_pdf,
  )
}

cat("\n改进版网络分析完成！\n")
cat("结果保存在: analysis_results/networks_improved/\n")
