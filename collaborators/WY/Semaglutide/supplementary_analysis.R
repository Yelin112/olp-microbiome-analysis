# ============================================================
# 补充分析：热图和功能预测
# ============================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 设置工作路径
base_path <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
setwd(base_path)

# ============================================================
# 1. 热图分析 - 主要类群在不同组的分布模式
# ============================================================

# 读取genus水平的相对丰度
genus_rel <- read.table("data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.g.relative.xls",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 读取metadata
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, row.names = 1)
metadata$SampleID_converted <- gsub("_", ".", rownames(metadata))

# 获取共同样本
common_samples <- intersect(metadata$SampleID_converted, colnames(genus_rel))
genus_rel_filtered <- genus_rel[, common_samples]
metadata_filtered <- metadata[metadata$SampleID_converted %in% common_samples, ]

# 选择Top 30丰度的属
genus_mean <- rowMeans(genus_rel_filtered)
top_genera <- names(sort(genus_mean, decreasing = TRUE))[1:30]
genus_top <- genus_rel_filtered[top_genera, ]

# 对数转换（避免0值问题）
genus_top_log <- log10(genus_top + 1e-6)

# 准备注释
ha <- HeatmapAnnotation(
  Site = metadata_filtered[colnames(genus_top), "Site"],
  Condition = metadata_filtered[colnames(genus_top), "Condition"],
  Treatment = metadata_filtered[colnames(genus_top), "Treatment"],
  col = list(
    Site = c("Oral" = "#E41A1C", "Gut" = "#377EB8", "Fecal" = "#4DAF4A", "Lung" = "#984EA3"),
    Condition = c("WT" = "#FF7F00", "DB" = "#FFFF33"),
    Treatment = c("Control" = "#A65628", "Semaglutide" = "#F781BF")
  ),
  annotation_name_side = "left"
)

# 绘制热图
pdf("analysis_results/figures/heatmap_top30_genera.pdf", width = 14, height = 10)
Heatmap(as.matrix(genus_top_log),
        name = "log10(Relative\nAbundance)",
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 6),
        col = colorRamp2(c(min(genus_top_log), 0, max(genus_top_log)), 
                        c("blue", "white", "red")))
dev.off()

# 按Site分别绘制热图
sites <- unique(metadata_filtered$Site)

for(site in sites) {
  site_samples <- metadata_filtered %>% 
    filter(Site == site) %>% 
    pull(SampleID_converted)
  
  if(length(site_samples) > 2) {
    genus_site <- genus_top[, site_samples]
    genus_site_log <- log10(genus_site + 1e-6)
    
    metadata_site <- metadata_filtered[site_samples, ]
    
    ha_site <- HeatmapAnnotation(
      Condition = metadata_site$Condition,
      Treatment = metadata_site$Treatment,
      col = list(
        Condition = c("WT" = "#FF7F00", "DB" = "#FFFF33"),
        Treatment = c("Control" = "#A65628", "Semaglutide" = "#F781BF")
      ),
      annotation_name_side = "left"
    )
    
    pdf(paste0("analysis_results/figures/heatmap_", site, "_genera.pdf"), 
        width = 10, height = 10)
    print(Heatmap(as.matrix(genus_site_log),
            name = "log10(RA)",
            top_annotation = ha_site,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            show_column_names = TRUE,
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            col = colorRamp2(c(min(genus_site_log), 0, max(genus_site_log)), 
                            c("blue", "white", "red")),
            column_title = paste("Top 30 Genera -", site)))
    dev.off()
  }
}

cat("\n热图分析完成！\n")

# ============================================================
# 2. 功能预测分析 (PICRUSt2)
# ============================================================

# 读取KO功能预测结果
ko_file <- "E:/semaglutide测序原始数据/扩增子/02.Result_X101SC22114623-Z01-J001_16S/result/05.FunctionPrediction/PICRUSt2/picrust2_out_pipeline/KO_predicted.tsv.gz"

if(file.exists(ko_file)) {
  ko_data <- read.table(gzfile(ko_file), header = TRUE, sep = "\t", 
                       row.names = 1, check.names = FALSE, comment.char = "")
  
  # 检查样本名
  ko_samples <- colnames(ko_data)
  common_ko_samples <- intersect(metadata$SampleID_converted, ko_samples)
  
  if(length(common_ko_samples) > 0) {
    ko_filtered <- ko_data[, common_ko_samples]
    
    # 计算每个KO的总丰度
    ko_sums <- rowSums(ko_filtered)
    
    # 选择Top 50丰度的KO
    top_ko <- names(sort(ko_sums, decreasing = TRUE))[1:50]
    ko_top <- ko_filtered[top_ko, ]
    
    # 归一化为相对丰度
    ko_rel <- sweep(ko_top, 2, colSums(ko_top), "/")
    
    # 对数转换
    ko_log <- log10(ko_rel + 1e-6)
    
    # 准备注释
    metadata_ko <- metadata_filtered[common_ko_samples, ]
    
    ha_ko <- HeatmapAnnotation(
      Site = metadata_ko$Site,
      Condition = metadata_ko$Condition,
      Treatment = metadata_ko$Treatment,
      col = list(
        Site = c("Oral" = "#E41A1C", "Gut" = "#377EB8", 
                "Fecal" = "#4DAF4A", "Lung" = "#984EA3"),
        Condition = c("WT" = "#FF7F00", "DB" = "#FFFF33"),
        Treatment = c("Control" = "#A65628", "Semaglutide" = "#F781BF")
      ),
      annotation_name_side = "left"
    )
    
    # 绘制KO热图
    pdf("analysis_results/figures/heatmap_top50_KO.pdf", width = 14, height = 12)
    Heatmap(as.matrix(ko_log),
            name = "log10(RA)",
            top_annotation = ha_ko,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            show_column_names = TRUE,
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 6),
            col = colorRamp2(c(min(ko_log), 0, max(ko_log)), 
                            c("blue", "white", "red")),
            column_title = "Top 50 KO Functions")
    dev.off()
    
    cat("\nKO功能预测分析完成！\n")
  } else {
    cat("\n警告：KO数据中没有找到匹配的样本名\n")
  }
} else {
  cat("\n警告：未找到KO预测文件\n")
}

# ============================================================
# 3. 特定类群的详细分析
# ============================================================

# 分析特定的感兴趣类群（如乳酸杆菌、双歧杆菌等）
interest_genera <- c("Lactobacillus", "Bifidobacterium", "Bacteroides", 
                    "Akkermansia", "Faecalibacterium")

genus_rel_long <- genus_rel_filtered %>%
  rownames_to_column("Genus") %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Abundance")

genus_rel_long <- genus_rel_long %>%
  left_join(metadata_filtered %>% 
             select(SampleID_converted, Site, Condition, Treatment) %>%
             rename(Sample = SampleID_converted), 
           by = "Sample")

# 筛选感兴趣的类群
interest_data <- genus_rel_long %>%
  filter(Genus %in% interest_genera)

if(nrow(interest_data) > 0) {
  # 为每个Site绘制感兴趣类群的丰度
  for(site in sites) {
    site_interest <- interest_data %>% filter(Site == site)
    
    if(nrow(site_interest) > 0) {
      site_interest$TreatGroup <- paste(site_interest$Condition, 
                                       site_interest$Treatment, sep = "_")
      
      p <- ggplot(site_interest, aes(x = TreatGroup, y = Abundance, fill = TreatGroup)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        facet_wrap(~ Genus, scales = "free_y", ncol = 2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(title = paste("感兴趣类群丰度 -", site),
             x = "Group", y = "Relative Abundance")
      
      ggsave(paste0("analysis_results/figures/interest_genera_", site, ".pdf"), 
             p, width = 10, height = 8)
    }
  }
  
  cat("\n感兴趣类群分析完成！\n")
} else {
  cat("\n警告：未找到指定的感兴趣类群\n")
}

# ============================================================
# 4. 组间差异的统计汇总
# ============================================================

# 汇总所有DESeq2结果
deseq_files <- list.files("analysis_results/tables", 
                         pattern = "DESeq2.*\\.csv", full.names = TRUE)

if(length(deseq_files) > 0) {
  deseq_summary <- data.frame()
  
  for(file in deseq_files) {
    file_name <- basename(file)
    
    # 提取比较信息
    site <- gsub("DESeq2_(.*)_.*\\.csv", "\\1", file_name)
    comparison <- gsub("DESeq2_.*_(.*)\\.csv", "\\1", file_name)
    
    # 读取结果
    res <- read.csv(file)
    
    # 筛选显著差异的类群
    sig_taxa <- res %>%
      filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
      mutate(Site = site, Comparison = comparison)
    
    deseq_summary <- rbind(deseq_summary, sig_taxa)
  }
  
  # 保存汇总结果
  write.csv(deseq_summary, 
           "analysis_results/tables/DESeq2_significant_summary.csv", 
           row.names = FALSE)
  
  # 统计每个比较中显著差异的类群数量
  sig_counts <- deseq_summary %>%
    group_by(Site, Comparison) %>%
    summarise(
      Total_Significant = n(),
      Increased = sum(log2FoldChange > 1),
      Decreased = sum(log2FoldChange < -1),
      .groups = 'drop'
    )
  
  write.csv(sig_counts, 
           "analysis_results/tables/DESeq2_significant_counts.csv", 
           row.names = FALSE)
  
  cat("\nDESeq2结果汇总完成！\n")
  print(sig_counts)
}

cat("\n所有补充分析完成！\n")







# ============================================================
# 跨部位协同变化菌群识别（属水平）
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)

# 创建输出文件夹
if(!dir.exists("analysis_results/cross_site")) {
  dir.create("analysis_results/cross_site", recursive = TRUE)
}

# ===== Step 1：收集所有Site的候选属 =====
cat("\n===== Step 1: 收集候选属 =====\n")

all_candidates <- list()
all_abundance <- list()
all_patterns <- list()

for(site in sites) {
  # 读取候选数据
  rds_file <- paste0("analysis_results/differential/", site, "_Genus_candidates_data.rds")
  sig_file <- paste0("analysis_results/differential/", site, "_Genus_Dunn_significant.csv")
  pattern_file <- paste0("analysis_results/differential/", site, "_pattern_classification.csv")
  
  if(file.exists(rds_file) && file.exists(sig_file)) {
    candidate_data <- readRDS(rds_file)
    sig_results <- read.csv(sig_file)
    
    # 提取显著差异的属
    sig_genera <- unique(sig_results$Feature)
    all_candidates[[site]] <- sig_genera
    
    # 保存丰度数据
    all_abundance[[site]] <- candidate_data$abundance
    
    # 读取变化模式
    if(file.exists(pattern_file)) {
      patterns <- read.csv(pattern_file)
      all_patterns[[site]] <- patterns
    }
    
    cat(paste0("  ", site, ": ", length(sig_genera), " 个显著属\n"))
  }
}

# ===== Step 2：交集分析 =====
cat("\n===== Step 2: 交集分析 =====\n")

# 找出在多个Site都出现的属
all_genera <- unique(unlist(all_candidates))
cat("  总共发现", length(all_genera), "个不重复的候选属\n")

# 统计每个属出现在几个Site
genus_count <- data.frame(
  Genus = all_genera,
  N_sites = sapply(all_genera, function(g) {
    sum(sapply(all_candidates, function(x) g %in% x))
  })
)

# 筛选：在>=2个Site出现的属（共同候选）
common_genera <- genus_count %>%
  filter(N_sites >= 2) %>%
  arrange(desc(N_sites))

cat("  在≥2个Site出现的属:", nrow(common_genera), "个\n")
cat("  在≥3个Site出现的属:", sum(common_genera$N_sites >= 3), "个\n")
cat("  在4个Site都出现的属:", sum(common_genera$N_sites == 4), "个\n")

# 保存
write.csv(genus_count,
         "analysis_results/cross_site/genus_site_count.csv",
         row.names = FALSE)

write.csv(common_genera,
         "analysis_results/cross_site/common_genera.csv",
         row.names = FALSE)

# ===== Step 3：变化方向一致性分析 =====
cat("\n===== Step 3: 变化方向分析 =====\n")

# 对每个共同属，提取在各Site的变化方向
change_direction <- data.frame(
  Genus = common_genera$Genus,
  stringsAsFactors = FALSE
)

for(site in sites) {
  if(!is.null(all_patterns[[site]])) {
    patterns <- all_patterns[[site]]
    
    # 计算变化方向：DB_Treated相对于DB_Control
    direction <- sapply(common_genera$Genus, function(g) {
      if(g %in% patterns$Genus) {
        row <- patterns[patterns$Genus == g, ]
        db <- row$Mean_DB
        tr <- row$Mean_Treated
        
        if(is.na(db) || is.na(tr)) return(NA)
        
        change <- tr - db
        
        # 分类：增加、减少、无变化
        if(abs(change) < 1e-6) {
          return("No_change")
        } else if(change > 0) {
          return("Increase")
        } else {
          return("Decrease")
        }
      } else {
        return(NA)
      }
    })
    
    change_direction[[site]] <- direction
  }
}

# 判断一致性
change_direction$Consistency <- apply(change_direction[, sites], 1, function(row) {
  valid <- row[!is.na(row) & row != "No_change"]
  if(length(valid) == 0) return("No_data")
  if(length(unique(valid)) == 1) {
    return(paste0("Consistent_", unique(valid)))
  } else {
    return("Inconsistent")
  }
})

# 统计
consistency_summary <- table(change_direction$Consistency)
cat("\n  一致性统计:\n")
print(consistency_summary)

# 保存
write.csv(change_direction,
         "analysis_results/cross_site/change_direction_matrix.csv",
         row.names = FALSE)

# 筛选一致性协同变化的属
consistent_genera <- change_direction %>%
  filter(Consistency %in% c("Consistent_Increase", "Consistent_Decrease"))

cat("\n  一致性协同变化属数:", nrow(consistent_genera), "\n")
cat("    - 一致性增加:", sum(consistent_genera$Consistency == "Consistent_Increase"), "\n")
cat("    - 一致性减少:", sum(consistent_genera$Consistency == "Consistent_Decrease"), "\n")

write.csv(consistent_genera,
         "analysis_results/cross_site/consistent_cochanging_genera.csv",
         row.names = FALSE)

# ===== Step 4：计算变化幅度 =====
cat("\n===== Step 4: 计算变化幅度 =====\n")

# 对一致性协同变化的属，计算log2FC
fold_change_matrix <- data.frame(
  Genus = consistent_genera$Genus,
  stringsAsFactors = FALSE
)

for(site in sites) {
  if(!is.null(all_patterns[[site]])) {
    patterns <- all_patterns[[site]]
    
    fc <- sapply(consistent_genera$Genus, function(g) {
      if(g %in% patterns$Genus) {
        row <- patterns[patterns$Genus == g, ]
        db <- row$Mean_DB
        tr <- row$Mean_Treated
        
        if(is.na(db) || is.na(tr) || db == 0) return(NA)
        
        log2((tr + 1e-6) / (db + 1e-6))
      } else {
        return(NA)
      }
    })
    
    fold_change_matrix[[paste0(site, "_log2FC")]] <- fc
  }
}

write.csv(fold_change_matrix,
         "analysis_results/cross_site/fold_change_matrix.csv",
         row.rows = FALSE)

# ===== 可视化 1：UpSet图（交集关系）=====
cat("\n===== 生成可视化 =====\n")

if(requireNamespace("UpSetR", quietly = TRUE)) {
  cat("  绘制UpSet图...\n")
  
  # 准备UpSet数据
  upset_data <- data.frame(
    Genus = all_genera,
    stringsAsFactors = FALSE
  )
  
  for(site in sites) {
    upset_data[[site]] <- as.integer(upset_data$Genus %in% all_candidates[[site]])
  }
  
  pdf("analysis_results/cross_site/upset_plot.pdf", width = 10, height = 6)
  upset(upset_data, 
        sets = sites,
        order.by = "freq",
        main.bar.color = "#e74c3c",
        sets.bar.color = "#3498db",
        text.scale = 1.5)
  dev.off()
}

# ===== 可视化 2：变化方向热图 =====
cat("  绘制变化方向热图...\n")

# 准备热图数据
direction_matrix <- change_direction %>%
  filter(Genus %in% consistent_genera$Genus) %>%
  select(Genus, all_of(sites))

# 转换为数值矩阵
direction_numeric <- data.frame(
  row.names = direction_matrix$Genus,
  stringsAsFactors = FALSE
)

for(site in sites) {
  direction_numeric[[site]] <- sapply(direction_matrix[[site]], function(x) {
    if(is.na(x)) return(0)
    if(x == "Increase") return(1)
    if(x == "Decrease") return(-1)
    return(0)
  })
}

direction_mat <- as.matrix(direction_numeric)

# 添加行注释（一致性类型）
row_anno <- rowAnnotation(
  Type = consistent_genera$Consistency,
  col = list(Type = c(
    "Consistent_Increase" = "#e74c3c",
    "Consistent_Decrease" = "#3498db"
  )),
  annotation_legend_param = list(
    Type = list(title = "Change Type")
  )
)

# 设置颜色
col_fun_direction <- colorRamp2(
  c(-1, 0, 1),
  c("#3498db", "white", "#e74c3c")
)

pdf("analysis_results/cross_site/direction_heatmap.pdf", width = 6, height = max(4, nrow(direction_mat) * 0.2))

ht <- Heatmap(
  direction_mat,
  name = "Direction",
  col = col_fun_direction,
  right_annotation = row_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_title = "Cross-site Consistent Co-changing Genera",
  rect_gp = gpar(col = "grey90", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Change\nDirection",
    at = c(-1, 0, 1),
    labels = c("Decrease", "No change", "Increase")
  )
)

draw(ht)
dev.off()

# ===== 可视化 3：平行坐标图（轨迹）=====
cat("  绘制平行坐标图...\n")

# 准备轨迹数据
trajectory_cross_site <- list()

for(g in consistent_genera$Genus) {
  for(site in sites) {
    if(!is.null(all_patterns[[site]]) && g %in% all_patterns[[site]]$Genus) {
      row <- all_patterns[[site]][all_patterns[[site]]$Genus == g, ]
      
      trajectory_cross_site[[length(trajectory_cross_site) + 1]] <- data.frame(
        Genus = g,
        Site = site,
        WT = row$Mean_WT,
        DB = row$Mean_DB,
        Treated = row$Mean_Treated,
        Type = consistent_genera$Consistency[consistent_genera$Genus == g],
        stringsAsFactors = FALSE
      )
    }
  }
}

trajectory_df <- do.call(rbind, trajectory_cross_site)

# 转为长格式
trajectory_long <- trajectory_df %>%
  pivot_longer(cols = c(WT, DB, Treated),
               names_to = "Group",
               values_to = "Abundance")

trajectory_long$Group <- factor(trajectory_long$Group, levels = c("WT", "DB", "Treated"))

# 分面绘制（按Type和Site）
p_parallel <- ggplot(trajectory_long, 
                     aes(x = Group, y = Abundance, group = Genus, color = Site)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(alpha = 0.7, size = 1.5) +
  facet_wrap(~ Type, scales = "free_y", ncol = 1) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right"
  ) +
  labs(
    title = "Cross-site Trajectories of Consistent Co-changing Genera",
    x = "",
    y = "Mean Relative Abundance",
    color = "Site"
  )

ggsave("analysis_results/cross_site/parallel_coordinates.pdf",
       p_parallel, width = 10, height = 8)

# ===== 可视化 4：变化幅度热图 =====
cat("  绘制变化幅度热图...\n")

fc_mat <- fold_change_matrix %>%
  select(ends_with("_log2FC")) %>%
  as.matrix()

rownames(fc_mat) <- fold_change_matrix$Genus
colnames(fc_mat) <- gsub("_log2FC", "", colnames(fc_mat))

# 移除全NA的行
fc_mat <- fc_mat[rowSums(!is.na(fc_mat)) > 0, ]

col_fun_fc <- colorRamp2(
  c(-2, 0, 2),
  c("#3498db", "white", "#e74c3c")
)

pdf("analysis_results/cross_site/fold_change_heatmap.pdf", width = 6, height = max(4, nrow(fc_mat) * 0.2))

ht_fc <- Heatmap(
  fc_mat,
  name = "log2FC",
  col = col_fun_fc,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_title = "Fold Change (DB_Treated vs DB_Control)",
  na_col = "grey90",
  rect_gp = gpar(col = "grey90", lwd = 0.5),
  heatmap_legend_param = list(
    title = "log2(Fold\nChange)",
    direction = "horizontal"
  )
)

draw(ht_fc, heatmap_legend_side = "bottom")
dev.off()

# ===== 生成总结报告 =====
cat("\n===== 生成总结报告 =====\n")

summary_text <- paste0(
  "========================================\n",
  "跨部位协同变化菌群分析总结\n",
  "========================================\n\n",
  
  "1. 候选属总数:\n",
  paste0("   - Oral: ", length(all_candidates[["Oral"]]), " 个\n"),
  paste0("   - Gut: ", length(all_candidates[["Gut"]]), " 个\n"),
  paste0("   - Fecal: ", length(all_candidates[["Fecal"]]), " 个\n"),
  paste0("   - Lung: ", length(all_candidates[["Lung"]]), " 个\n"),
  paste0("   - 总计（去重）: ", length(all_genera), " 个\n\n"),
  
  "2. 交集分析:\n",
  paste0("   - 在≥2个Site出现: ", nrow(common_genera), " 个\n"),
  paste0("   - 在≥3个Site出现: ", sum(common_genera$N_sites >= 3), " 个\n"),
  paste0("   - 在4个Site都出现: ", sum(common_genera$N_sites == 4), " 个\n\n"),
  
  "3. 一致性协同变化:\n",
  paste0("   - 一致性增加: ", sum(consistent_genera$Consistency == "Consistent_Increase"), " 个\n"),
  paste0("   - 一致性减少: ", sum(consistent_genera$Consistency == "Consistent_Decrease"), " 个\n"),
  paste0("   - 不一致: ", sum(change_direction$Consistency == "Inconsistent"), " 个\n\n"),
  
  "4. 关键发现:\n",
  "   这些一致性协同变化的属在多个取样部位都呈现相同的治疗响应模式，\n",
  "   提示它们可能是Semaglutide治疗的全局性响应标志物。\n\n",
  
  "5. 输出文件:\n",
  "   - genus_site_count.csv: 所有属的Site出现次数\n",
  "   - common_genera.csv: 共同候选属\n",
  "   - change_direction_matrix.csv: 变化方向矩阵\n",
  "   - consistent_cochanging_genera.csv: 一致性协同变化属\n",
  "   - fold_change_matrix.csv: 变化幅度矩阵\n",
  "   - upset_plot.pdf: 交集关系图\n",
  "   - direction_heatmap.pdf: 变化方向热图\n",
  "   - fold_change_heatmap.pdf: 变化幅度热图\n",
  "   - parallel_coordinates.pdf: 跨Site轨迹图\n"
)

writeLines(summary_text, "analysis_results/cross_site/SUMMARY.txt")
cat(summary_text)

cat("\n跨部位协同变化分析完成！\n")
cat("结果保存在: analysis_results/cross_site/\n")