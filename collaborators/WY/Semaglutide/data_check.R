# ============================================================
# 数据完整性检查脚本
# 运行主分析前使用此脚本检查数据文件和样本名
# ============================================================

# 设置工作路径
base_path <- "E:/打工人/Zeng/2023-3 OLP Microbiome/分析/help_others/WY/Semaglutide"
setwd(base_path)

cat("=====================================\n")
cat("数据完整性检查\n")
cat("=====================================\n\n")

# ============================================================
# 1. 检查必需文件是否存在
# ============================================================

cat("1. 检查必需文件...\n")

files_to_check <- list(
  "Metadata" = "metadata.txt",
  "ASV相对丰度" = "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.relative.xls",
  "门水平相对丰度" = "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.p.relative.xls",
  "纲水平相对丰度" = "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.c.relative.xls",
  "目水平相对丰度" = "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.o.relative.xls",
  "科水平相对丰度" = "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.f.relative.xls",
  "属水平相对丰度" = "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.g.relative.xls",
  "纲水平绝对丰度" = "data/02.ASVanalysis/Taxa_abundance/Evenabs/asv_table.c.absolute.xls",
  "分类信息" = "data/02.ASVanalysis/Seq_taxonomy/seq_taxonomy_qza/taxonomy.tsv"
)

all_files_exist <- TRUE

for (name in names(files_to_check)) {
  file_path <- files_to_check[[name]]
  exists <- file.exists(file_path)

  if (exists) {
    cat(sprintf("  ✓ %s: %s\n", name, file_path))
  } else {
    cat(sprintf("  ✗ %s: %s [文件不存在]\n", name, file_path))
    all_files_exist <- FALSE
  }
}

if (!all_files_exist) {
  cat("\n警告：部分必需文件不存在，请检查路径！\n")
} else {
  cat("\n所有必需文件都存在！\n")
}

# ============================================================
# 2. 检查metadata格式和内容
# ============================================================

cat("\n2. 检查metadata文件...\n")

if (file.exists("metadata.txt")) {
  metadata <- read.table(
    "metadata.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = 1
  )

  cat(sprintf("  样本数: %d\n", nrow(metadata)))
  cat(sprintf("  列数: %d\n", ncol(metadata)))
  cat("  列名:", paste(colnames(metadata), collapse = ", "), "\n")

  # 检查必需的列
  required_cols <- c("Site", "Condition", "Treatment")
  missing_cols <- setdiff(required_cols, colnames(metadata))

  if (length(missing_cols) > 0) {
    cat("  ✗ 缺少必需的列:", paste(missing_cols, collapse = ", "), "\n")
  } else {
    cat("  ✓ 所有必需的列都存在\n")
  }

  # 显示各因素的水平
  cat("\n  因素水平统计:\n")
  cat("    Site:", paste(unique(metadata$Site), collapse = ", "), "\n")
  cat(
    "    Condition:",
    paste(unique(metadata$Condition), collapse = ", "),
    "\n"
  )
  cat(
    "    Treatment:",
    paste(unique(metadata$Treatment), collapse = ", "),
    "\n"
  )

  # 各组样本数
  cat("\n  各组样本数:\n")
  group_counts <- metadata %>%
    group_by(Site, Condition, Treatment) %>%
    summarise(Count = n(), .groups = 'drop')
  print(as.data.frame(group_counts))

  # 检查样本名格式
  cat("\n  样本名示例（前5个）:\n")
  cat("   ", paste(head(rownames(metadata), 5), collapse = ", "), "\n")

  # 检查样本名是否包含下划线
  has_underscore <- any(grepl("_", rownames(metadata)))
  cat(sprintf("  样本名包含下划线: %s\n", ifelse(has_underscore, "是", "否")))
} else {
  cat("  ✗ metadata.txt文件不存在！\n")
}

# ============================================================
# 3. 检查丰度数据文件和样本名匹配
# ============================================================

cat("\n3. 检查丰度数据和样本名匹配...\n")

if (
  file.exists(
    "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.relative.xls"
  )
) {
  asv_data <- read.table(
    "data/02.ASVanalysis/Taxa_abundance/Relative/asv_table.relative.xls",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    nrows = 5
  )

  asv_samples <- colnames(asv_data)
  cat(sprintf("  ASV表中的样本数: %d\n", length(asv_samples)))
  cat("  样本名示例（前5个）:\n")
  cat("   ", paste(head(asv_samples, 5), collapse = ", "), "\n")

  # 检查样本名格式
  has_dot <- any(grepl("\\.", asv_samples))
  cat(sprintf("  样本名包含点号: %s\n", ifelse(has_dot, "是", "否")))

  # 检查metadata和ASV表的样本名匹配
  if (exists("metadata")) {
    # 转换metadata中的样本名（下划线转点号）
    metadata_samples_converted <- gsub("_", ".", rownames(metadata))

    # 找出共同样本
    common_samples <- intersect(metadata_samples_converted, asv_samples)
    only_in_metadata <- setdiff(metadata_samples_converted, asv_samples)
    only_in_asv <- setdiff(asv_samples, metadata_samples_converted)

    cat(sprintf("\n  共同样本数: %d\n", length(common_samples)))
    cat(sprintf("  仅在metadata中的样本数: %d\n", length(only_in_metadata)))
    cat(sprintf("  仅在ASV表中的样本数: %d\n", length(only_in_asv)))

    if (length(only_in_metadata) > 0) {
      cat("  仅在metadata中的样本（前10个）:\n")
      cat("   ", paste(head(only_in_metadata, 10), collapse = ", "), "\n")
    }

    if (length(only_in_asv) > 0) {
      cat("  仅在ASV表中的样本（前10个）:\n")
      cat("   ", paste(head(only_in_asv, 10), collapse = ", "), "\n")
    }

    if (length(common_samples) == nrow(metadata)) {
      cat("  ✓ 所有metadata样本都在ASV表中找到匹配！\n")
    } else {
      cat("  ⚠ 部分metadata样本在ASV表中未找到，可能影响分析\n")
    }
  }
}

# ============================================================
# 4. 检查分类信息文件
# ============================================================

cat("\n4. 检查分类信息文件...\n")

if (
  file.exists("data/02.ASVanalysis/Seq_taxonomy/seq_taxonomy_qza/taxonomy.tsv")
) {
  taxonomy <- read.table(
    "data/02.ASVanalysis/Seq_taxonomy/seq_taxonomy_qza/taxonomy.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    stringsAsFactors = FALSE,
    comment.char = "",
    nrows = 10
  )

  cat(sprintf("  ASV数（仅读取前10行用于检查）\n"))
  cat("  分类示例（前3个）:\n")

  for (i in 1:min(3, nrow(taxonomy))) {
    cat(sprintf(
      "    %s: %s\n",
      rownames(taxonomy)[i],
      substr(taxonomy$taxonomy[i], 1, 80)
    ))
  }

  # 检查分类层级
  sample_tax <- strsplit(taxonomy$taxonomy[1], ";")[[1]]
  cat(sprintf("\n  分类层级数: %d\n", length(sample_tax)))
  cat("  层级示例:", paste(sample_tax, collapse = " > "), "\n")
} else {
  cat("  ✗ taxonomy.tsv文件不存在！\n")
}

# ============================================================
# 5. 检查功能预测文件
# ============================================================

cat("\n5. 检查功能预测文件...\n")

picrust_path <- "E:/semaglutide测序原始数据/扩增子/02.Result_X101SC22114623-Z01-J001_16S/result/05.FunctionPrediction/PICRUSt2/picrust2_out_pipeline"

picrust_files <- c(
  "KO" = file.path(picrust_path, "KO_predicted.tsv.gz"),
  "COG" = file.path(picrust_path, "COG_predicted.tsv.gz"),
  "EC" = file.path(picrust_path, "EC_predicted.tsv.gz")
)

for (name in names(picrust_files)) {
  file_path <- picrust_files[[name]]
  exists <- file.exists(file_path)

  if (exists) {
    cat(sprintf("  ✓ %s预测文件存在\n", name))
  } else {
    cat(sprintf("  ⚠ %s预测文件不存在: %s\n", name, file_path))
  }
}

# ============================================================
# 6. 检查R包安装情况
# ============================================================

cat("\n6. 检查必需的R包...\n")

required_packages <- list(
  "CRAN" = c(
    "ggplot2",
    "vegan",
    "dplyr",
    "tidyr",
    "ape",
    "RColorBrewer",
    "gridExtra",
    "reshape2",
    "scales",
    "ggpubr"
  ),
  "Bioconductor" = c("phyloseq", "DESeq2", "ComplexHeatmap", "microbiome")
)

all_packages_installed <- TRUE

cat("\n  CRAN包:\n")
for (pkg in required_packages$CRAN) {
  installed <- require(pkg, character.only = TRUE, quietly = TRUE)
  if (installed) {
    cat(sprintf("    ✓ %s\n", pkg))
  } else {
    cat(sprintf("    ✗ %s [未安装]\n", pkg))
    all_packages_installed <- FALSE
  }
}

cat("\n  Bioconductor包:\n")
for (pkg in required_packages$Bioconductor) {
  installed <- require(pkg, character.only = TRUE, quietly = TRUE)
  if (installed) {
    cat(sprintf("    ✓ %s\n", pkg))
  } else {
    cat(sprintf("    ✗ %s [未安装]\n", pkg))
    all_packages_installed <- FALSE
  }
}

# ============================================================
# 7. 总结
# ============================================================

cat("\n=====================================\n")
cat("检查总结\n")
cat("=====================================\n")

if (all_files_exist && all_packages_installed) {
  cat("✓ 所有检查通过！可以运行主分析脚本了。\n")
} else {
  cat("⚠ 发现一些问题，请根据上述提示进行修复：\n")

  if (!all_files_exist) {
    cat("  - 部分数据文件不存在，请检查文件路径\n")
  }

  if (!all_packages_installed) {
    cat("  - 部分R包未安装，请运行以下命令安装：\n")
    cat("\n")
    cat("    # 安装CRAN包\n")
    cat("    install.packages(c('ggplot2', 'vegan', 'dplyr', 'tidyr', 'ape',\n")
    cat(
      "                       'RColorBrewer', 'gridExtra', 'reshape2', 'scales', 'ggpubr'))\n"
    )
    cat("\n")
    cat("    # 安装Bioconductor包\n")
    cat("    if (!requireNamespace('BiocManager', quietly = TRUE))\n")
    cat("        install.packages('BiocManager')\n")
    cat(
      "    BiocManager::install(c('phyloseq', 'DESeq2', 'ComplexHeatmap', 'microbiome'))\n"
    )
  }
}

cat("\n检查完成！\n")
