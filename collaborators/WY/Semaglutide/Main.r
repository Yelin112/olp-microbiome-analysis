setwd(
  'E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\help_others\\WY\\Semaglutide'
)


# 1. 定义需要加载的包名向量
my_packages <- c(
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
)

# 2. 批量加载
# lapply 会返回一个列表（每个元素是 TRUE/FALSE），
# 使用 invisible() 可以隐藏这个不必要的返回结果，让控制台更整洁。
invisible(lapply(my_packages, library, character.only = TRUE))


source('data_check.R')

source('microbiome_analysis.R')

source('help_others\WY\Semaglutide\diff_analysis.R')

source('help_others\WY\Semaglutide\RNA_transcriptome_analysis.R')





# 计算菌群紊乱挽救效应 -------------------------------------------------------------
df <- read.csv('help_others\\WY\\Semaglutide\\analysis_results\\dysbiosis\\dysbiosis_scores_Oral.csv')

library(dplyr)
library(ggplot2)
library(tidyr)

# 1. 计算各组均值
group_mean <- df %>%
  group_by(TreatGroup) %>%
  summarise(
    MedianCLV_mean = mean(MedianCLV),
    ShannonJSDScore_mean = mean(ShannonJSDScore),
    CloudScore_mean = mean(CloudScore)
  )

# 2. 计算恢复比例（以平均值为基准）
recovery_group <- data.frame(
  Index = c("MedianCLV","ShannonJSDScore","CloudScore"),
  Recovery_Percent = c(
    (group_mean$MedianCLV_mean[group_mean$TreatGroup=="DB_Treated"] - 
       group_mean$MedianCLV_mean[group_mean$TreatGroup=="DB_Control"]) /
    (group_mean$MedianCLV_mean[group_mean$TreatGroup=="WT_Control"] - 
       group_mean$MedianCLV_mean[group_mean$TreatGroup=="DB_Control"]) * 100,
    
    (group_mean$ShannonJSDScore_mean[group_mean$TreatGroup=="DB_Treated"] - 
       group_mean$ShannonJSDScore_mean[group_mean$TreatGroup=="DB_Control"]) /
    (group_mean$ShannonJSDScore_mean[group_mean$TreatGroup=="WT_Control"] - 
       group_mean$ShannonJSDScore_mean[group_mean$TreatGroup=="DB_Control"]) * 100,
    
    (group_mean$CloudScore_mean[group_mean$TreatGroup=="DB_Treated"] - 
       group_mean$CloudScore_mean[group_mean$TreatGroup=="DB_Control"]) /
    (group_mean$CloudScore_mean[group_mean$TreatGroup=="WT_Control"] - 
       group_mean$CloudScore_mean[group_mean$TreatGroup=="DB_Control"]) * 100
  )
)

recovery_group





# 先计算 DB_Control 和 WT_Control 的组均值
control_mean <- df %>%
  filter(TreatGroup %in% c("DB_Control","WT_Control")) %>%
  group_by(TreatGroup) %>%
  summarise(
    MedianCLV_mean = mean(MedianCLV),
    ShannonJSDScore_mean = mean(ShannonJSDScore),
    CloudScore_mean = mean(CloudScore)
  )

DBControl_mean <- control_mean$MedianCLV_mean[control_mean$TreatGroup=="DB_Control"]
WTControl_mean <- control_mean$MedianCLV_mean[control_mean$TreatGroup=="WT_Control"]

# 对 DB_Treated 样本计算恢复比例
recovery_sample <- df %>%
  filter(TreatGroup=="DB_Treated") %>%
  mutate(
    MedianCLV_Recovery = (MedianCLV - DBControl_mean)/(WTControl_mean - DBControl_mean) * 100
  )

# 同理可以计算其他指标
DBControl_Shannon <- control_mean$ShannonJSDScore_mean[control_mean$TreatGroup=="DB_Control"]
WTControl_Shannon <- control_mean$ShannonJSDScore_mean[control_mean$TreatGroup=="WT_Control"]

recovery_sample <- recovery_sample %>%
  mutate(
    ShannonJSDScore_Recovery = (ShannonJSDScore - DBControl_Shannon)/
                               (WTControl_Shannon - DBControl_Shannon) * 100
  )

DBControl_Cloud <- control_mean$CloudScore_mean[control_mean$TreatGroup=="DB_Control"]
WTControl_Cloud <- control_mean$CloudScore_mean[control_mean$TreatGroup=="WT_Control"]

recovery_sample <- recovery_sample %>%
  mutate(
    CloudScore_Recovery = (CloudScore - DBControl_Cloud)/(WTControl_Cloud - DBControl_Cloud) * 100
  )

# 查看每个样本恢复比例
recovery_sample[, c("SampleID_converted", "MedianCLV_Recovery", "ShannonJSDScore_Recovery", "CloudScore_Recovery")]




library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

#==============================
# 1. 读取数据
#==============================
# df <- read_csv("dysbiosis_scores_Oral.csv") %>%
#   filter(Site == "Oral")

# 设置分组顺序
df <- df %>%
  mutate(
    TreatGroup = factor(
      TreatGroup,
      levels = c("WT_Control", "WT_Treated", "DB_Control", "DB_Treated")
    )
  )

#==============================
# 2. 计算各组统计量（MedianCLV）
#==============================
summary_median <- df %>%
  group_by(TreatGroup) %>%
  summarise(
    n = n(),
    mean = mean(MedianCLV, na.rm = TRUE),
    sd = sd(MedianCLV, na.rm = TRUE),
    se = sd / sqrt(n)
  )

print(summary_median)

#==============================
# 3. 重点三组（WT_Control, DB_Control, DB_Treated）
#==============================
focus_df <- df %>%
  filter(TreatGroup %in% c("WT_Control", "DB_Control", "DB_Treated")) %>%
  mutate(
    TreatGroup = factor(
      TreatGroup,
      levels = c("WT_Control", "DB_Control", "DB_Treated")
    )
  )

focus_summary <- focus_df %>%
  group_by(TreatGroup) %>%
  summarise(
    n = n(),
    mean = mean(MedianCLV, na.rm = TRUE),
    sd = sd(MedianCLV, na.rm = TRUE),
    se = sd / sqrt(n)
  )

print(focus_summary)

#==============================
# 4. 计算组均值恢复比例
#   Recovery = (DB_Treated - DB_Control) / (WT_Control - DB_Control) * 100
#==============================
WT_control_mean <- focus_summary$mean[focus_summary$TreatGroup == "WT_Control"]
DB_control_mean <- focus_summary$mean[focus_summary$TreatGroup == "DB_Control"]
DB_treated_mean <- focus_summary$mean[focus_summary$TreatGroup == "DB_Treated"]

recovery_group <- (DB_treated_mean - DB_control_mean) /
  (WT_control_mean - DB_control_mean) * 100

recovery_group_df <- data.frame(
  Index = "MedianCLV",
  RecoveryPercent = recovery_group
)

print(recovery_group_df)

#==============================
# 5. 计算样本级恢复比例
#==============================
recovery_sample <- df %>%
  filter(TreatGroup == "DB_Treated") %>%
  mutate(
    RecoveryPercent = (MedianCLV - DB_control_mean) /
      (WT_control_mean - DB_control_mean) * 100
  )

print(recovery_sample %>% select(SampleID_converted, MedianCLV, RecoveryPercent))






# 假设 recovery_sample 已经包含 SampleID_converted 和 RecoveryPercent
# 检验恢复比例是否显著大于 0
t.test(recovery_sample$RecoveryPercent, mu = 0, alternative = "greater")

# 检验恢复比例是否显著接近 100
t.test(recovery_sample$RecoveryPercent, mu = 100)

# 或者非参数检验
wilcox.test(recovery_sample$RecoveryPercent, mu = 0, alternative = "greater")
wilcox.test(recovery_sample$RecoveryPercent, mu = 100)








#==============================
# 6. 可视化 1：
#   全部分组原始分布图（箱线图 + 散点 + 均值）
#==============================
p1 <- ggplot(df, aes(x = TreatGroup, y = MedianCLV, fill = TreatGroup)) +
  geom_boxplot(width = 0.55, alpha = 0.45, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.8, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "white") +
  scale_fill_manual(values = c(
    "WT_Control" = "#4DBBD5",
    "WT_Treated" = "#00A087",
    "DB_Control" = "#E64B35",
    "DB_Treated" = "#F39B7F"
  )) +
  labs(
    title = "MedianCLV distribution across groups",
    x = NULL,
    y = "MedianCLV"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#==============================
# 7. 可视化 2：
#   重点趋势图（WT_Control -> DB_Control -> DB_Treated）
#   突出“疾病变化 + 干预挽救”
#==============================
p2 <- ggplot(focus_summary, aes(x = TreatGroup, y = mean, group = 1)) +
  geom_line(linewidth = 0.8, color = "grey40") +
  geom_point(aes(color = TreatGroup), size = 3.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12, linewidth = 0.7) +
  geom_text(aes(label = round(mean, 3)), vjust = -1, size = 3.8) +
  scale_color_manual(values = c(
    "WT_Control" = "#4DBBD5",
    "DB_Control" = "#E64B35",
    "DB_Treated" = "#F39B7F"
  )) +
  labs(
    title = "Trend of MedianCLV: WT_Control vs DB_Control vs DB_Treated",
    x = NULL,
    y = "Mean MedianCLV"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#==============================
# 8. 可视化 3：
#   组均值恢复比例图
#==============================
p3 <- ggplot(recovery_group_df, aes(x = Index, y = RecoveryPercent)) +
  geom_col(fill = "#7E6148", width = 0.55, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_hline(yintercept = 100, linetype = 2, color = "grey50") +
  geom_text(aes(label = paste0(round(RecoveryPercent, 1), "%")), vjust = -0.5, size = 4.2) +
  ylim(min(0, recovery_group_df$RecoveryPercent) - 10,
       max(110, recovery_group_df$RecoveryPercent + 10)) +
  labs(
    title = "Group-level rescue effect of DB_Treated",
    x = NULL,
    y = "Recovery (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#==============================
# 9. 可视化 4：
#   样本级恢复比例图
#==============================
p4 <- ggplot(recovery_sample, aes(x = "DB_Treated", y = RecoveryPercent)) +
  geom_boxplot(width = 0.35, fill = "#F39B7F", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3, color = "#B24745") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_hline(yintercept = 100, linetype = 2, color = "grey50") +
  geom_text(aes(label = SampleID_converted), hjust = -0.15, size = 3.2) +
  labs(
    title = "Sample-level rescue effect in DB_Treated",
    x = NULL,
    y = "Recovery (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#==============================
# 10. 拼图展示
#==============================
final_plot <- (p1 | p2) / (p3 | p4)

final_plot



t.test(MedianCLV ~ TreatGroup,
            data = focus_df %>% filter(TreatGroup %in% c("WT_Control", "DB_Control")))

wilcox.test(MedianCLV ~ TreatGroup,
            data = focus_df %>% filter(TreatGroup %in% c("DB_Control", "DB_Treated")))

wilcox.test(MedianCLV ~ TreatGroup,
            data = focus_df %>% filter(TreatGroup %in% c("WT_Control", "DB_Treated")))