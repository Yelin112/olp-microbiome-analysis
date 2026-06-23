# 检查data_HC结构
library(readxl)

# 读取数据
data_HC <- read_excel(
  "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\患者管理\\Clinics_Normal.xlsx"
)

# 显示数据结构
cat("data_HC的结构:\n")
str(data_HC)

cat("\n\n列名:\n")
print(names(data_HC))

cat("\n\n数据类型:\n")
print(sapply(data_HC, class))

cat("\n\n前几行数据:\n")
print(head(data_HC))

# 识别字符串列
string_cols <- sapply(data_HC, function(x) is.character(x) || is.factor(x))
cat("\n\n字符串/因子列:\n")
print(names(data_HC)[string_cols])

# 显示字符串列的唯一值
cat("\n\n字符串列的示例值:\n")
for (col in names(data_HC)[string_cols]) {
  if (is.character(data_HC[[col]])) {
    cat(paste0("\n列 '", col, "' 的唯一值 (前10个):\n"))
    print(head(unique(data_HC[[col]]), 10))
  }
}
