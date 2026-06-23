差异基因分析首先对原始的readcount进行标准化(normalization)，主要是对测序深度的校正。然后通过统计学模型进行假设检验概率(p-value)的计算，最后进行多，重假设检验校正(BH)，得到FDR值(错误发现率，padj是其常见形式，以下皆使用padj代表FDR)。

1.deglist/

	*/*_deg.xls 每个比较组合的差异分析列表
		gene_id：基因编号
		sample：各样本标准化后的readcount值
		group：各组readcount均值，若组名的表头是sample.1，表示的是为了和sample区别
		log2FoldChange：处理组与对照组基因表达水平的比值，再以2为底取对数
		pvalue：显著性检验的p值
		padj：多重假设检验校正后的p值
		gene_name：基因名称
		gene_chr：基因所在的染色体名称
		gene_start：基因在染色体的起始位置
		gene_end：基因在染色体的终止位置
		gene_strand：基因所在染色体的正负链信息
		gene_length：基因长度，基因起始到终止所有exon非重叠区域的总和
		gene_biotype：基因类型，如编码蛋白基因，长链非编码基因等
		gene_description：基因功能描述
		gene_tf_family：基因转录因子家族注释

	*/*_deg_all.xls 每个比较组合差异显著性基因列表
	*/*_deg_up.xls 每个比较组合差异显著性上调基因列表
	*/*_deg_down.xls 每个比较组合差异显著性下调基因列表
	*/*_volcano.pdf(svg/png) 每个比较组合差异基因火山图，pdf/svg/png三种格式

    diff_stat.xls：各个比较组合的差异基因数目统计及筛选阈值
		compare：比较组合名称
		all：该比较组合差异基因总数
		up：该比较组合差异基因上调的数目
		down：该比较组合差异基因下调的数目
		threshold：该比较组合筛选差异的阈值
