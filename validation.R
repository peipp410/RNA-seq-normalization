library(edgeR)
counts = read.table("./data/small_counts.txt", sep = '\t', header = T)
rownames(counts) = counts$Gene
counts = counts[, -1]
dge = DGEList(counts = counts)
dgelist_norm = calcNormFactors(dge, method = 'RLE')

library(DESeq2)
colData = data.frame(row.names = colnames(counts), 
                     group_list = c(0, 1, 0, 1))
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = colData,
                             design = ~group_list)
dds <- estimateSizeFactors(dds)
dds$sizeFactor

library(preprocessCore)
quantile.normalized = normalize.quantiles(as.matrix(counts))
