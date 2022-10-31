# RNA-seq normalization algorithms comparison

Multiple normalization algorithms have been developed to alleviate different sources of biased while analyzing RNA-seq data. Due to some issues left over from the past, most of the packages are written in R. Therefore, it's of vital necessity to implement them in Python. Here, I will implement these algorithms in Python to help me get profound understandings of these algorithms.

## Load the packages

```python
import pandas as pd
import TranscriptomeNormalization as tm
```

## Initialization

You can initialize with the sample data

```python
sample = tm.CountsMatrix()
sample.load_sample()
```

Or you can initialize with the count matrix loaded by your self.

```python
counts = pd.read_table("./data/small_counts.txt", sep='\t', index_col=0)
sample = tm.CountsMatrix(counts)
```

Each column of the matrix represents a sample and each row represents a gene. The gene name should be set as the index of this data frame.

## FPKM/RPKM/TPM

These algorithms applies to inter-sample normalization. You can refer to [StatQuest: RPKM, FPKM and TPM, Clearly Explained!!!](https://www.youtube.com/watch?v=TTUrtCY2k-w&list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp&index=6) for the details.

```python
sample.inter_sample(method="TPM", paired=False, length=None)
```

1. method: FPKM or RPKM or TPM, default: TPM
2. paired: whether this is a paired-end library, default: False. When True, we have FPKM = RPKM/2.
3. length: A list containing the lengths of all the genes corresponding to the index.

## Quantile Normalization

The same as the `normalize.quantiles` function in R package `preprocessCore`. However, that's not the `upperquartile`  method in the `edgeR` function `calcNormFactors`.

```Python
sample.quantile()
```

## Median of Ratios Method

The same as the `estimateSizeFactors` function in R package `DESeq2`. However, that's not the `RLE`  method in the `calcNormFactors` function of `edgeR`. You can refer to [StatQuest: DESeq2, part 1, Library Normalization](https://www.youtube.com/watch?v=UFB993xufUU&list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp&index=11) for the details.

```python
sample.median_of_ratios()
```

## Trimmed Mean of M-values

The same as the `TMM` method in the `calcNormFactors` function of `edgeR`. You can refer to [StatQuest: edgeR, part 1, Library Normalization](https://www.youtube.com/watch?v=Wdt6jdi-NQo&list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp&index=12) for the details.

```python
sample.tmm(ref_column=None)
```

ref_column: The reference column used for normalization. You can specific one column, or just let the algorithm to decide which one to be the reference column.

## Report the normalized count matrix

```python
sample.get_normalized()
```

