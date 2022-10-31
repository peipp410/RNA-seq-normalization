import pandas as pd
import TranscriptomeNormalization as tm

counts = pd.read_table("./data/small_counts.txt", sep='\t', index_col=0)
sample = tm.CountsMatrix()
sample.load_sample()
# sample.inter_sample(method="TPM", length=[1000, 10000, 2000, 9000, 3000, 8000, 4000, 7000, 5000, 6000])
sample.tmm()
print(sample.get_normalized())

