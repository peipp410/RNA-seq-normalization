import numpy as np
import pandas as pd
import math


class CountsMatrix(object):

    def __init__(self, raw=None):
        self.__raw = raw
        self.__normalized = None

    def load_sample(self):
        self.__raw = pd.read_table("./data/small_counts.txt", sep='\t', index_col=0)

    def get_raw(self):
        return self.__raw

    def get_normalized(self):
        if self.__normalized is None:
            print("Warning: Normalized count matrix does not exist!\n")
        return self.__normalized

    def inter_sample(self, method="TPM", paired=False, length=None):
        if len(length) != self.__raw.shape[0]:
            print("Warning: Please check the lengths of all the genes!\n")
        df = self.__raw.copy()
        if method == 'TPM':
            for i in range(len(length)):
                df.iloc[i, :] = df.iloc[i, :]/(length[i]/1000)
            df = df/(df.sum()/1000000)
        if method == 'FPKM' or method == 'RPKM':
            df = df / (df.sum() / 1000000)
            for i in range(len(length)):
                df.iloc[i, :] = df.iloc[i, :]/(length[i]/1000)
        if paired:
            df = df/2
        self.__normalized = df

    def log2_trans(self):
        self.__normalized = np.log(self.__raw+1)/np.log(2)

    # reference: https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe/48249980#48249980
    def quantile(self):
        df = self.__raw.copy()
        dic = {}
        for col in df:
            dic.update({col: sorted(df[col])})
        sorted_df = pd.DataFrame(dic)
        rank = sorted_df.mean(axis=1).tolist()
        for col in df:
            t = np.searchsorted(np.sort(df[col]), df[col])
            df[col] = [rank[i] for i in t]
        self.__normalized = df
        # rank_mean = counts.stack().groupby(counts.rank(method='first').stack().astype(int)).mean()
        # counts.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    def median_of_ratios(self):
        # df = pd.DataFrame(self.__raw.values.T, columns=self.__raw.index, index=self.__raw.columns)
        df = self.__raw.copy()
        for i in range(df.shape[0]):
            df.iloc[i, :] = df.iloc[i, :] / np.exp(np.mean(np.log(df.iloc[i, :])))
        factor = []
        for col in df:
            arr = df[col].values
            arr = arr[~np.isinf(arr)]
            factor.append(np.median(arr))
        sf = pd.DataFrame(data=factor, index=df.columns)
        print("The size factors: \n")
        print(sf)
        for j in range(df.shape[1]):
            df.iloc[:, j] = df.iloc[:, j] / factor[j]
        self.__normalized = df

    def tmm(self, ref_column=None):
        if ref_column is None:
            cpm = self.__raw / self.__raw.sum()
            f75 = cpm.quantile(0.75)
            f75 = abs(f75 - f75.mean())
            ind = f75.index.tolist()
            f75 = f75.tolist()
            refcol = ind[f75.index(min(f75))]
        else:
            refcol = ref_column
        print("The reference column is %s.\n" % refcol)
        lib_size = self.__raw.sum()
        factor = []
        df = self.__raw.copy()
        for col in df:
            obs = df[col].values
            ref = df[refcol].values
            lib_size_obs = lib_size.loc[col, ]
            lib_size_ref = lib_size.loc[refcol, ]
            m = np.log((obs/lib_size_obs) / (ref/lib_size_ref))
            a = (np.log(obs/lib_size_obs) + np.log(ref/lib_size_ref)) / 2
            var = (lib_size_obs-obs)/lib_size_obs/obs + (lib_size_ref-ref)/lib_size_ref/ref
            a_cutoff = -1e10
            fin = ~np.isinf(m) & ~np.isinf(a) & (a > a_cutoff)
            m = m[fin]
            a = a[fin]
            var = var[fin]
            logratio_trim = 0.3  # m
            sum_trim = 0.05  # a
            n = a.shape[0]
            m_low = math.floor(n * logratio_trim) + 1
            m_high = n+1-m_low
            a_low = math.floor(n * sum_trim) + 1
            a_high = n+1-a_low
            m_rank = np.searchsorted(np.sort(m), m) + 1
            a_rank = np.searchsorted(np.sort(a), a) + 1
            keep = (m_rank >= m_low) & (m_rank <= m_high) & (a_rank >= a_low) & (a_rank <= a_high)
            f = np.sum(m[keep]/var[keep]) / np.sum(1/var[keep])
            f = np.exp(f)
            if np.isnan(f):
                f = 1
            factor.append(f)
        factor = np.array(factor)
        factor = factor / np.exp(np.mean(np.log(factor)))
        sf = pd.DataFrame(data=factor, index=df.columns)
        print("The size factors: \n")
        print(sf)
        for col in df:
            df[col] = df[col] / (lib_size.loc[col, ] * sf.loc[col, 0]) * 1e6
        self.__normalized = df


