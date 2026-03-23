import pandas as pd
import scipy.stats as stats
import numpy as np
import dansy
from typing import Literal
from collections import Counter
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multitest

class EnrichmentAnalysis():

    def __init__(self, df:pd.DataFrame, ngrams:list, mapper:dict, method:Literal['bh','by'] = 'bh'):
        
        self.data = df
        self.ngrams = ngrams
        self.mapper = mapper
        self.method = method
        self.res = ngram_enrichment(df, ngrams, mapper, method)
        self.ngram2idx = dict(zip(self.res.ngram, self.res.index))

    def get_enrichment(self, ngram):
        return self.res.loc[ngram]
    
    def plot(self,catOI, **kwargs):
        waterfallplot_enrichment(self.res, cat = catOI, **kwargs)
        


def ngram_enrichment(df:pd.DataFrame, ngrams:list, columns:dict, fdr_method:str = 'bh'):
    '''
    Performs enrichment analysis using the Fisher's exact test on a collection of domain architectures for n-grams of interest with both the nominal p-value and a p-adjusted values returned.

    Parameters:
    -----------
    - df: pandas DataFrame
        DataFrame that contains the dataset of interest
    - ngrams: list
        List of n-grams to look for enrichment of
    - columns: dict
        Mapping dictionary with keys 'domain_architecture' and 'category' that tells which columns have the domain architectures and which are the categories to define groups
    - fdr_method: str (Optional)
        The false_discovery_control method from the scipy stats module. If not provided, defaults to Benjamini-Hochberg
    
    Returns:
    --------
    - res: pandas DataFrame
        DataFrame with columns for each category that represents the nominal p-value and an ammended category with _q that represents the adjusted p-value.
    '''

    if 'category' not in columns.keys() or 'domain_architecture' not in columns.keys():
        raise ValueError('The column mapping dict is missing at least one of the required keys.')
    
    # Make sure there are no duplicates in the ngrams provided
    ngrams = sorted(set(ngrams))

    comps = columns['category']
    archs = columns['domain_architecture']
    df_c = df.dropna(subset=archs)
    df_g = df_c.groupby(comps)
    M = len(df_c)

    # Get the counts for each of the categories
    cats_N = Counter(df_c[comps])

    # Run through each of the ngrams of interest and calculate the enrichment using a one-tailed hypergeometric test (aka Fisher's exact)
    raw = []
    for ngram in tqdm(ngrams):
        t = pd.Series(index = [k for k in cats_N.keys()], dtype = float, name = ngram) # Using the object dtype to prevent a FutureWarning from pandas
        n = sum([ngram in arch for arch in df[archs].values])
        t['n'] = n
        for cat, N in cats_N.items():
            temp = df_g.get_group(cat)
            x = sum([ngram in arch for arch in temp[archs].values])
            t[cat] = stats.hypergeom.sf(x-1,M,n,N)
        raw.append(t)

    res = pd.concat(raw, axis = 1).T
    
    for cat in cats_N.keys():
        res[cat+'_q'] = stats.false_discovery_control(res[cat], method = fdr_method)

    res.reset_index(inplace=True)
    res.rename(columns={'index':'ngram'},inplace=True)

    return res

def waterfallplot_enrichment(dataOI:pd.DataFrame,cat, plot_adjusted = True, **kwargs):

    if plot_adjusted:
        col = cat+'_q'
    else:
        col = cat
    
    data = dataOI.copy()
    data[col] = data[col].apply(lambda x: -1*np.log10(x))
    data = data.sort_values(col, ascending=False)

    sns.scatterplot(data, x = range(len(data)), y = col,edgecolor = None, color = 'deepskyblue', **kwargs)
    plt.xlabel('Rank')
    if plot_adjusted:
        plt.ylabel(r'$log_{10}$ p-adjusted')
    else:
        plt.ylabel(r'$log_{10}$ p')

def add_top_ngram_annotations(ax,df:pd.DataFrame, d:dansy.dansy,cat:str,plot_adjusted = True, n:int = 10, yoffset:list = None):
    
    xoffset = 20
    if yoffset is None:
        yoffset = [3]*(n+1)
    
    if plot_adjusted:
        col = cat+'_q'
    else:
        col = cat
    
    data = df.copy()
    data[col] = data[col].apply(lambda x: -1*np.log10(x))
    data = data.sort_values(col, ascending=False,ignore_index= True)

    for i,ngram_info in data[0:n].iterrows():
        ax.annotate(d.return_legible_ngram(ngram_info['ngram']),
                    (i, ngram_info[col]),
                    (xoffset,yoffset[i]),
                    textcoords='offset points',
                    arrowprops={'arrowstyle':"-",
                                'relpos':(0,0.5)},
                    )