from __future__ import division
import os
import sys
import pandas as pd
from numpy import mean, std
from math import sqrt

def cohen_d(x,y):
    return (
        mean(x) - mean(y)) / sqrt((std(x, ddof=1) ** 2 + std(y, ddof=1) ** 2) / 2.0
    )

# Metabolomics d values
raw_metabolomics_url = os.path.abspath("./data/lung_tumor_pr000305_st000390/measurements_cohensD.txt")
metabolomics = pd.read_csv(
    raw_metabolomics_url,
    sep='\t',
    index_col=0)
metabolomics.index.name = None
metabolomics = metabolomics.T
metabolomics_tumor = metabolomics.loc[metabolomics['group'] == 'tumor'].T
metabolomics_normal = metabolomics.loc[metabolomics['group'] == 'normal'].T

# malate
cohen_d(
    metabolomics_tumor.loc['malic acid'].astype(float).values,
    metabolomics_normal.loc['malic acid'].astype(float).values)

# xanthine
cohen_d(
    metabolomics_tumor.loc['xanthine'].astype(float).values,
    metabolomics_normal.loc['xanthine'].astype(float).values)

# glutamate
cohen_d(
    metabolomics_tumor.loc['glutamic acid'].astype(float).values,
    metabolomics_normal.loc['glutamic acid'].astype(float).values)

# glyceric acid
cohen_d(
    metabolomics_tumor.loc['glyceric acid'].astype(float).values,
    metabolomics_normal.loc['glyceric acid'].astype(float).values)

# 3-Phosphoglyceric acid
cohen_d(
    metabolomics_tumor.loc['3-phosphoglycerate'].astype(float).values,
    metabolomics_normal.loc['3-phosphoglycerate'].astype(float).values)

# citrate
cohen_d(
    metabolomics_tumor.loc['citric acid'].astype(float).values,
    metabolomics_normal.loc['citric acid'].astype(float).values)

# succinate (0.43302502	0.102575479)
cohen_d(
    metabolomics_tumor.loc['succinic acid'].astype(float).values,
    metabolomics_normal.loc['succinic acid'].astype(float).values)

# Hypoxanthine
cohen_d(
    metabolomics_tumor.loc['hypoxanthine'].astype(float).values,
    metabolomics_normal.loc['hypoxanthine'].astype(float).values)

# Spermidine
cohen_d(
    metabolomics_tumor.loc['spermidine'].astype(float).values,
    metabolomics_normal.loc['spermidine'].astype(float).values)

# 5'-methylthioadenosine
cohen_d(
    metabolomics_tumor.loc['5’-deoxy-5’-methylthioadenosine'].astype(float).values,
    metabolomics_normal.loc['5’-deoxy-5’-methylthioadenosine'].astype(float).values)

# Glutamine
cohen_d(
    metabolomics_tumor.loc['glutamine'].astype(float).values,
    metabolomics_normal.loc['glutamine'].astype(float).values)

# alpha-ketoglutarate
cohen_d(
    metabolomics_tumor.loc['alpha ketoglutaric acid'].astype(float).values,
    metabolomics_normal.loc['alpha ketoglutaric acid'].astype(float).values)

# Volcano plot
processed_metabolomics_url = os.path.abspath("./data/lung_tumor_pr000305_st000390/lung_tumor_vs_normal_measurements.txt")
metabolomics_processed = pd.read_csv(
    processed_metabolomics_url,
    sep='\t',
    index_col=0)

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.scatter(
    metabolomics_processed['log2fc'].values,
    -1 * (np.log10(metabolomics_processed['p-adj'].values + 0.000001)),
    c='black',
    edgecolors='grey')
plt.hlines(
    y=1.3,
    xmin=-2,
    xmax=2,
    linestyles='dashed',
    colors='red')
plt.vlines(
    x=0.68,
    ymin=0,
    ymax=6.25,
    linestyles='dashed',
    colors='red')
plt.vlines(
    x=-0.68,
    ymin=0,
    ymax=6.25,
    linestyles='dashed',
    colors='red')
plt.xlabel(
    xlabel='log$_2$(fold change)')
plt.ylabel(
    ylabel='-log$_1$$_0$(p-adj)')
plt.savefig(
    os.path.abspath("./data/lung_tumor_pr000305_st000390/lung_tumor_vs_normal_volcano.png"),
    bbox_inches='tight',
    dpi=1200)
