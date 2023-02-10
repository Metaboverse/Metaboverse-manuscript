from __future__ import division
import os
import sys
import pandas as pd
import numpy as np
from numpy import mean, std
from math import sqrt
from scipy.stats import ttest_ind, sem, f
import scipy.stats as st
import xpressplot as xp
import gzip

__path__ = os.getcwd()

def mean_confidence_interval(data, confidence=0.95):
    _interval = st.t.interval(
        alpha=confidence,
        df=len(data)-1,
        loc=np.mean(data),
        scale=st.sem(data)
    )
    _interval = ["{:e}".format(e) for e in _interval]

    return np.mean(data), _interval

def mean_confidence_interval_norm(data, confidence=0.95):
    _interval = st.norm.interval(
        alpha=confidence,
        loc=np.mean(data),
        scale=st.sem(data)
    )
    _interval = ["{:e}".format(e) for e in _interval]
    return np.mean(data), _interval

def cohen_d(x,y):
    return (
        mean(x) - mean(y)) / sqrt((std(x, ddof=1) ** 2 + std(y, ddof=1) ** 2) / 2.0
    )

def confidence_overlap(set1, set2, id, confidence=0.95):

    conf1 = mean_confidence_interval(
        set1.loc[id].values.tolist(),
        confidence=confidence)
    conf2 = mean_confidence_interval(
        set2.loc[id].values.tolist(),
        confidence=confidence)

    conf1 = [float(x) for x in conf1[1]]
    conf2 = [float(x) for x in conf2[1]]

    def getOverlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))

    _overlap = getOverlap(conf1, conf2)
    _diff = max(conf1 + conf2) - min(conf1 + conf2)

    if _overlap > 0:
        print('Confidence intervals (' + str(confidence) + ') OVERLAP')
        print(
            'Overlap: ' + str(_overlap)
            + ' (' + str(round((_overlap / _diff) * 100, 2)) + '%)')
        print('    Interval 1: ' + str(conf1))
        print('    Interval 2: ' + str(conf2))
    else:
        print('Confidence intervals (' + str(confidence) + ') DO NOT overlap')

def g_decompress(
        path,
        file,
        output):

    g_open = open(os.path.join(path, output), "wb")
    with gzip.open(os.path.join(path, file), "rb") as f:
        g_data = f.read()
    g_open.write(g_data)
    g_open.close()

def read_table(
        url,
        sep='\t',
        index_col=0,
        header='infer',
        low_memory=False,
        compression='infer',
        comment='#'):
    """Read tab-delimited table
    """

    data = pd.read_csv(
        url,
        sep=sep,
        index_col=index_col,
        header=header,
        low_memory=low_memory,
        compression=compression,
        comment=comment)

    return data


# RNA d values
rna_url = os.path.abspath(os.path.join(".", "data", "sce_mct1_omics", "sce_mct1_03hr_counts_diffx.tsv"))
rna_de = pd.read_csv(
    rna_url,
    sep='\t',
    index_col=0)

mct1_rnaseq_tpm = read_table(
    url=os.path.join(
        __path__,
        "mct1_analysis",
        "data",
        "mct1_rnaseq_data",
        "sce_mct1_deduped_tpm_threshold25.tsv")
)

gtf = pd.read_csv(
    os.path.join(
        __path__,
        "mct1_analysis",
        "data",
        "analysis_lists",
        "Saccharomyces_cerevisiae.R64-1-1.103.gtf.gz"),
    sep='\t',
    comment='#',
    low_memory=False,
    header=None)
orig_name_label='gene_id "'
orig_name_location=0
new_name_label='gene_name "'
new_name_location=1
gtf_genes = gtf.loc[gtf[2] == 'gene']
gtf_genes['original'] = gtf[8].str.split(';').str[orig_name_location]
gtf_genes['new'] = gtf[8].str.split(';').str[new_name_location]
gtf_genes['original'] = gtf_genes['original'].map(lambda x: x.lstrip(str(orig_name_label)).rstrip('"').lstrip('"').rstrip(' '))
gtf_genes['new'] = gtf_genes['new'].map(lambda x: x.lstrip(str(new_name_label)).rstrip('"').rstrip(' '))
gtf_genes = gtf_genes[['original','new']].copy()

gene_dict = {}
for index, row in gtf_genes.iterrows():
    if row[1] == 'source "sgd':
        gene_dict[row[0]] = row[0]
    else:
        gene_dict[row[0]] = row[1]

rna_de_renamed = rna_de.copy()
rna_de_renamed['new'] = rna_de_renamed.index.to_series().map(gene_dict).fillna(
    rna_de_renamed.index.to_series())
rna_de_renamed = rna_de_renamed.set_index('new')
rna_de_renamed.index.name = None

mct1_rnaseq_tpm_renamed = mct1_rnaseq_tpm.copy()
mct1_rnaseq_tpm_renamed['new'] = mct1_rnaseq_tpm_renamed.index.to_series().map(
        gene_dict
    ).fillna(
        mct1_rnaseq_tpm_renamed.index.to_series())
mct1_rnaseq_tpm_renamed = mct1_rnaseq_tpm_renamed.set_index('new')
mct1_rnaseq_tpm_renamed.index.name = None

"""
14251X4	    WT
14251X6	    mct1del
14251X10	WT
14251X12	mct1del
14251X16	WT
14251X18	mct1del
14251X22	WT
14251X24	mct1del
"""
rna_mct1 = mct1_rnaseq_tpm_renamed[[
    '14251X6',
    '14251X12',
    '14251X18',
    '14251X24']]
rna_wt = mct1_rnaseq_tpm_renamed[[
    '14251X4',
    '14251X10',
    '14251X16',
    '14251X22']]


rna_list = [
    'GLK1',
    'PGI1',
    'PFK1', 'PFK2',
    'FBP1',
    'FBA1', 'TPI1',
    'TDH1', 'TDH2', 'TDH3',
    'PGK1',
    'GPM1',
    'ENO1', 'ENO2',
    'PYK2', 'CDC19',
    'PYC1', 'PYC2',
    'CIT1', 'CIT2', 'CIT3',
    'ACO1', 'ACO2',
    'IDH1', 'IDH2',
    'LPD1', 'KGD1', 'KGD2',
    'LSC1', 'LSC2',
    'SDH1', 'SDH2', 'SDH3', 'SDH4',
    'FUM1',
    'MDH1',
    'GDH1', 'GDH2', 'GLN1',
    'ALT1', 'ALT2',
    'DLD1', 'DLD2', 'DLD3',
    'MAE1',
    'AAT1',
    'CTP1', 'DIC1'
]

for x in rna_list:
    try:
        print(x)
        print('fold change: '+ str(np.log2(
            (sum(rna_mct1.loc[x].values))
            / (sum(rna_wt.loc[x].values))))

        )
        print('Cohen\'s d: '+ str(
            cohen_d(
                rna_mct1.loc[x].values,
                rna_wt.loc[x].values)
            )
        )
        confidence_overlap(rna_mct1, rna_wt, x, confidence=0.95)
        print(
            rna_de_renamed.loc[x]
        )
        print()
        print()
    except:
        print('Skipping ' + x + '...')
        print()
        print()

for x in rna_list:
    try:
        fc = str(round(np.log2(
            (sum(rna_mct1.loc[x].values))
            / (sum(rna_wt.loc[x].values)))
            , 2))
        p = str("{:.2e}".format(rna_de_renamed.loc[x]['FDR']))
        d = str(round(
            cohen_d(
                rna_mct1.loc[x].values,
                rna_wt.loc[x].values)
            , 2))

        print(x, ' & 3 hr & ', fc, ' & ', p, ' & ', d, '\\\\ \hline')
    except:
        pass









# Proteomics d values
raw_proteomics_url = os.path.abspath(os.path.join(".", "mct1_analysis", "data", "mct1_proteomics_data", "proteomics_values.txt"))

proteomics = pd.read_csv(
    raw_proteomics_url,
    sep='\t',
    index_col=0)
proteomics_mct1 = proteomics[[
    'mct1_1',
    'mct1_2',
    'mct1_3']]
proteomics_wt = proteomics[[
    'WT_1',
    'WT_2',
    'WT_3']]

protein_list = [
    'GLK1',
    'PGI1',
    'PFK1', 'PFK2',
    'FBP1',
    'FBA1', 'TPI1',
    'TDH1', 'TDH2', 'TDH3',
    'PGK1',
    'GPM1',
    'ENO1', 'ENO2',
    'PYK2', 'CDC19',
    'PYC1', 'PYC2',
    'CIT1', 'CIT2', 'CIT3',
    'ACO1', 'ACO2',
    'IDH1', 'IDH2',
    'LPD1', 'KGD1', 'KGD2',
    'LSC1', 'LSC2',
    'SDH1', 'SDH2', 'SDH3', 'SDH4',
    'FUM1',
    'MDH1',
    'GDH1', 'GDH2', 'GLN1',
    'ALT1', 'ALT2',
    'DLD1', 'DLD2', 'DLD3',
    'MAE1',
    'AAT1',
    'CTP1', 'DIC1'
]

for x in protein_list:
    print(x)
    print('fold change: '+ str(np.log2(
        (sum(proteomics_mct1.loc[x].values))
        / (sum(proteomics_wt.loc[x].values))))

    )
    print('Cohen\'s d: '+ str(
        cohen_d(
            proteomics_mct1.loc[x].values,
            proteomics_wt.loc[x].values)
        )
    )
    confidence_overlap(proteomics_mct1, proteomics_wt, x, confidence=0.95)
    print(
        ttest_ind(
            proteomics_mct1.loc[x].values,
            proteomics_wt.loc[x].values)
        )
    print()
    print()

for x in protein_list:
    try:
        fc = str(round(np.log2(
            (sum(proteomics_mct1.loc[x].values))
            / (sum(proteomics_wt.loc[x].values)))
            , 2))
        p = str("{:.2e}".format(
            ttest_ind(
                proteomics_mct1.loc[x].values,
                proteomics_wt.loc[x].values)[1]
            ))
        d = str(round(
            cohen_d(
                proteomics_mct1.loc[x].values,
                proteomics_wt.loc[x].values)
            , 2))

        print(x.capitalize(), '(' + x + ') & 12 hr & ', fc, ' & ', p, ' & ', d, '\\\\ \hline')
    except:
        pass










raw_metabolomics_url = os.path.abspath("./mct1_analysis/data/mct1_metabolomics_data/mct1_metabolomics_allValues.txt")
raw_metabolomics_url = os.path.abspath("./mct1_analysis/data/mct1_metabolomics_data/mct1_metabolomics_allValues.txt")
metabolomics_180min = pd.read_csv(
    raw_metabolomics_url,
    sep='\t',
    index_col=0)

metabolomics_mct1_180min = metabolomics_180min[[
    'mct1KO at 3hr minR+Leu 1',
    'mct1KO at 3hr minR+Leu 2',
    'mct1KO at 3hr minR+Leu 3',
    'mct1KO at 3hr minR+Leu 4',
    'mct1KO at 3hr minR+Leu 5',
    'mct1KO at 3hr minR+Leu 6'
]]
metabolomics_wt_180min = metabolomics_180min[[
    'WT at 3hr minR+Leu 1',
    'WT at 3hr minR+Leu 2',
    'WT at 3hr minR+Leu 3',
    'WT at 3hr minR+Leu 4',
    'WT at 3hr minR+Leu 5'
]]

metabolomics_wt_180min.index.tolist()

gc_metabolomics = [
    ' Glucose 6-phosphate',
    ' Fructose-6-phosphate ',
    ' fructose-1,6-diphosphate',
    ' 3-Phosphoglyceric acid ',
    ' Phosphoenolpyruvate ',
    ' Pyruvic acid ',
    ' Glyoxylic acid ',
    ' Citric acid ',
    ' Isocitric acid ',
    ' Succinic acid ',
    ' Fumaric acid ',
    ' D-Malic acid ',
    ' L-Alanine',
    ' L-Lactic acid ',
    ' L-Aspartic acid',
    ' L-Glutamic acid',
    ' 2-Hydroxyglutaric acid',
]

for x in gc_metabolomics:
    print(x)
    print('fold change: '+ str(np.log2(
        (sum(metabolomics_mct1_180min.loc[x].values))
        / (sum(metabolomics_wt_180min.loc[x].values))))

    )
    print('Cohen\'s d: '+ str(
        cohen_d(
            metabolomics_mct1_180min.loc[x].values,
            metabolomics_wt_180min.loc[x].values)
        )
    )
    confidence_overlap(
        metabolomics_mct1_180min,
        metabolomics_wt_180min,
        x,
        confidence=0.95)
    print(
        ttest_ind(
            metabolomics_mct1_180min.loc[x].values,
            metabolomics_wt_180min.loc[x].values)
        )
    print()
    print()




for x in gc_metabolomics:
    try:
        fc = str(round(np.log2(
            (sum(metabolomics_mct1_180min.loc[x].values))
            / (sum(metabolomics_wt_180min.loc[x].values)))
            , 2))
        p = str("{:.2e}".format(
            ttest_ind(
                metabolomics_mct1_180min.loc[x].values,
                metabolomics_wt_180min.loc[x].values)[1]
            ))
        d = str(round(
            cohen_d(
                metabolomics_mct1_180min.loc[x].values,
                metabolomics_wt_180min.loc[x].values)
            , 2))

        print(x[1:-1].capitalize(), '() & 3 hr & ', fc, ' & ', p, ' & ', d, '\\\\ \hline')
    except:
        pass







# CTP1 Metabolomics
raw_metabolomics_ctp1 = os.path.abspath(os.path.join(".", "mct1_analysis", "data", "ctp1_metabolomics_data", "ctp1_metabolomics_allValues.txt"))

metabolomics = pd.read_csv(
    raw_metabolomics_ctp1,
    sep='\t',
    index_col=0)

ctp1_sr_wt_ev = metabolomics[[
    'WT_SR_EV1',
    'WT_SR_EV2',
    'WT_SR_EV3']]

ctp1_sr_wt_a1 = metabolomics[[
    'WT_SR_A1',
    'WT_SR_A2',
    'WT_SR_A3']]

ctp1_sr_mct_ev = metabolomics[[
    'MCT_SR_EV1',
    'MCT_SR_EV2',
    'MCT_SR_EV3']]

ctp1_sr_mct_a1 = metabolomics[[
    'MCT_SR_A1',
    'MCT_SR_A2',
    'MCT_SR_A3']]


lc_metabolomics = [
    'Glucose',
    'F6P',
    'F16BP',
    'Pyruvate',
    'CoA',
    'Citrate',
    'a-KG',
    'Succinate',
    'Fumarate',
    'Malate',
    'Glutamate',
    'Glutamine',
    'Aspartate',
    'Alanine'
]

print('mct1-del vs WT')
for x in lc_metabolomics:
    print(x)
    print('fold change: '+ str(np.log2(
        (sum(ctp1_sr_mct_ev.loc[x].values))
        / (sum(ctp1_sr_wt_ev.loc[x].values))))

    )
    print('Cohen\'s d: '+ str(
        cohen_d(
            ctp1_sr_mct_ev.loc[x].values,
            ctp1_sr_wt_ev.loc[x].values)
        )
    )
    confidence_overlap(
        ctp1_sr_mct_ev,
        ctp1_sr_wt_ev,
        x,
        confidence=0.95)
    print(
        ttest_ind(
            ctp1_sr_mct_ev.loc[x].values,
            ctp1_sr_wt_ev.loc[x].values)
        )
    print()
    print()

for x in lc_metabolomics:
    try:
        fc = str(round(np.log2(
            (sum(ctp1_sr_mct_ev.loc[x].values))
            / (sum(ctp1_sr_wt_ev.loc[x].values)))
            , 2))
        p = str("{:.2e}".format(
            ttest_ind(
                ctp1_sr_mct_ev.loc[x].values,
                ctp1_sr_wt_ev.loc[x].values)[1]
            ))
        d = str(round(
            cohen_d(
                ctp1_sr_mct_ev.loc[x].values,
                ctp1_sr_wt_ev.loc[x].values)
            , 2))

        print(x.capitalize(), '() & 12 hr & ', fc, ' & ', p, ' & ', d, '\\\\ \hline')
    except:
        pass






print('mct1-del+CTP1 vs mct1-del')
for x in lc_metabolomics:
    print(x)
    print('fold change: '+ str(np.log2(
        (sum(ctp1_sr_mct_a1.loc[x].values))
        / (sum(ctp1_sr_mct_ev.loc[x].values))))

    )
    print('Cohen\'s d: '+ str(
        cohen_d(
            ctp1_sr_mct_a1.loc[x].values,
            ctp1_sr_mct_ev.loc[x].values)
        )
    )
    confidence_overlap(
        ctp1_sr_mct_a1,
        ctp1_sr_mct_ev,
        x,
        confidence=0.95)
    print(
        ttest_ind(
            ctp1_sr_mct_a1.loc[x].values,
            ctp1_sr_mct_ev.loc[x].values)
        )
    print()
    print()


for x in lc_metabolomics:
    try:
        fc = str(round(np.log2(
            (sum(ctp1_sr_mct_a1.loc[x].values))
            / (sum(ctp1_sr_mct_ev.loc[x].values)))
            , 2))
        p = str("{:.2e}".format(
            ttest_ind(
                ctp1_sr_mct_a1.loc[x].values,
                ctp1_sr_mct_ev.loc[x].values)[1]
            ))
        d = str(round(
            cohen_d(
                ctp1_sr_mct_a1.loc[x].values,
                ctp1_sr_mct_ev.loc[x].values)
            , 2))

        print(x.capitalize(), '() & 12 hr & ', fc, ' & ', p, ' & ', d, '\\\\ \hline')
    except:
        pass
