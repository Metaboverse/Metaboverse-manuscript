"""_summary_
Perform sensitivity analysis of Metaboverse using two test datasets presented in the manuscript

"""
import os
import sys
import pandas as pd
import numpy as np
import random

__path__ = os.getcwd()

DROPOUT_N = 10
DROPOUT_PERCENTS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

### Take random samples from human lung metabolomics dataset

# Create output directory
lung_output = os.path.join(__path__, "lung_random_samples")
if not os.path.exists(lung_output):
    os.makedirs(lung_output)
    
# Read metabolomics and generate drop-outs
lung_metabolomics = pd.read_csv(
    os.path.join(__path__, "..", "data", "lung_tumor_pr000305_st000390", "lung_tumor_vs_normal_measurements.txt"),
    sep="\t",
    index_col=0
)

seed_counter = 0
for i in DROPOUT_PERCENTS:
    for j in range(DROPOUT_N):
        np.random.seed(seed_counter)
        seed_counter += 1
        lung_metabolomics.sample(frac=1-i).to_csv(
            os.path.join(lung_output, "lung_metabolomics_dropout{0}percent_rep{1}.txt".format(int(i*100), j)),
            sep="\t"
        )

### Take random samples from yeast mct1 multi-omics dataset

# Create output directory
mct1_output = os.path.join(__path__, "mct1_random_samples")
if not os.path.exists(mct1_output):
    os.makedirs(mct1_output)

# Read metabolomics and generate drop-outs
mct1_metabolomics = pd.read_csv(
    os.path.join(__path__, "..", "data", "sce_mct1_omics", "mct1_12hr_metabolomics.txt"),
    sep="\t",
    index_col=0
)

seed_counter = 0
for i in DROPOUT_PERCENTS:
    for j in range(DROPOUT_N):
        np.random.seed(seed_counter)
        seed_counter += 1
        mct1_metabolomics.sample(frac=1-i).to_csv(
            os.path.join(mct1_output, "mct1_12hr_metabolomics_dropout{0}percent_rep{1}.txt".format(int(i*100), j)),
            sep="\t"
        )

# Read proteomics and generate drop-outs
mct1_proteomics = pd.read_csv(
    os.path.join(__path__, "..", "data", "sce_mct1_omics", "proteomics_mct1_12hr.txt"),
    sep="\t",
    index_col=0
)

seed_counter = 0
for i in DROPOUT_PERCENTS:
    for j in range(DROPOUT_N):
        np.random.seed(seed_counter)
        seed_counter += 1
        mct1_proteomics.sample(frac=1-i).to_csv(
            os.path.join(mct1_output, "mct1_12hr_proteomics_dropout{0}percent_rep{1}.txt".format(int(i*100), j)),
            sep="\t"
        )

# Read transcriptomics and generate drop-outs
mct1_transcriptomics = pd.read_csv(
    os.path.join(__path__, "..", "data", "sce_mct1_omics", "sce_mct1_12hr_counts_diffx.txt"),
    sep="\t",
    index_col=0
)

seed_counter = 0
for i in DROPOUT_PERCENTS:
    for j in range(DROPOUT_N):
        np.random.seed(seed_counter)
        seed_counter += 1
        mct1_transcriptomics.sample(frac=1-i).to_csv(
            os.path.join(mct1_output, "mct1_12hr_transcriptomics_dropout{0}percent_rep{1}.txt".format(int(i*100), j)),
            sep="\t"
        )