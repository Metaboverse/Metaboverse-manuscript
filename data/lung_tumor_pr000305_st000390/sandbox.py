
# Survival curve for TCGA LUAD tumor data
# Input data is from Wang, 2018, Sci Data (data type 2, quantile normalized, but not batch corrected)
try:
    import tarfile
    clin_path = os.path.abspath("./data/lung_tumor_pr000305_st000390")
    clin_url = os.path.join(clin_path, "clinical.project-TCGA-LUAD.2020-10-13.tar.gz")
    tar = tarfile.open(clin_url, "r:gz")
    tar.extractall(path=clin_path)
    tar.close()
except:
    print("File not found")

clinical = pd.read_csv(
     os.path.join(clin_path, "clinical.tsv"),
    header=0,
    sep='\t')
clinical.set_index("case_submitter_id", drop=False, inplace=True)
clinical.index.name = None

clinical_brief = clinical[[
    'case_id',
    'case_submitter_id',
    'project_id',
    'age_at_index',
    'days_to_birth',
    'days_to_death',
    'treatment_or_therapy',
    'gender',
    'ethnicity',
    'race',
    'vital_status',
    'year_of_birth',
    'age_at_diagnosis',
    'days_to_best_overall_response',
    'days_to_diagnosis',
    'days_to_last_follow_up',
    'days_to_last_known_disease_status',
    'days_to_recurrence',
    'tumor_grade',
    'tumor_stage',
    'year_of_diagnosis',
]]
clinical_brief_uniq = clinical_brief.drop_duplicates(subset='case_id')
clinical_brief_uniq_T = clinical_brief_uniq.T

data = pd.read_csv(
    os.path.abspath("./data/lung_tumor_pr000305_st000390/luad-rsem-fpkm-tcga-t.txt.gz"),
    compression='gzip',
    header=0,
    index_col=0,
    sep='\t',
    quotechar='"',
    error_bad_lines=False)
data.index.name = None
data = data.drop('Entrez_Gene_Id', axis=1)
data.columns = data.columns.str.split('-').str[:3].str.join('-')
data_dedup = data.loc[:,~data.columns.duplicated()]

combined = data_dedup.append(clinical_brief_uniq_T, sort=False)
combined_T = combined.T

from lifelines import KaplanMeierFitter

combined_T_c = combined_T.copy()
combined_T_c['days_to_death'] = combined_T_c['days_to_death'].replace("'--", '5050')
combined_T_c["vital_status_binary"] = [1 if x == "Dead" else 0 for x in combined_T_c["vital_status"]]
combined_T_c = combined_T_c.dropna()

combined_T_c['days_to_death'] = combined_T_c['days_to_death'].astype(float)

kmf = KaplanMeierFitter()
T = combined_T_c["days_to_death"]
E = combined_T_c["vital_status_binary"]
kmf.fit(T, event_observed=E)
kmf.survival_function_.plot()
plt.title('Survival function of LUAD patients')

name = "GLYCTK"
ax = plt.subplot(111)
gene = (combined_T_c[name] >= combined_T_c[name].median())
kmf.fit(T[gene], event_observed=E[gene], label=' '.join([name,"high"]), alpha=0.05)
kmf.plot(ax=ax)
kmf.fit(T[~gene], event_observed=E[~gene], label=' '.join([name,"low"]), alpha=0.05)
kmf.plot(ax=ax)
sample_number = str(len(T))
plt.title(' '.join([name, "expression and LUAD for", sample_number, "samples"]))
