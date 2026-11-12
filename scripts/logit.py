import re
import sys

import numpy as np
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf


matrix_argv = sys.argv[1]    # Dissimilarity matrix - input file
metadata_argv = sys.argv[2]  # Metadata - input file
crude_argv = sys.argv[3]     # Crude Logit model summary - output file 
full_argv = sys.argv[4]      # Full Logit model summary - output file


df = pd.read_csv(matrix_argv,
                 sep="\t", index_col=0)
mdf = pd.read_csv(metadata_argv,
                  sep="\t", index_col=False, comment="#")

mask = np.triu(np.ones(df.shape), k=1).astype(bool)
df_unique = df.where(mask).stack().reset_index().dropna() \
    .reset_index(drop=True)
df_unique.columns = ["sample1", "sample2", "Dissimilarity"]
df_unique.sample2 = df_unique.sample2.astype("int64")

cols = ["sample-id", "host_subject_id"]  # sample ID, patient ID
df_unique = df_unique \
    .merge(mdf[cols], how="left", left_on="sample1", right_on="sample-id") \
    .rename(columns={"host_subject_id": "patient1"})
del df_unique["sample-id"]
df_unique = df_unique \
    .merge(mdf[cols], how="left", left_on="sample2", right_on="sample-id") \
    .rename(columns={"host_subject_id": "patient2"})
del df_unique["sample-id"]
df_unique = df_unique[["patient1", "patient2", "Dissimilarity"]]

func = lambda x: x if x.patient1 == x.patient2 else None
df_pairs = df_unique.apply(func, axis=1) \
    .dropna().reset_index(drop=True) \
    .rename(columns={"patient1": "PatientID"})
del df_pairs["patient2"]
cols = ["host_subject_id", "long_survival"]
full = ["age_cat", "sex", "ana_location", "surgery_year", "asa_cat",
        "differentiation_grade", "stage_tnm", "radical_surgery"]
df_pairs = df_pairs \
    .merge(mdf[cols + full], how="left", left_on="PatientID", right_on="host_subject_id") \
    .drop_duplicates().reset_index(drop=True)
del df_pairs["host_subject_id"]


if "aitchison" in matrix_argv:
    scale = 10
else:
    scale = 0.1


# Crude model

y = df_pairs.long_survival
X = df_pairs.Dissimilarity

X = sm.add_constant(X)
X.Dissimilarity = X.Dissimilarity / scale  # See https://github.com/ctmrbio/crc-survival/blob/main/ipynb/08-Paired%20Beta-Survival.ipynb

model = sm.Logit(y, X)
result = model.fit()

with open(crude_argv, "w", encoding="utf-8") as f:
    print(result.summary2(), file=f)


# Full model

full_ = ("age_cat + sex + ana_location + surgery_year + asa_cat"
         " + differentiation_grade + stage_tnm + radical_surgery")

cols = ["long_survival", "Dissimilarity"]
df_pairs_copy = df_pairs.copy()
df_pairs_copy["long_survival"] = df_pairs_copy["long_survival"].astype(int)
df_pairs_copy["Dissimilarity"] /= scale
model = smf.logit(f"long_survival ~ Dissimilarity + {full_}", data=df_pairs_copy[cols + full])
result = model.fit(disp=False)

with open(full_argv, "w", encoding="utf-8") as f:
    print(result.summary2(), file=f)

