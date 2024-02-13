#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 12:48:21 2024

@author: unicornshin5103
"""

# Download Colon CPTAC Proteomics database

import cptac
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from statistics import mean
import math
import seaborn as sns


cptac.list_datasets()
cptac.download("Colon", version = "latest")
colon=cptac.Colon()
colon.list_data()
proteomics = colon.get_proteomics()

# Set each name of row and index in CPTAC proteomics database
samples = proteomics.index
proteins = proteomics.columns


### Match tumor and normal samples
tumor_patients = [k for k in samples if not ".N" in k]
tumor = proteomics.loc[tumor_patients]

normal_patients = [s for s in samples if ".N" in s]
normal = proteomics.loc[normal_patients]

search = '.N'
for i, word in enumerate(normal_patients):
    if search in word: 
        normal_patients[i] = word.strip(search)
        
overlap_sample = list(set(tumor_patients).intersection(normal_patients))

normal_patients = [x + search for x in overlap_sample]

tumor_match = tumor.loc[overlap_sample]
tumor_96 = ['tumor']*96
normal_96 = ['normal']*96
tumor_match.insert(0, 'type', tumor_96) 
normal_match = normal.loc[normal_patients]
normal_match.insert(0, 'type', normal_96)



## 1.Collagen subtypes

# Extract collagen subtypes from CPTAC proteomics database
COL = [i for i in proteins if "COL" in i]
not_collagen = {'COL4A3BP','COLEC10', 'COLEC12', 'COLGALT1', 'PCOLCE'}
collagen_subtypes = [i for i in COL if i not in not_collagen]

# Extract profiles of collagen subtypes each at tumor and normal tissue, excluding any 'NAN' value and storing as excel file
tumor_COL = tumor_match.loc[:, collagen_subtypes]
tumor_COL.insert(0, 'type', tumor_96)
normal_COL = normal_match.loc[:,collagen_subtypes]
normal_COL.insert(0, 'type', normal_96)
entire_COL = pd.concat([tumor_COL, normal_COL], axis = 0)
entire_COL_dropNA = entire_COL.dropna(axis = 1)
#entire_COL_dropNA.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/colon_collagen.xlsx", sheet_name = 'log2 transformed)')

# Calculate relative collagen enrichment in the tumor region to the normal
tumor_COL_dropNA = entire_COL_dropNA[entire_COL_dropNA['type'].str.contains('tumor')]
normal_COL_dropNA = entire_COL_dropNA[entire_COL_dropNA['type'].str.contains('normal')]
tumor_COL_df = tumor_COL_dropNA.drop('type', axis=1)
normal_COL_df = normal_COL_dropNA.drop('type', axis=1)
collagen_DEP_value = tumor_COL_df.values - normal_COL_df.values
collagen_DEP = pd.DataFrame(index = tumor_COL_df.index, columns = tumor_COL_df.columns, data = collagen_DEP_value)
collagen_DEP_raw = 2**collagen_DEP #convert log2-form into raw-value
index1 = collagen_DEP_raw.index.tolist()
index1.append('mean')
collagen_DEP_raw = collagen_DEP_raw.reindex(index1)
collagen_DEP_raw.loc['mean'] = collagen_DEP_raw.mean()
index2 = collagen_DEP_raw.index.tolist()
index2.append('log2_mean')
collagen_DEP_raw = collagen_DEP_raw.reindex(index2)
collagen_DEP_raw.loc['log2_mean'] = np.log2(collagen_DEP_raw.loc['mean'])
#collagen_DEP_raw.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/Collagen_DEP_raw.xlsx", sheet_name = 'collagen_no_NA')

# Select depleted collagen subtypes (log2_mean < -1.0)
collagen_DEP_raw_T = collagen_DEP_raw.T # exchange row and index
collagen_DEP_deplete = collagen_DEP_raw_T[collagen_DEP_raw_T.log2_mean < -1.0]
depleted_collagen = collagen_DEP_deplete.index.to_list()

# Extract the depleted collagen profile of tumor and normal region for later two-tailed paired t-test by Excel
tumor_COL_df_raw = 2**tumor_COL_df
tumor_depleted_collagen = tumor_COL_df_raw.loc[:,depleted_collagen]
tumor_depleted_collagen.insert(0, 'type', tumor_96)
normal_COL_df_raw = 2**normal_COL_df
normal_depleted_collagen = normal_COL_df_raw.loc[:,depleted_collagen]
normal_depleted_collagen.insert(0, 'type', normal_96)
normal_depleted_collagen.index = tumor_depleted_collagen.index
entire_depleted_collagen = pd.concat([tumor_depleted_collagen, normal_depleted_collagen], axis = 1)
#entire_depleted_collagen.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/depleted collagen.xlsx", sheet_name = 'depleted collagen')


## 2.Remodeler subtypes

# Refer to 'Reactome Pathway Browser' to select remodelr subtypes reactive to collagen subtypes significant for both protein and mRNA level (type I, III, and V)
 #<Collagen: Remodeler matches>
 #type I collagen: MMP1, MMP2, MMP8, MMP13, MMP14, MMP15, PRSS2
 #type III collagen: MMP1, MMP8, MMP9, MMP10, MMP13, MMP14, MMP15
 #type V collagen: MMP2, MMP9, MMP10
 
collagen_remodeler = ['MMP1', 'MMP2', 'MMP8', 'MMP9', 'MMP10', 'MMP13', 'MMP14', 'MMP15', 'PRSS2']
remodeler = list(set(collagen_remodeler).intersection(proteins))
tumor_remodelers = tumor_match.loc[:,remodeler]
tumor_remodelers.insert(0, 'type', tumor_96)
normal_remodelers = normal_match.loc[:,remodeler]
normal_remodelers.insert(0, 'type', normal_96)

# Extract profiles of remodeler subtypes each at tumor and normal tissue, excluding any 'NAN' value and storing as excel file
entire_remodelers = pd.concat([tumor_remodelers, normal_remodelers], axis = 0)
entire_remodelers_dropNA = entire_remodelers.dropna(axis = 1)
#entire_remodelers_dropNA.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/colon_collagen_remodelers.xlsx", sheet_name = 'log2 transformed)')

# Calculate relative remodeler enrichment in the tumor region to the normal
tumor_remodelers_dropNA = entire_remodelers_dropNA[entire_remodelers_dropNA['type'].str.contains('tumor')]
normal_remodelers_dropNA = entire_remodelers_dropNA[entire_COL_dropNA['type'].str.contains('normal')]
tumor_remodelers_df = tumor_remodelers_dropNA.drop('type', axis=1)
normal_remodelers_df = normal_remodelers_dropNA.drop('type', axis=1)
remodelers_DEP_value = tumor_remodelers_df.values - normal_remodelers_df.values
remodelers_DEP = pd.DataFrame(index = tumor_remodelers_df.index, columns = tumor_remodelers_df.columns, data = remodelers_DEP_value)
remodelers_DEP_raw = 2**remodelers_DEP #convert log2-form into raw-value
index3 = remodelers_DEP_raw.index.tolist()
index3.append('mean')
remodelers_DEP_raw = remodelers_DEP_raw.reindex(index3)
remodelers_DEP_raw.loc['mean'] = remodelers_DEP_raw.mean()
index4 = remodelers_DEP_raw.index.tolist()
index4.append('log2_mean')
remodelers_DEP_raw = remodelers_DEP_raw.reindex(index4)
remodelers_DEP_raw.loc['log2_mean'] = np.log2(remodelers_DEP_raw.loc['mean'])
#remodelers_DEP_raw.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/Remodelers_DEP_raw.xlsx", sheet_name = 'remodelers_no_NA')

# Select enriched remodeler subtypes (log2_mean > 0.5)
remodelers_DEP_raw_T = remodelers_DEP_raw.T # exchange row and index
remodelers_DEP_enrich = remodelers_DEP_raw_T[remodelers_DEP_raw_T.log2_mean > 0.5]
enriched_remodelers = remodelers_DEP_enrich.index.to_list()

# Extract the enriched profile of tumor and normal region for later two-tailed paired t-test by Excel
tumor_remodelers_df_raw = 2**tumor_remodelers_df
tumor_enriched_remodelers = tumor_remodelers_df_raw.loc[:,enriched_remodelers]
tumor_enriched_remodelers.insert(0, 'type', tumor_96)
normal_remodelers_df_raw = 2**normal_remodelers_df
normal_enriched_remodelers = normal_remodelers_df_raw.loc[:,enriched_remodelers]
normal_enriched_remodelers.insert(0, 'type', normal_96)
normal_enriched_remodelers.index = tumor_enriched_remodelers.index
entire_enriched_remodelers = pd.concat([tumor_enriched_remodelers, normal_enriched_remodelers], axis = 1)
#entire_enriched_remodelers.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/enriched_remodelers.xlsx", sheet_name = 'depleted collagen')


## 3.Integrin subtypes

# Extract profiles of integrin subtypes each at tumor and normal tissue, excluding any 'NAN' value and storing as excel file
ITG = [i for i in proteins if "ITG" in i]
tumor_ITG = tumor_match.loc[:, ITG]
tumor_ITG.insert(0, 'type', tumor_96)
normal_ITG = normal_match.loc[:, ITG]
normal_ITG.insert(0, 'type', normal_96)
entire_ITG = pd.concat([tumor_ITG, normal_ITG], axis = 0)
entire_ITG_dropNA = entire_ITG.dropna(axis = 1)

# Calculate relative integrin enrichment in the tumor region to the normal
tumor_ITG_dropNA = entire_ITG_dropNA[entire_ITG_dropNA['type'].str.contains('tumor')]
normal_ITG_dropNA = entire_ITG_dropNA[entire_ITG_dropNA['type'].str.contains('normal')]
tumor_ITG_df = tumor_ITG_dropNA.drop('type', axis=1)
normal_ITG_df = normal_ITG_dropNA.drop('type', axis=1)
integrin_DEP_value = tumor_ITG_df.values - normal_ITG_df.values
integrin_DEP = pd.DataFrame(index = tumor_ITG_df.index, columns = tumor_ITG_df.columns, data = integrin_DEP_value)
integrin_DEP_raw = 2**integrin_DEP #convert log2-form into raw-value
index5 = integrin_DEP_raw.index.tolist()
index5.append('mean')
integrin_DEP_raw = integrin_DEP_raw.reindex(index5)
integrin_DEP_raw.loc['mean'] = integrin_DEP_raw.mean()
index6 = integrin_DEP_raw.index.tolist()
index6.append('log2_mean')
integrin_DEP_raw = integrin_DEP_raw.reindex(index6)
integrin_DEP_raw.loc['log2_mean'] = np.log2(integrin_DEP_raw.loc['mean'])
#integrin_DEP_raw.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/Integrin_DEP_raw.xlsx", sheet_name = 'integrin_no_NA')


# Select enriched integrin subtypes (log2_mean > 0.5)
integrin_DEP_raw_T = integrin_DEP_raw.T # exchange row and index
integrin_DEP_enrich = integrin_DEP_raw_T[integrin_DEP_raw_T.log2_mean > 0.5]
enriched_integrin = integrin_DEP_enrich.index.to_list()

# Extract the enriched integrin profile of tumor and normal region for later two-tailed paired t-test by Excel
tumor_ITG_df_raw = 2**tumor_ITG_df
tumor_enriched_integrin = tumor_ITG_df_raw.loc[:,enriched_integrin]
tumor_enriched_integrin.insert(0, 'type', tumor_96)
normal_ITG_df_raw = 2**normal_ITG_df
normal_enriched_integrin = normal_ITG_df_raw.loc[:,enriched_integrin]
normal_enriched_integrin.insert(0, 'type', normal_96)
normal_enriched_integrin.index = tumor_enriched_integrin.index
entire_enriched_integrin = pd.concat([tumor_enriched_integrin, normal_enriched_integrin], axis = 1)
#entire_enriched_integrin.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/enriched_integrin.xlsx", sheet_name = 'enriched_integrin')



## 4. Enriched phosphoproteins while their unmodified version depleted in the tumor region

# Extract profiles of phosphoproteins each at tumor and normal tissue, excluding any 'NAN' value and storing as excel file
phosphoproteomics = colon.get_phosphoproteomics()
phosphoproteomics = phosphoproteomics.dropna(axis=1)
phospho_samples = phosphoproteomics.index
phosphoproteins = phosphoproteomics.columns
phospho_tumor = phosphoproteomics.loc[tumor_patients]
phospho_normal = phosphoproteomics.loc[normal_patients]
phospho_tumor_match = phospho_tumor.loc[overlap_sample]
phospho_normal_match = phospho_normal.loc[normal_patients]
tumor_phospho = phospho_tumor_match.loc[:, phosphoproteins]
tumor_phospho.insert(0, 'type', tumor_96)
normal_phospho = phospho_normal_match.loc[:, phosphoproteins]
normal_phospho.insert(0, 'type', normal_96)
entire_phospho = pd.concat([tumor_phospho, normal_phospho], axis = 0)
entire_phospho_dropNA = entire_phospho.dropna(axis = 1)


# 4-1. Calculate relative phosphoprotein enrichment in the tumor region to the normal
tumor_phospho_dropNA = entire_phospho_dropNA[entire_phospho_dropNA['type'].str.contains('tumor')]
normal_phospho_dropNA = entire_phospho_dropNA[entire_phospho_dropNA['type'].str.contains('normal')]
tumor_phospho_df = tumor_phospho_dropNA.drop('type', axis=1)
normal_phospho_df = normal_phospho_dropNA.drop('type', axis=1)
phosphoprotein_DEP_value = tumor_phospho_df.values - normal_phospho_df.values
phosphoprotein_DEP = pd.DataFrame(index = tumor_phospho_df.index, columns = tumor_phospho_df.columns, data = phosphoprotein_DEP_value)
phosphoprotein_DEP_raw = 2**phosphoprotein_DEP #convert log2-form into raw-value
index7 = phosphoprotein_DEP_raw.index.tolist()
index7.append('mean')
phosphoprotein_DEP_raw = phosphoprotein_DEP_raw.reindex(index7)
phosphoprotein_DEP_raw.loc['mean'] = phosphoprotein_DEP_raw.mean()
index8 = phosphoprotein_DEP_raw.index.tolist()
index8.append('log2_mean')
phosphoprotein_DEP_raw = phosphoprotein_DEP_raw.reindex(index8)
phosphoprotein_DEP_raw.loc['log2_mean'] = np.log2(phosphoprotein_DEP_raw.loc['mean'])
#phosphoprotein_DEP_raw.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/Phosphoprotein_DEP_raw.xlsx", sheet_name = 'phosphoprotein_no_NA')


# Select enriched phosphoproteins (log2_mean > 1.0)
phosphoprotein_DEP_raw_T = phosphoprotein_DEP_raw.T # exchange row and index
phosphoprotein_DEP_enrich = phosphoprotein_DEP_raw_T[phosphoprotein_DEP_raw_T.log2_mean > 1.0]
enriched_phosphoprotein = phosphoprotein_DEP_enrich.index.to_list()

# Extract the enriched phosphoprotein profile of tumor and normal region for later two-tailed paired t-test by Excel
tumor_phospho_df_raw = 2**tumor_phospho_df
tumor_enriched_phosphoprotein = tumor_phospho_df_raw.loc[:,enriched_phosphoprotein]
tumor_enriched_phosphoprotein.insert(0, 'type', tumor_96)
normal_phospho_df_raw = 2**normal_phospho_df
normal_enriched_phosphoprotein = normal_phospho_df_raw.loc[:,enriched_phosphoprotein]
normal_enriched_phosphoprotein.insert(0, 'type', normal_96)
normal_enriched_phosphoprotein.index = tumor_enriched_phosphoprotein.index
entire_enriched_phosphoprotein = pd.concat([tumor_enriched_phosphoprotein, normal_enriched_phosphoprotein], axis = 1)
#entire_enriched_phosphoprotein.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/enriched_phosphoprotein.xlsx", sheet_name = 'enriched_phosphoprotein')

# 4-2. Calculate relative enrichment of unmodified version of phosphoproteins enriched in the tumor region 
enriched_phospho = phosphoprotein_DEP_enrich.index.get_level_values(0)
tumor_unphospho = tumor_match.loc[:, enriched_phospho]
tumor_unphospho.insert(0, 'type', tumor_96)
normal_unphospho = normal_match.loc[:, enriched_phospho]
normal_unphospho.insert(0, 'type', normal_96)
entire_unphospho = pd.concat([tumor_unphospho, normal_unphospho], axis = 0)
entire_unphospho_dropNA = entire_unphospho.dropna(axis = 1)
tumor_unphospho_dropNA = entire_unphospho_dropNA[entire_unphospho_dropNA['type'].str.contains('tumor')]
normal_unphospho_dropNA = entire_unphospho_dropNA[entire_unphospho_dropNA['type'].str.contains('normal')]
tumor_unphospho_df = tumor_unphospho_dropNA.drop('type', axis=1)
normal_unphospho_df = normal_unphospho_dropNA.drop('type', axis=1)
unmodified_DEP_value = tumor_unphospho_df.values - normal_unphospho_df.values

unmodified_DEP = pd.DataFrame(index = tumor_unphospho_df.index, columns = tumor_unphospho_df.columns, data = unmodified_DEP_value)
unmodified_DEP_raw = 2**unmodified_DEP #convert log2-form into raw-value
index9 = unmodified_DEP_raw.index.tolist()
index9.append('mean')
unmodified_DEP_raw = unmodified_DEP_raw.reindex(index9)
unmodified_DEP_raw.loc['mean'] = unmodified_DEP_raw.mean()
index10 = unmodified_DEP_raw.index.tolist()
index10.append('log2_mean')
unmodified_DEP_raw = unmodified_DEP_raw.reindex(index10)
unmodified_DEP_raw.loc['log2_mean'] = np.log2(unmodified_DEP_raw.loc['mean'])
#unmodified_DEP_raw.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/Unmodified_phosphoproteins_DEP_raw.xlsx", sheet_name = 'unmodified_no_NA')


# Select depleted unmodified version of the enriched phosphoproteins (log2_mean < 0)
unmodified_DEP_raw_T = unmodified_DEP_raw.T # exchange row and index
unmodified_DEP_deplete = unmodified_DEP_raw_T[unmodified_DEP_raw_T.log2_mean < 0]
depleted_unmodified = unmodified_DEP_deplete.index.to_list()

# Extract the depleted unmodified version of the enriched phosphoprotein of tumor and normal region for later two-tailed paired t-test by Excel
tumor_unphospho_df_raw = 2**tumor_unphospho_df
tumor_depleted_unmodified = tumor_unphospho_df_raw.loc[:,depleted_unmodified]
tumor_depleted_unmodified.insert(0, 'type', tumor_96)
normal_unphospho_df_raw = 2**normal_unphospho_df
normal_depleted_unmodified = normal_unphospho_df_raw.loc[:,depleted_unmodified]
normal_depleted_unmodified.insert(0, 'type', normal_96)
normal_depleted_unmodified.index = tumor_depleted_unmodified.index
entire_depleted_unmodified = pd.concat([tumor_depleted_unmodified, normal_depleted_unmodified], axis = 1)
#entire_depleted_unmodified.to_excel("/Users/unicornshin5103/Desktop/Journal_Submission/Raw_data/depleted_unmodified.xlsx", sheet_name = 'depleted_unmodified')









