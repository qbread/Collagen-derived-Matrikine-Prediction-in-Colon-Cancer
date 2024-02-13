#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 13:04:02 2023

@author: unicornshin5103
"""
import pandas as pd
import numpy as np
import os

#COL1A1
COL1A1_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/MMP-14/COL1A1/COL1A1_RNASeq.csv")
COL1A1_normal = COL1A1_TNM.iloc[1::2,:]
COL1A1_tumor = COL1A1_TNM.iloc[0::2,:]
COL1A1_RNA_matched = np.concatenate((COL1A1_tumor, COL1A1_normal), axis = 1)
COL1A1_RNA_matched = pd.DataFrame(COL1A1_RNA_matched)
COL1A1_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/ECM_project/Matrikines/MMP-14/COL1A1/COL1A1_RNA_matched.xlsx")

#COL1A2
COL1A2_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/MMP-14/COL1A2/COL1A2_RNASeq.csv")
COL1A2_normal = COL1A2_TNM.iloc[1::2,:]
COL1A2_tumor = COL1A2_TNM.iloc[0::2,:]
COL1A2_RNA_matched = np.concatenate((COL1A2_tumor, COL1A2_normal), axis = 1)
COL1A2_RNA_matched = pd.DataFrame(COL1A2_RNA_matched)
COL1A2_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/MMP-14/COL1A2/COL1A2_RNA_matched.xlsx")

#COL3A1
COL3A1_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/COL3A1/COL3A1_RNASeq.csv")
COL3A1_normal = COL3A1_TNM.iloc[1::2,:]
COL3A1_tumor = COL3A1_TNM.iloc[0::2,:]
COL3A1_RNA_matched = np.concatenate((COL3A1_tumor, COL3A1_normal), axis = 1)
COL3A1_RNA_matched = pd.DataFrame(COL3A1_RNA_matched)
COL3A1_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/COL3A1/COL3A1_RNA_matched.xlsx")

#COL5A3
COL5A3_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/COL5A3/COL5A3_RNASeq.csv")
COL5A3_normal = COL5A3_TNM.iloc[1::2,:]
COL5A3_tumor = COL5A3_TNM.iloc[0::2,:]
COL5A3_RNA_matched = np.concatenate((COL5A3_tumor, COL5A3_normal), axis = 1)
COL5A3_RNA_matched = pd.DataFrame(COL5A3_RNA_matched)
COL5A3_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/COL5A3/COL5A3_RNA_matched.xlsx")

#COL6A1
COL6A1_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/COL6A1/COL6A1_RNASeq.csv")
COL6A1_normal = COL6A1_TNM.iloc[1::2,:]
COL6A1_tumor = COL6A1_TNM.iloc[0::2,:]
COL6A1_RNA_matched = np.concatenate((COL6A1_tumor, COL6A1_normal), axis = 1)
COL6A1_RNA_matched = pd.DataFrame(COL6A1_RNA_matched)
COL6A1_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/COL6A1/COL6A1_RNA_matched.xlsx")

#COL6A2
COL6A2_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/COL6A2/COL6A2_RNASeq.csv")
COL6A2_normal = COL6A2_TNM.iloc[1::2,:]
COL6A2_tumor = COL6A2_TNM.iloc[0::2,:]
COL6A2_RNA_matched = np.concatenate((COL6A2_tumor, COL6A2_normal), axis = 1)
COL6A2_RNA_matched = pd.DataFrame(COL6A2_RNA_matched)
COL6A2_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/COL6A2/COL6A2_RNA_matched.xlsx")

#COL6A3
COL6A3_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/COL6A3/COL6A3_RNASeq.csv")
COL6A3_normal = COL6A3_TNM.iloc[1::2,:]
COL6A3_tumor = COL6A3_TNM.iloc[0::2,:]
COL6A3_RNA_matched = np.concatenate((COL6A3_tumor, COL6A3_normal), axis = 1)
COL6A3_RNA_matched = pd.DataFrame(COL6A3_RNA_matched)
COL6A3_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/COL6A3/COL6A3_RNA_matched.xlsx")

#COL14A1
COL14A1_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/COL14A1/COL14A1_RNASeq.csv")
COL14A1_normal = COL14A1_TNM.iloc[1::2,:]
COL14A1_tumor = COL14A1_TNM.iloc[0::2,:]
COL14A1_RNA_matched = np.concatenate((COL14A1_tumor, COL14A1_normal), axis = 1)
COL14A1_RNA_matched = pd.DataFrame(COL14A1_RNA_matched)
COL14A1_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/COL14A1/COL14A1_RNA_matched.xlsx")

#MMP-9
MMP9_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/MMP-9/MMP9_RNASeq.csv")
MMP9_normal = MMP9_TNM.iloc[1::2,:]
MMP9_tumor = MMP9_TNM.iloc[0::2,:]
MMP9_RNA_matched = np.concatenate((MMP9_tumor, MMP9_normal), axis = 1)
MMP9_RNA_matched = pd.DataFrame(MMP9_RNA_matched)
MMP9_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/MMP-9/MMP9_RNA_matched.xlsx")

#MMP-14
MMP14_TNM = pd.read_csv("/Users/unicornshin5103/Desktop/Matrikines/MMP-14/MMP14_RNASeq.csv")
MMP14_normal = MMP14_TNM.iloc[1::2,:]
MMP14_tumor = MMP14_TNM.iloc[0::2,:]
MMP14_RNA_matched = np.concatenate((MMP14_tumor, MMP14_normal), axis = 1)
MMP14_RNA_matched = pd.DataFrame(MMP14_RNA_matched)
MMP14_RNA_matched.to_excel("/Users/unicornshin5103/Desktop/Matrikines/MMP-14/MMP14_RNA_matched.xlsx")

