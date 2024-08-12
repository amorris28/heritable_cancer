#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:23:07 2024

@author: andrewmorris
"""
import pandas as pd

#%% Remove comment rows
sumstats = pd.read_csv('raw_data/sumstats/PGS000067_hmPOS_GRCh37.txt.gz', comment = '#', sep = '\t')
sumstats.head()


#%% Remove ancillary columns
sumstats.columns
sumstats = sumstats.loc[:, 'rsID':'effect_weight']
sumstats.head()

#%% Add case and control numbers manually from paper methods
n_case = 18868
n_control = 12879
sumstats = sumstats.assign(n_cases = n_case,
                           n_controls = n_control)
#%% Export processed sumstats
sumstats.to_csv('derived_data/sumstats/PRACTICAL_PCa_Seibert2018.txt.gz', sep = "\t", index = False, na_rep = "NA")