#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:19:50 2024

@author: andrewmorris
"""

import pandas as pd

#%%

huynhle2021 = pd.read_csv('derived_data/sumstats/PRACTICAL_PCa_Huynh-Le2021.GRCh37.gz', comment = '#', sep = '\t')
karunamuni2021 = pd.read_csv('derived_data/sumstats/PRACTICAL_PCa_Karunamuni2021Africa.GRCh37.gz', comment = '#', sep = '\t')


#%%

kar2021africa = pd.concat([huynhle2021, karunamuni2021])

kar2021africa.to_csv('derived_data/sumstats/PRACTICAL_PCa_Karunamuni2021Africa_all.GRCh37.gz', sep = "\t", index = False, na_rep = "NA")