#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:38:01 2024

@author: andrewmorris
"""

yaml = [
  "cleansumstats_metafile_kind: minimal",
  "path_sumStats: PRACTICAL_PCa_2018.txt.gz",
  "stats_Model: cox",
  "col_CHR: chr_name",
  "col_POS: chr_position",
  "col_EffectAllele: effect_allele",
  "col_OtherAllele: other_allele",
  "col_BETA: effect_weight",
  "col_SNP: rsID",
  "col_CaseN: n_cases",
  "col_ControlN: n_controls",
  "col_P: P",
  "col_Z: Z"]

with open('derived_data/sumstats/PRACTICAL_PCa_Seibert2018.yaml', 'w') as outfile:
  outfile.write('\n'.join(str(i) for i in yaml))