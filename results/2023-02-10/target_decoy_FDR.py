#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:47:16 2023

@author: ptruong
"""

import pandas as pd
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
# https://www.researchgate.net/figure/Target-and-decoy-scores-are-sorted-and-iteratively-different-thresholds-can-then-be_fig3_283567450
# https://proteomesci.biomedcentral.com/articles/10.1186/s12953-021-00179-7
#(a) Sort target results on score/ p -value/ e -value from best to worst hit. 
#(b) Sort decoy results on score/ p -value/ e -value from best to worst hit. 
#(c) For every target score as threshold, count D and T , and the number of decoys and targets above threshold. 
#(d) Calculate FDR using Eq. 1. This algorithm explains simple FDR using Kall’s method. Other methods differ 
#    in definition of false-positive counts, so should be accordingly calculated. 
#(e) This FDR is also the q -value for the PSM serving as score threshold, and other PSMs with same score. 
#(f) A tabulated list of number of targets and corresponding q -value can be used to create an ROC plot. 3


directory = "../../data/pxd025560/"
file = directory+"pxd025560.report.tsv"

df = pd.read_csv(file, sep = "\t")
df = df.sort_values(by = "PEP")
df.sort_values("CScore", ascending = False, inplace = True)
df["decoy_hit"] = df["Protein.Ids"].str.contains("Random")
df["target_hit"] = ~df["Protein.Ids"].str.contains("Random")
df["decoy_cumsum"] = df.decoy_hit.cumsum()
df["target_cumsum"] = df.target_hit.cumsum()

df["fdrTargetDecoy"] = df["decoy_cumsum"]/df["target_cumsum"]




