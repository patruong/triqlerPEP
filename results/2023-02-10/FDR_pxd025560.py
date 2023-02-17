#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 12:29:37 2023

@author: ptruong
"""

import pandas as pd
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline


df = pd.read_csv("pxd025560.report.tsv", sep = "\t")
df = df.sort_values(by = "PEP")

# calculate PEP from Q.Value
x = df.PEP.reset_index().drop("index", axis = 1)
y = df["Q.Value"].reset_index().drop("index", axis = 1)
f = InterpolatedUnivariateSpline(x, y, k=1)
dfdx = f.derivative()
deriv = pd.DataFrame(dfdx(x), columns = ["PEP_calculated"])
pep_calc = deriv.rename({})

# Calculate q-value from PEP
PEP = df["PEP"].reset_index().drop("index", axis = 1).reset_index()
PEP_sum = PEP.cumsum()
PEP["Q_calculated"] = PEP_sum.PEP / (PEP_sum.index+1)
df.reset_index(inplace = True)
res = pd.concat([df["Protein.Ids"], df["Stripped.Sequence"], PEP, pep_calc, df.iloc[:,df.columns.str.contains("Q.Value")]], axis = 1)
res

df = pd.read_csv("pxd025560.report.tsv", sep = "\t")
df = df.sort_values(by = "PEP")


