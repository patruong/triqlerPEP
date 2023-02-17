#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:54:48 2023

@author: ptruong
"""

import pandas as pd
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import pymc as pm
import pytensor.tensor as pt

import os
os.chdir("/hdd_14T/data/PXD004684/lukas")
os.chdir("/home/ptruong/git/triqler")


def convert_diann_to_triqler(filename):
    df = pd.read_csv(filename, sep = "\t", 
                     usecols = ["Run", "Precursor.Charge", "PEP",
                                "Precursor.Normalised", "Stripped.Sequence",
                                "Protein.Ids"])

    #df = df[df["Q.Value"] < fdr_max]

    df["condition"] = df.Run.map(lambda x: x.split("-")[0][0])



    df.rename({"Run":"run", "Precursor.Charge":"charge", "PEP":"PSM_id_PEP",
               "Precursor.Normalised":"intensity", "Stripped.Sequence":"peptide",
               "Protein.Ids":"proteins"}, axis = 1, inplace = True)
    df = df[["run", "condition", "charge", "PSM_id_PEP", "intensity", "peptide", "proteins"]]
    #df["searchScore"] = -np.log(df["searchScore"])
    df = df.dropna()
    #df.to_csv(output, sep = "\t", index = False)
    return df



triqler = pd.read_csv("triqler_input.csv", sep = "\t")
triqler.columns
triqler[["peptide"]]

data = convert_diann_to_triqler("diann_report.tsv")
delta_offset = 10e-16
data.intensity = np.log2(data.intensity + delta_offset) #we add a small delta to avoid -inf problems
data.intensity.min()
data.intensity.max()


data.columns
data

PSM_PEP = 

# parameters for empirical prior
data.intensity.mean()
data.intensity.std()


coords = {"protein": data["condition"].to_numpy(),
          "condition": data["proteins"].to_numpy()
}


#c_grn 
data.PSM_id_PEP


def sigmoid_matthew(x, mu, sigma): # This is Matthew's logit function, How do we infer mu and sigma from this function and why do we even need these? with pyMC - this is the sigmoid function
    return 0.5 + 0.5*np.tanh((x - mu)/sigma)

def sigmoid(x):
    return 0.5 + 0.5*np.tanh(x)

def logit(p):
    return np.log(p / (1-p))

def expit(p):
    #https://math.stackexchange.com/questions/3816925/how-to-adjust-logit-functions-input-domain
    return np.e**(logit(p)) / (1 - np.e**(logit(p)))


test_data = data.PSM_id_PEP.sample(1000)

with pm.Model(coords=coords) as model:
    
    normaldist = pm.Normal(name = "prior", mu=data.intensity.mean(), sigma=data.intensity.std())
    t_grn = pm.Bernoulli(name = "t_grn, Spurious peptide id", p = test_data) # should it be logit_p or p parameter for pm.Bernoulli?
    m_grn = pm.Bernoulli(name = "m_grn, Missingness", p = expit(test_data))
    
    #m_grn = pm.Bernoulli()
    #idata = pm.sample()

#az.plot_trace(idata)
  
pm.model_to_graphviz(model)


####### Plotting

# Exponential, half-cauchy, Gamma, Weibull or Chi-squared fits missingness
with pm.Model():
    x = pm.Exponential('x', lam=10)
    
    idata = pm.sample()
    
az.plot_trace(idata)


# Import matplotlib, numpy and math
import matplotlib.pyplot as plt
import numpy as np
import math
  
x = np.linspace(-10, 10, 1001)
#z = 1/(1 + np.exp(-x)) # 0.5 off-set to let all intensities at 0 have 100 % missing values
z = 0.5 + 0.5*np.tanh(x)
z = 0.5 + 0.5*np.tanh((x-x.mean())/x.std())
#t = pd.DataFrame([x,z]).T
#t[t[0] == 0]

plt.plot(x, 1-z)
plt.xlabel("x")
plt.ylabel("Sigmoid(X)")
  
plt.show()
