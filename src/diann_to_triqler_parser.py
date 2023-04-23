#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 14:22:30 2023

@author: ptruong
"""


filename = "/hdd_14T/data/PXD022992/pxd022992.report.tsv"


import numpy as np
import pandas as pd
import argparse

def diann_to_triqler(filename, qvalue_treshold = 1.00):
    df = pd.read_csv(filename, sep = "\t", 
                     usecols = ["Run", "Precursor.Charge", "Q.Value", "PEP",
                                "Precursor.Normalised", "Stripped.Sequence",
                                "Protein.Ids"])
    df = df[df["Q.Value"] < 1.00]

    df["condition"] = df.Run.map(lambda x: x.split("_")[2]).map({"7951":"metastatic",
                                                                  "SH4":"metastatic",
                                                                  "HTB69":"metastatic",
                                                                  "SK":"primary",
                                                                  "A375":"primary",
                                                                  "G361":"primary"})



    df.rename({"Run":"run", "Precursor.Charge":"charge", "PEP":"PEP",
               "Precursor.Normalised":"intensity", "Stripped.Sequence":"peptide",
               "Protein.Ids":"proteins"}, axis = 1, inplace = True)
    df = df[["run", "condition", "charge", "PEP", "intensity", "peptide", "proteins"]]
    df["searchScore"] = -np.log(df["searchScore"])
    df = df.dropna()
    df.to_csv(output, sep = "\t", index = False)
    return df_triq

def main(input_file, output, qCol = "q"):
    df_triq = diann_to_triqler(input_file)
    #df_triq.intensity = np.log(df_triq.intensity) #if required to log intensity.
    df_triq.to_csv(output, sep = "\t", index = False)

parser = argparse.ArgumentParser(
    description='Converts diann to Triqler input format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name.')

parser.add_argument('--qCol', type=str, default = "q",
                    help='q value column to use.')

parser.add_argument('--output', type=str, default = "msqrob2_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
qCol = args.qCol
output = args.output

if __name__ == "__main__":
    main(input_file, output)

