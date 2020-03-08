#!/usr/bin/python

import pandas as pd

gff1 = pd.read_csv("~/phd/RNASeq/combined_gff_files/version_8/escherichia_5-6_merged.gff", sep="\t")

gff2 = pd.read_csv("~/phd/RNASeq/combined_gff_files/version_8/escherichia_5-4_merged.gff", sep="\t")


gff1 = gff1.loc[ : , gff1.columns[:15]]
gff2 = gff2.loc[ : , gff1.columns[:15]]


gff = pd.concat([gff1, gff2])

gff = gff.reset_index()

for i in range(0, len(gff.index) - 1):
    start_i = gff.loc[i, 'start']
    print(start_i)