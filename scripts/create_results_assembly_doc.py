#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
import sys
from os.path import exists
import pandas as pd
mega=sys.argv[1]
cap=sys.argv[2]
out=sys.argv[3]
fileout = open(sys.argv[3],'w')
# Start filling context dic to map values into template
context = {}

if exists(sys.argv[1]):
        with open(sys.argv[1],'r') as asb:
            data=pd.read_csv(sys.argv[1], sep="\t", header=None)
            context["n_cont_mega"]=len(data.index)
            c=data.min()
            context["min_mega"]=c.values[1]
            c=data.max()
            context["max_mega"]=c.values[1]
            c=(data.mean())
            context["mean_mega"]=int(c.values[0])
            c=data.median()
            context["med_mega"]=int(c.values[0])
else:
    print (sys.argv[1]+"do not exit in current dir")
    sys.exit()
#

if exists(sys.argv[2]):
        with open(sys.argv[2],'r') as asb:
            data=pd.read_csv(sys.argv[2], sep="\t", header=None)
            context["n_cont_cap"]=len(data.index)
            c=data.min()
            context["min_cap"]=c.values[1]
            c=data.max()
            context["max_cap"]=c.values[1]
            c=(data.mean())
            context["mean_cap"]=int(c.values[0])
            c=data.median()
            context["med_cap"]=int(c.values[0])
else:
    print (sys.argv[2]+"do not exit in current dir")
    sys.exit()

template = """
# SeekViralReads results for  **_"Aseembly"_**

---
#Assemby megahit/cap3

| variable | value |
| --- | --- |
| Megahit number of contigs  | {n_cont_mega} |
| Megahit min contig length  | {min_mega} bp  |
| Megahit max contig length  | {max_mega} bp  |
| Megahit mean contig length | {mean_mega} bp |
| Megahit median contig length | {med_mega} bp  |

| CAP3 number of contigs  | {n_cont_cap} |
| CAP3 min contig length  | {min_cap} bp  |
| CAP3 max contig length  | {max_cap} bp  |
| CAP3 mean contig length | {mean_cap} bp |
| CAP3 median contig length | {med_cap} bp  |
"""
fileout.write(template.format(**context))
