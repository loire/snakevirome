#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
import sys
from os.path import exists
# import pandas as pd
# import locale
# locale.setlocale(locale.LC_ALL, '')
basename=sys.argv[1]
fileout = open("logs/"+sys.argv[1]+"_report.md",'w')
# Start filling context dic to map values into template
context = {}
context["basename"]=basename
# ("{0:,d}".format(
# First part Cutadapt1
# parse log
if exists("logscutadapt/"+sys.argv[1]+"_1_cut1.log"):
    with open("logscutadapt/"+sys.argv[1]+"_1_cut1.log",'r') as cut1f:
        for line in cut1f:
            if "Command line parameters" in line:
                c =line.split(":")
                context["CutAdapt1CMD"]="cutadapt "+c[1]
                continue
            if "Total reads processed" in line:
                c = line.replace(',','').split(":")
                context["CutA_1_total"]=int(c[1])
                context["CutA_1_total_sum"] =int(c[1])
                continue
            if "Reads with adapters" in line:
                c = line.replace('  ','').split(":")
                context["CutA_1_R1_adap"] = c[1]
                continue
            if "Reads written (passing" in line:
                c=line.split(":")
                context["CutA_1_written"] = c[1]
                c = line.replace(',','').replace('(',':').split(':')
                context["CutA_1_written_sum"] =int(c[2])
                continue
            if "Total basepairs processed" in line:
                c = line.split(":")
                context["CutA_1_LR1"] = c[1]
                line = next(cut1f)
                c = line.split(":")
                context["CutA_1_LR1A"] = c[1]
                break
else:
    print ("Error, file"+sys.argv[1]+"_1_cut1.log do not exit in current dir")
    sys.exit()

if exists("logscutadapt/"+sys.argv[1]+"_2_cut1.log"):
    with open("logscutadapt/"+sys.argv[1]+"_2_cut1.log",'r') as cut1f:
        for line in cut1f:
            if "Command line parameters" in line:
                c =line.split(":")
                context["CutAdapt2CMD"]="cutadapt "+c[1]
                continue
            if "Total reads processed" in line:
                c = line.replace(',','').split(":")
                context["CutA_2_total"]=int(c[1])
                context["CutA_2_total_sum"] =int(c[1])
                context["CutA_1_2_total"]=int(context["CutA_1_total_sum"])+int(context["CutA_2_total_sum"])
                continue
            if "Reads with adapters" in line:
                c = line.split(":")
                context["CutA_1_R2_adap"] = c[1]
                continue
            if "Reads written (passing" in line:
                c=line.split(":")
                context["CutA_2_written"] = c[1]
                c = line.replace(',','').replace('(',':').split(':')
                context["CutA_2_written_sum"] =int(c[2])
                context["CutA_1_2_written"]=int(context["CutA_1_written_sum"])+int(context["CutA_2_written_sum"])
                context["CutA_1_2_written_per_h"]=(context["CutA_1_2_written"]*100)/(context["CutA_1_2_total"])
                continue
            if "Total basepairs processed" in line:
                c = line.split(":")
                context["CutA_1_LR2"] = c[1]
                line = next(cut1f)
                c = line.split(":")
                context["CutA_1_LR2A"] = c[1]
                break
else:
    print ("Error, file"+sys.argv[1]+"_2_cut1.log do not exit in current dir")
    sys.exit()

# Cutadapt2
# parse log
if exists("logscutadapt/"+sys.argv[1]+"_1_cut2.log"):
    with open("logscutadapt/"+sys.argv[1]+"_1_cut2.log",'r') as cut2f:
        for line in cut2f:
            if "Command line parameters" in line:
                c =line.split(":")
                context["CutAdapt1CMD2"]="cutadapt "+c[1]
                continue
            if "Total reads processed" in line:
                c = line.split(":")
                context["CutA_2_1_total"]=c[1]
                continue
            if "Reads with adapter" in line:
                c = line.split(":")
                context["CutA_2_R1_adap"] = c[1]
                continue
            if "Reads written (passing" in line:
                c=line.split(":")
                context["CutA_2_1_written"] = c[1]
                continue
            if "Total basepairs processed" in line:
                c = line.split(":")
                context["CutA_2_LR1"] = c[1]
                line = next(cut2f)
                line = next(cut2f)
                c = line.split(":")
                context["CutA_2_LR1A"] = c[1]
                break
else:
    print ("Error, file "+sys.argv[1]+"_1_cut2.log do not exit in current dir")
    sys.exit()

if exists("logscutadapt/"+sys.argv[1]+"_2_cut2.log"):
    with open("logscutadapt/"+sys.argv[1]+"_2_cut2.log",'r') as cut2f:
        for line in cut2f:
            if "Command line parameters" in line:
                c =line.split(":")
                context["CutAdapt2CMD2"]="cutadapt "+c[1]
                continue
            if "Total reads processed" in line:
                c = line.split(":")
                context["CutA_2_2_total"]=c[1]
                continue
            if "Reads with adapter" in line:
                c = line.split(":")
                context["CutA_2_R2_adap"] = c[1]
                continue
            if "Reads written (passing" in line:
                c=line.split(":")
                context["CutA_2_2_written"] = c[1]
                continue
            if "Total basepairs processed" in line:
                c = line.split(":")
                context["CutA_2_LR2"] = c[1]
                line = next(cut2f)
                line = next(cut2f)
                c = line.split(":")
                context["CutA_2_LR2A"] = c[1]
                break
else:
    print ("Error, file"+sys.argv[1]+"_2_cut2.log do not exit in current dir")
    sys.exit()

#log host mapping widows and pairs
if exists("logsMapHost/"+sys.argv[1]+"_bwa_pairs_on_hist.log"):
    with open("logsMapHost/"+sys.argv[1]+"_bwa_pairs_on_hist.log","r") as bwaf:
        for line in bwaf:
            if "CMD" in line:
                c = line.split(":")
                context["bwaCMDp"]=c[1]
                break
else:
    print ("Error,"+sys.argv[1]+"_bwa_pairs_on_hist.log do not exit in current dir")
    sys.exit()

if exists("logsMapHost/"+sys.argv[1]+"_bwa_widows_on_hist.log"):
    with open("logsMapHost/"+sys.argv[1]+"_bwa_widows_on_hist.log","r") as bwaf:
        for line in bwaf:
            if "CMD" in line:
                c = line.split(":")
                context["bwaCMDw"]=c[1]
                break
else:
    print ("Error,"+sys.argv[1]+"_bwa_widows_on_hist.log do not exit in current dir")
    sys.exit()

if exists("logs_contaminent/Stats_mapping_contaminent_pair_"+sys.argv[1]+".txt"):
    with open("logs_contaminent/Stats_mapping_contaminent_pair_"+sys.argv[1]+".txt") as samfp:
        for line in samfp:
            if "duplicates" in line:
                line = next(samfp)
                c = line.split("+")
                context["bwaMappedp"] = str(int(c[0]))
                d = c[1].split("(")
                e = d[1].split(':')
                context["bwaMappedpcp"] = e[0]
                continue
            if "supplementary" in line :
                c = line.split("+")
                context["bwaInputReadsuplep"] = str(int(c[0]))
                continue
            if "in total" in line:
                c = line.split("+")
                context["bwaInputReadsp"] = int(c[0])
                continue
    context["bwaInputReadsp"] = str(int(context["bwaInputReadsp"])-int(context["bwaInputReadsuplep"]))
    context["bwaOutputp"]=str(int(context["bwaInputReadsp"])-(int(context["bwaMappedp"])))
else:
    print ("logs_contaminent/Stats_mapping_contaminent_pair_"+sys.argv[1]+".txt do not exit in current dir")
    sys.exit()


if exists("logs_contaminent/Stats_mapping_contaminent_widows_"+sys.argv[1]+".txt"):
    with open("logs_contaminent/Stats_mapping_contaminent_widows_"+sys.argv[1]+".txt") as samfw:
        for line in samfw:
            if "duplicates" in line:
                line = next(samfw)
                c = line.split("+")
                context["bwaMappedw"] = str(int(c[0]))
                d = c[1].split("(")
                e = d[1].split(':')
                context["bwaMappedpcw"] = e[0]
                continue
            if "supplementary" in line :
                c = line.split("+")
                context["bwaInputReadsuplew"] = str(int(c[0]))
                continue
            if "in total" in line:
                c = line.split("+")
                context["bwaInputReadsw"] = str(int(c[0]))
                continue
    context["bwaInputReadsw"] = str(int(context["bwaInputReadsw"])-int(context["bwaInputReadsuplew"]))
    context["bwaOutputw"]=str(int(context["bwaInputReadsw"])-(int(context["bwaMappedw"])))
else:
    print ("logs_contaminent/Stats_mapping_contaminent_widows_"+sys.argv[1]+".txt do not exit in current dir")
    sys.exit()


if exists("logsFLASH/"+sys.argv[1]+"_flash.log"):
    with open("logsFLASH/"+sys.argv[1]+"_flash.log",'r') as flf :
        for line in flf:
            if "Total pairs" in line:
                c = line.split(":")
                context["flashInput"] = c[1]
                line = next(flf)
                c = line.split(":")
                context["flashCombined"] = c[1]
                line = next(flf)
                c = line.split(":")
                context["flashUncombined"] = c[1]
                line = next(flf)
                c = line.split(":")
                context["flashCombinedpc"] = c[1]
                context["flashOutput"] = str(int(context["flashCombined"])+ 2*int(context["flashUncombined"]))
                break

else:
    print (sys.argv[1]+"_flash.log do not exit in current dir")
    sys.exit()


if exists('logs/'+sys.argv[1]+"_Aedesdataset_stats_mapping_assembly.txt"):
    with open("logs/"+sys.argv[1]+"_Aedesdataset_stats_mapping_assembly.txt",'r') as samf:
        for line in samf:
            if "QC-passed" in line:
                c = line.split("+")
                context["num_seq"] = c[0]
                line = next(samf)
                line = next(samf)
                c = line.split("+")
                context["supl_assem"] = c[0]
                line = next(samf)
                line = next(samf)
                c = line.replace(' ',',').replace(',','(').split('(')
                c1=str(c[0])+' ('+str(c[5])+')'
                context["mapped"] = c1
                break
else:
    print (sys.argv[1]+"_Aedesdataset_stats_mapping_assembly.txt do not exit in current dir")
    sys.exit()
context["num_seq"]=int(context["num_seq"])-int(context["supl_assem"])

if exists("logs/"+sys.argv[1]+"_Aedesdataset_blast_ntvir.txt"):
    with open("logs/"+sys.argv[1]+"_Aedesdataset_blast_ntvir.txt",'r') as f:
        for line in f:
            if "Unmapped" in line:
                c1 = int(line.split(":")[1])
                line = next(f)
                c2 = int(line.split(":")[1])
                c=c1+c2
                context["to_blast"]=c
                continue
            if "total" in line:
                c=line.replace(' ','').split("total\n")
                context["blast_hits"]=int(c[0])
                c=round((context["blast_hits"]*100)/(context["to_blast"]),2)
                context["blast_hits_h"]=c
                line = next(f)
                c = line.split("\n")
                context["to_blast_contigs"]=int(c[0])
                line = next(f)
                c = line.split(" ")
                context["blast_hits_contigs"]=int(c[0])
                c=round((context["blast_hits_contigs"]*100)/(context["to_blast_contigs"]),2)
                context["blast_hits_h_contigs"]=c
                break
else:
    print (sys.argv[1]+ "_Aedesdataset_blast_ntvir.txt do not exit in current dir")
    sys.exit()




template = """
# SeekViralReads results for sample **_"{basename}"_**

---
# Adapters removal by cutadapt

## command line:
```bash
{CutAdapt1CMD}
{CutAdapt2CMD}
```

| variable | value |
| --- | --- |
| Total read pairs processsed R1| {CutA_1_total} |
| Total read pairs processsed R2| {CutA_2_total} |
| R1 with adapters | {CutA_1_R1_adap}|
| R2 with adapters | {CutA_1_R2_adap} |
| R1 total length before | {CutA_1_LR1} |
| R1 total length after | {CutA_1_LR1A} |
| R2 total length before | {CutA_1_LR2} |
| R2 total length after | {CutA_1_LR2A} |
| Reads 1 written | {CutA_1_written} |
| Reads 2 written | {CutA_2_written} |
| Reads pairs written | {CutA_1_2_written} ({CutA_1_2_written_per_h}%) |


## Quality trimming by cutadapt

### command line:
```bash
{CutAdapt1CMD2}
{CutAdapt2CMD2}
```

| variable | value |
| --- | --- |
| Total reads pairs processsed | {CutA_2_1_total} |
| Total reads pairs processsed2 | {CutA_2_2_total} |
| R1 with adapters | {CutA_2_R1_adap} |
| R2 with adapters | {CutA_2_R2_adap} |
| R1 total length before2 | {CutA_2_LR1} |
| R1 total length after2 | {CutA_2_LR1A} |
| R2 total length before2 | {CutA_2_LR2} |
| R2 total length after2 | {CutA_2_LR2A} |
| Reads 1 written2 | {CutA_2_1_written} |
| Reads 2 written2 | {CutA_2_2_written} |


## Mapping pairs on contaminent sequences

### command line:
```bash
{bwaCMDp}
```

| variable | value |
| --- | --- |
| Reads input p  | {bwaInputReadsp} |
| Supplementary p | {bwaInputReadsuplep} |
| Properly Mapped p | {bwaMappedp} |
| Properly Mapped %  p| {bwaMappedpcp} |
| Reads pair written  p| {bwaOutputp} |


## Mapping widows on contaminent sequences

### command line:
```bash
{bwaCMDw}
```

| variable | value |
| --- | --- |
| Reads input w  | {bwaInputReadsw} |
| Supplementary w | {bwaInputReadsuplew}|
| Properly Mapped w| {bwaMappedw} |
| Properly Mapped % w| {bwaMappedpcw} |
| Reads pair written w| {bwaOutputw} |


## Merging overlapping reads with flash

| variable | value |
| --- | --- |
| Reads input f | {flashInput}  |
| Combined pairs | {flashCombined} |
| Percent combined | {flashCombinedpc} |
| Uncombined pairs | {flashUncombined} |
| Number of output sequences | {flashOutput} |


##Mapping sequences on assembly

| variable | value |
| --- | --- |
| Reads input a | {num_seq} |
| Supplementary a | {supl_assem} |
| Reads mapped | {mapped} |


##Blast on nt_vir

| variable | value |
| --- | --- |
| Reads input b | {to_blast} |
| Reads with blast hit | {blast_hits} ({blast_hits_h}%) |
| Contigs input | {to_blast_contigs} |
| Contigs with blast hit | {blast_hits_contigs} ({blast_hits_h_contigs}%) |
"""

fileout.write(template.format(**context))


# Assemby megahit/cap3
#
# | variable | value |
# | --- | --- |
# | Megahit number of contigs  | {n_cont_mega} |
# | Megahit min contig length  | {min_mega} bp  |
# | Megahit max contig length  | {max_mega} bp  |
# | Megahit mean contig length | {mean_mega} bp |
# | Megahit median contig length | {med_mega} bp  |
#
# | CAP3 number of contigs  | {n_cont_cap} |
# | CAP3 min contig length  | {min_cap} bp  |
# | CAP3 max contig length  | {max_cap} bp  |
# | CAP3 mean contig length | {mean_cap} bp |
# | CAP3 median contig length | {med_cap} bp  |
