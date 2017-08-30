#!/usr/local/bioinfo/python/2.7.9_build2/bin/python
# -*- coding: utf-8 -*-
import sys
import os
from os.path import exists

path=sys.argv[1]
context = {}
fileout = open("logs/Aedes_report_merge.md",'w')
for file in os.listdir(path):
    if file.endswith("report.md"):
        with open(sys.argv[1]+file,'r') as log_report :
            for line in log_report:
                if "Total read pairs processsed R1" in line:
                    c =line.replace(',','').split("|")
                    if "CutA_1_total" not in context:
                        context["CutA_1_total"]=int(c[2])
                    else:
                        context["CutA_1_total"]=int(context["CutA_1_total"])+int(c[2])

                if "Total read pairs processsed R2" in line:
                    c =line.replace(',','').split("|")
                    if "CutA_2_total" not in context:
                        context["CutA_2_total"]=int(c[2])
                    else:
                        context["CutA_2_total"]=int(context["CutA_2_total"])+int(c[2])


                if "R1 with adapters" in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_1_R1_adap" not in context:
                        context["CutA_1_R1_adap"]=int(c[2])
                    else:
                        context["CutA_1_R1_adap"]=int(context["CutA_1_R1_adap"])+int(c[2])


                if "R2 with adapters" in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_1_R2_adap" not in context:
                        context["CutA_1_R2_adap"]=int(c[2])
                    else:
                        context["CutA_1_R2_adap"]=int(context["CutA_1_R2_adap"])+int(c[2])


                if "R1 total length before" in line:
                    c =line.replace(' bp','').replace(',','').split("|")
                    if "CutA_1_LR1" not in context:
                        context["CutA_1_LR1"]=int(c[2])
                    else:
                        context["CutA_1_LR1"]=int(context["CutA_1_LR1"])+int(c[2])


                if "R1 total length after " in line:
                    c =line.replace(' bp','|').replace(',','').split("|")
                    if "CutA_1_LR1A" not in context:
                        context["CutA_1_LR1A"]=int(c[2])
                    else:
                        context["CutA_1_LR1A"]=int(context["CutA_1_LR1A"])+int(c[2])


                if "R2 total length before" in line:
                    c =line.replace(' bp','').replace(',','').split("|")
                    if "CutA_1_LR2" not in context:
                        context["CutA_1_LR2"]=int(c[2])
                    else:
                        context["CutA_1_LR2"]=int(context["CutA_1_LR1"])+int(c[2])


                if "R2 total length after" in line:
                    c =line.replace(' bp','|').replace(',','').split("|")
                    if "CutA_1_LR2A" not in context:
                        context["CutA_1_LR2A"]=int(c[2])
                    else:
                        context["CutA_1_LR2A"]=int(context["CutA_1_LR1"])+int(c[2])


                if "Reads 1 written " in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_1_written" not in context:
                        context["CutA_1_written"]=int(c[2])
                    else:
                        context["CutA_1_written"]=int(context["CutA_1_written"])+int(c[2])

                if "Reads 2 written " in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_2_written" not in context:
                        context["CutA_2_written"]=int(c[2])
                    else:
                        context["CutA_2_written"]=int(context["CutA_2_written"])+int(c[2])

                if "Reads pairs written" in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_1_2_written" not in context:
                        context["CutA_1_2_written"]=int(c[2])
                    else:
                        context["CutA_1_2_written"]=int(context["CutA_1_2_written"])+int(c[2])

                if "Total reads pairs processsed " in line:
                    c =line.replace(',','').split("|")
                    if "CutA_2_1_total" not in context:
                        context["CutA_2_1_total"]=int(c[2])
                    else:
                        context["CutA_2_1_total"]=int(context["CutA_2_1_total"])+int(c[2])


                if "Total reads pairs processsed2" in line:
                    c =line.replace(',','').split("|")
                    if "CutA_2_2_total" not in context:
                        context["CutA_2_2_total"]=int(c[2])
                    else:
                        context["CutA_2_2_total"]=int(context["CutA_2_2_total"])+int(c[2])


                if "R1 total length before2" in line:
                    c =line.replace(' bp','').replace(',','').split("|")
                    if "CutA_2_LR1" not in context:
                        context["CutA_2_LR1"]=int(c[2])
                    else:
                        context["CutA_2_LR1"]=int(context["CutA_2_LR1"])+int(c[2])



                if "R1 total length after2" in line:
                    c =line.replace(' bp','|').replace(',','').split("|")
                    if "CutA_2_LR1A" not in context:
                        context["CutA_2_LR1A"]=int(c[2])
                    else:
                        context["CutA_2_LR1A"]=int(context["CutA_2_LR1A"])+int(c[2])

                if "R2 total length before2" in line:
                    c =line.replace(' bp','').replace(',','').split("|")
                    if "CutA_2_LR2" not in context:
                        context["CutA_2_LR2"]=int(c[2])
                    else:
                        context["CutA_2_LR2"]=int(context["CutA_2_LR2"])+int(c[2])

                if "R2 total length after" in line:
                    c =line.replace(' bp','|').replace(',','').split("|")
                    if "CutA_2_LR2A" not in context:
                        context["CutA_2_LR2A"]=int(c[2])
                    else:
                        context["CutA_2_LR2A"]=int(context["CutA_2_LR2A"])+int(c[2])

                if "Reads 1 written2" in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_2_1_written" not in context:
                        context["CutA_2_1_written"]=int(c[2])
                    else:
                        context["CutA_2_1_written"]=int(context["CutA_2_1_written"])+int(c[2])

                if "Reads 2 written2" in line:
                    c =line.replace(' (','|').replace(',','').split("|")
                    if "CutA_2_2_written" not in context:
                        context["CutA_2_2_written"]=int(c[2])
                    else:
                        context["CutA_2_2_written"]=int(context["CutA_2_2_written"])+int(c[2])

                if "Reads input p" in line:
                    c =line.split("|")
                    if "bwaInputReadsp" not in context:
                        context["bwaInputReadsp"]=int(c[2])
                    else:
                        context["bwaInputReadsp"]=int(context["bwaInputReadsp"])+int(c[2])

                if "Supplementary p" in line:
                    c =line.split("|")
                    if "bwaInputReadsuplp" not in context:
                        context["bwaInputReadsuplp"]=int(c[2])
                    else:
                        context["bwaInputReadsuplp"]=int(context["bwaInputReadsuplp"])+int(c[2])

                if "Properly Mapped p" in line:
                    c =line.split("|")
                    if "bwaMappedp" not in context:
                        context["bwaMappedp"]=int(c[2])
                    else:
                        context["bwaMappedp"]=int(context["bwaMappedp"])+int(c[2])

                if "Reads input w" in line:
                    c =line.split("|")
                    if "bwaInputReadsw" not in context:
                        context["bwaInputReadsw"]=int(c[2])
                    else:
                        context["bwaInputReadsw"]=int(context["bwaInputReadsw"])+int(c[2])

                if "Supplementary w" in line:
                    c =line.split("|")
                    if "bwaInputReadsuplw" not in context:
                        context["bwaInputReadsuplw"]=int(c[2])
                    else:
                        context["bwaInputReadsuplw"]=int(context["bwaInputReadsuplw"])+int(c[2])

                if "Properly Mapped w" in line:
                    c =line.split("|")
                    if "bwaMappedw" not in context:
                        context["bwaMappedw"]=int(c[2])
                    else:
                        context["bwaMappedw"]=int(context["bwaMappedw"])+int(c[2])

                if " Reads pair written  p" in line:
                    c =line.split("|")
                    if "bwaOutputp" not in context:
                        context["bwaOutputp"]=int(c[2])
                    else:
                        context["bwaOutputp"]=int(context["bwaOutputp"])+int(c[2])


                if "| Reads pair written w" in line:
                    c =line.split("|")
                    if "bwaOutputw" not in context:
                        context["bwaOutputw"]=int(c[2])
                    else:
                        context["bwaOutputw"]=int(context["bwaOutputw"])+int(c[2])

                if "Reads input f" in line:
                    c =line.split("|")
                    if "flashInput" not in context:
                        context["flashInput"]=int(c[2])
                    else:
                        context["flashInput"]=int(context["flashInput"])+int(c[2])


                if "Combined pairs" in line:
                    c =line.split("|")
                    if "flashCombined" not in context:
                        context["flashCombined"]=int(c[2])
                    else:
                        context["flashCombined"]=int(context["flashCombined"])+int(c[2])


                if "Uncombined pairs" in line:
                    c =line.split("|")
                    if "flashUncombined" not in context:
                        context["flashUncombined"]=int(c[2])
                    else:
                        context["flashUncombined"]=int(context["flashUncombined"])+int(c[2])


                if "Number of output sequences" in line:
                    c =line.split("|")
                    if "flashOutput" not in context:
                        context["flashOutput"]=int(c[2])
                    else:
                        context["flashOutput"]=int(context["flashOutput"])+int(c[2])


                if "Reads input a" in line:
                    c =line.split("|")
                    if "num_seq" not in context:
                        context["num_seq"]=int(c[2])
                    else:
                        context["num_seq"]=int(context["num_seq"])+int(c[2])

                if "Supplementary a" in line:
                    c =line.split("|")
                    if "ass_supl" not in context:
                        context["ass_supl"]=int(c[2])
                    else:
                        context["ass_supl"]=int(context["ass_supl"])+int(c[2])

                if " Reads mapped" in line:
                    c =line.replace(' (','|').split("|")
                    if "mapped" not in context:
                        context["mapped"]=int(c[2])
                    else:
                        context["mapped"]=int(context["mapped"])+int(c[2])


                if "Reads input b" in line:
                    c =line.split("|")
                    if "to_blast" not in context:
                        context["to_blast"]=int(c[2])
                    else:
                        context["to_blast"]=int(context["to_blast"])+int(c[2])


                if "Reads with blast hit" in line:
                    c =line.replace(' (','|').split("|")
                    if "blast_hits" not in context:
                        context["blast_hits"]=int(c[2])
                    else:
                        context["blast_hits"]=int(context["blast_hits"])+int(c[2])


                if "Contigs input" in line:
                    c =line.split("|")
                    if "to_blast_contigs" not in context:
                        context["to_blast_contigs"]=int(c[2])
                    else:
                        break

                if "Contigs with blast hit" in line:
                    c =line.split("|")
                    if "blast_hits_contigs" not in context:
                        context["blast_hits_contigs"]=c[2]
                    else:
                        break

context["flashCombinedpc"]=round((int(context["flashCombined"])*100)/(int(context["flashInput"])),2)
context["bwaMappedpcp"]=round((int(context["bwaMappedp"])*100)/(int(context["bwaInputReadsp"])+int(context["bwaInputReadsuplp"])),2)
context["bwaMappedpcw"]=round((int(context["bwaMappedw"])*100)/(int(context["bwaInputReadsw"])+int(context["bwaInputReadsuplw"])),2)
context["mappedpc"]=round((int(context["mapped"])*100)/(int(context["num_seq"])+int(context["ass_supl"])),2)
context["blast_hitspc"]=round((int(context["blast_hits"])*100)/(int(context["to_blast"])),2)

template = """
# SeekViralReads results for sample **_"Aedes"_**

---

# Adapters removal by cutadapt

## command line:


| variable | value |
| --- | --- |
| Total read pairs processsed R1| {CutA_1_total} |
| Total read pairs processsed R2| {CutA_2_total} |
| R1 with adapters |{CutA_1_R1_adap}|
| R2 with adapters | {CutA_1_R2_adap} |
| R1 total length before | {CutA_1_LR1} pb |
| R1 total length after | {CutA_1_LR1A} pb |
| R2 total length before | {CutA_1_LR2} pb |
| R2 total length after | {CutA_1_LR2A} pb |
| Reads 1 written | {CutA_1_written} |
| Reads 2 written | {CutA_2_written} |
| Reads pairs written | {CutA_1_2_written} |

## Quality trimming by cutadapt

### command line:



| variable | value |
| --- | --- |
| Total reads pairs processsed | {CutA_2_1_total} |
| Total reads pairs processsed2 | {CutA_2_2_total} |
| R1 total length before2 | {CutA_2_LR1} pb |
| R1 total length after2 | {CutA_2_LR1A} pb |
| R2 total length before2 | {CutA_2_LR2} pb |
| R2 total length after2 | {CutA_2_LR2A} pb |
| Reads 1 written2 | {CutA_2_1_written} |
| Reads 2 written2 | {CutA_2_2_written} |


## Mapping pairs on contaminent sequences

| variable | value |
| --- | --- |
| Reads input | {bwaInputReadsp} |
| Supplementary | {bwaInputReadsuplp} |
| Properly Mapped | {bwaMappedp} ({bwaMappedpcp}%) |
| Reads pair written | {bwaOutputp} |


## Mapping widows on contaminent sequences

| variable | value |
| --- | --- |
| Reads input | {bwaInputReadsw} |
| Supplementary | {bwaInputReadsuplw} |
| Properly Mapped | {bwaMappedw} ({bwaMappedpcw}%) |
| Reads written | {bwaOutputw} |


## Merging overlapping reads with flash

| variable | value |
| --- | --- |
| Reads input | {flashInput}  |
| Combined pairs | {flashCombined} ({flashCombinedpc}%) |
| Uncombined pairs | {flashUncombined} |
| Number of output sequences | {flashOutput}  |


##Mapping sequences on assembly
| variable | value |
| --- | --- |
| Reads input | {num_seq} |
| Supplementary | {ass_supl} |
| Reads mapped | {mapped} ({mappedpc}%)|



##Blast on nt_vir
| variable | value |
| --- | --- |
| Reads input b  | {to_blast} |
| Reads with blast hit | {blast_hits} ({blast_hitspc}%) |
| Contigs input  | {to_blast_contigs} |
| Contigs with blast hit | {blast_hits_contigs} |
"""


fileout.write(template.format(**context))
