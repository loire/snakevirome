import sys, os
from collections import Counter
import pandas as pd
import csv

#Take 3 input Blast resum , Mat files , lineage
#Output two files Coverage  of contigs for each samples and coverage for lineages.

#Files with the blast infos for each seq/contigs. Output of the rule Join_seq_acc_taxo
Seq_hit_info_file= sys.argv[1]
seq_hit_info = open(Seq_hit_info_file, 'r')

#Path to the directory CountsMapping_all/. Outputs of the rule Quantify_contigs_coverage
Path_file_mat= sys.argv[2]
file_mat= os.listdir(Path_file_mat)

#File with the taxonomic lineage for each tax_ids. Output of the rule get_lineage_from_taxids
Taxonomy_file=sys.argv[3]
Taxonomy=open(Taxonomy_file, 'r')

output_coverge_contigs= sys.argv[4]
output_coverge_taxonomy=sys.argv[5]

### Produce a dic with Tax_ids as key and list of contig/seq corresponding as value
Dic_merge_taxids={}
for hits in seq_hit_info:
	Dic_best_hit={}
	list_hit=hits.split()
	if list_hit[0] not in Dic_best_hit :
		Dic_best_hit[list_hit[0]]=list_hit[5]
	else :
		continue
	if Dic_best_hit[list_hit[0]] not in Dic_merge_taxids:
		Dic_merge_taxids[Dic_best_hit[list_hit[0]]]=[list_hit[0]]
	elif list_hit[0] not in Dic_merge_taxids[Dic_best_hit[list_hit[0]]] :
		Dic_merge_taxids[Dic_best_hit[list_hit[0]]].append(list_hit[0])
	else:
		continue
del Dic_best_hit
del list_hit
seq_hit_info.close()
###Produce a csv file with the coverage of contigs and sequences for each samples.
###Take in input all the mat files from the  rule Quantify_contigs_coverage, map and unmapped files from the mapping of all cleaned sequences on contigs.

#Dictionnary with Sample(file-name) as key and subdic with seq/contig id as key an count as value. For sequence the widows count for 1 and the paired sequence have the same id and the value is 2
Dic_count={}
for files in file_mat:
	if files.endswith(".mat"):
		f=open(Path_file_mat+files, 'r')
		files_name=files.split("_R_1_2_counts_")[0]
		for line in f :
			line=line.replace('\n','').split()
			if files_name not in Dic_count:
				Dic_count[files_name]={}
				Dic_count[files_name].update({line[1]:line[0]})
			else:
				Dic_count[files_name].update({line[1]:line[0]})
	f.close()
data_frame = pd.DataFrame(Dic_count).fillna(0)
data_frame.to_csv(output_coverge_contigs)
del Dic_count
del line

# Dictionnary with seq/contig id as key and list of count values for each samples .
Dic_count_by_samples={}
data_frame_to_str=data_frame.to_csv()
list_count=data_frame_to_str.split("\n")
samples=list_count[0].split(",")
for line in list_count[1:]:
    liste_count=line.split(',')
    Dic_count_by_samples[liste_count[0]]=[]
    for count in liste_count[1:]:
        Dic_count_by_samples[liste_count[0]].append(int(count))
del list_count
del data_frame_to_str

Dic_count_taxids={}
for tax_ids in Dic_merge_taxids:
    Dic_count_taxids[tax_ids]=[]
    for ids in Dic_merge_taxids[tax_ids]:
        ids=ids.split('/', 1)[0] #Remove the /1 or /2 from the ids of sequences
        Dic_count_taxids[tax_ids].append(Dic_count_by_samples[ids])
del Dic_count_by_samples
#Produce th final dic with the sum of count from each list
Sum_count_per_taxids={}
for tax_ids in Dic_count_taxids :
	l=[]
	for listes in Dic_count_taxids[tax_ids]:
		l.append(listes)
		Sum_count_per_taxids[tax_ids]=list(map(sum, zip(*l)))
del Dic_count_taxids
del l
#Merge information of coverage by tax_ids and lineage to give coverage for each sample for taxonomy
ranks=Taxonomy.readline().replace('\n','').split("|")
lineage_coverage=[ranks+ samples[1:]]

for lineages in Taxonomy:
	list_lineage=lineages.replace('\n','').split('|')
	lineage_coverage.append(list_lineage + Sum_count_per_taxids[list_lineage[0]])
with open(output_coverge_taxonomy, "w") as f:
    writer = csv.writer(f)
    writer.writerows(lineage_coverage)
Taxonomy.close()
