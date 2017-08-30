# #!/bin/python

import csv
from ete3 import NCBITaxa
import sys
from collections import Counter

ncbi = NCBITaxa()

def get_desired_ranks(taxids, desired_ranks):
    lineage = ncbi.get_lineage(taxids)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    lineage2name = ncbi.translate_to_names(lineage)
    Dic_lineage2name=dict(zip(lineage, lineage2name))
    for rank in ranks2lineage :
        ranks2lineage[rank]=Dic_lineage2name[ranks2lineage[rank]]
    ranks2lineage['tax_id']=taxids
    return {'{}'.format(rank): ranks2lineage.get(rank, '<not_present>') for rank in desired_ranks}
    print({'{}'.format(rank): ranks2lineage.get(rank, '<not_present>') for rank in desired_ranks})

def main(taxids, desired_ranks, path):
    with open(path, 'w') as csvfile:
        fieldnames = ['{}'.format(rank) for rank in desired_ranks]
        writer = csv.DictWriter(csvfile, delimiter='|', fieldnames=fieldnames)
        writer.writeheader()
        for taxid in taxids:
            writer.writerow(get_desired_ranks(taxid, desired_ranks))
    csvfile.close()


if __name__ == '__main__':
    taxids=sys.stdin.read().split(',')
    desired_ranks = ['tax_id', 'superkingdom', 'class', 'order', 'family', 'subfamily','genus', 'species', 'no rank']
    path = sys.argv[1]
    main(taxids, desired_ranks, path)
