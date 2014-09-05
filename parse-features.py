#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

gbank = SeqIO.parse(open(sys.argv[1], "r"), "genbank")
output_handle=open("test.fa","w")

CDSs = []
n = 1
for genome in gbank :
    for feature in genome.features:
        if(feature.type == "CDS"):
            ID = feature.qualifiers['db_xref'][0] + "_" + str(n)
            desc = feature.qualifiers['locus_tag'][0]
            if feature.qualifiers.has_key('gene'):
                gene = feature.qualifiers['gene'][0]
            else:
                gene = "no_gene_name"
            seq = feature.extract(genome.seq)
            record = SeqRecord(seq, id=ID, description=desc+"_"+gene)
            CDSs.append(record)
            n += 1

SeqIO.write(CDSs, output_handle, "fasta")
output_handle.close()
