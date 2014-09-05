#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

import sys
Entrez.email = "adina@anl.gov"

gbank = SeqIO.parse(open(sys.argv[1], "r"), "genbank")
output_handle=open("test.fa","w")

info = []
n = 1

def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def get_sequence_from_entrez(id_gen, strand, start, end):
    handle = Entrez.efetch(db="nucleotide", id=id_gen, rettype="fasta", retmode="text")
    record2 = SeqIO.parse(handle, "fasta")
    for gb_record2 in record2:
        if strand == 1:
            seq = gb_record2.seq[start:end]
        elif strand == -1:
            seq = gb_record2.seq[start:end]
        return seq

for genome in gbank :
    id = genome.name
    for feature in genome.features:
        if(feature.type == "CDS"):
            if feature.qualifiers.has_key('coded_by'):
                if feature.qualifiers['coded_by'][0].startswith("complement"):
                    dat_pre = feature.qualifiers['coded_by'][0].split("complement")[1][1:-1]
                    dat = dat_pre.split(':')
                    orig_genome = dat[0]
                    coord = dat[1].split('..')
                    start = int(coord[0]) - 1
                    end = int(coord[1])
                    strand = -1
                    nuc_seq = get_sequence_from_entrez(orig_genome, strand, start, end)
                    print ">" + orig_genome + "_" +"compl" + "_" + str(start) + "_" + str(end) + "_" + id
                    print nuc_seq
                else:
                    dat = feature.qualifiers['coded_by'][0].split(':')
                    orig_genome = dat[0]
                    coord = dat[1].split('..')
                    start = int(coord[0]) - 1
                    end = int(coord[1]) 
                    strand = 1
                    nuc_seq = get_sequence_from_entrez(orig_genome, strand, start, end)
                    print ">" + orig_genome + "_" + str(start) + "_" + str(end) + "_" + id
                    print nuc_seq

