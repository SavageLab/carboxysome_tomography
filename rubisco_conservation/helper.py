import pandas as pd
from Bio import SeqIO, AlignIO, Seq
import numpy as np

def parse_uclust(infile,fasta,outfasta,outfile):
    header = ['Type','Cluster','Size','%Id','Strand','Qlo','Tlo','Alignment','Query','Target']
    uclust = pd.read_csv(infile, sep='\t', names=header, index_col=False)
    uclust.loc[uclust['Target'] == '*','Target'] = uclust.loc[uclust['Target'] == '*','Query']
    uclust = uclust[uclust['Type'] !='S']
    centroids = uclust[uclust['Type']=='C']
    c_list = centroids.Target.values
    c_list = [c.split(" ")[0] for c in c_list]

    sequences = []
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in c_list:
            sequences.append(record)

    seq2 = []
    seq3 = []
    for i,record in enumerate(sequences):
        if not record.id in seq2: 
            seq2.append(record.id)
            seq3.append(record)
    with open(outfasta, "w") as output_handle:
        SeqIO.write(seq3, output_handle, "fasta")