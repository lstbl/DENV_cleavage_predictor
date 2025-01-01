#this program will take a gff file and fasta file of all protein coding sequences and output a fasta file with only
#the longest isoform protein coding sequences
# python protein_longest_isoforms.py <gff_file> <fasta_file_with_protein_sequence_from_NCBI>

import gzip
import re
from collections import defaultdict
import sys
#parse gff file to get refseq proteins
def return_gff_protein_dic(gff_file):
    if gff_file[-3:] == ".gz":
        openfn = gzip.open
    else:
        openfn = open
    with openfn(gff_file) as fh:
        gene_dic = defaultdict(set)
        for line in fh:
            if line.startswith('#'): continue
            line = line.strip().split('\t')
            if line[1] == "BestRefSeq" and line[2] == "CDS":
                proteinID = re.findall('protein_id=(NP_[^\n;]*)',line[-1])
                geneID = re.findall('GeneID:(\d+)',line[-1])
                if proteinID and geneID:
                    proteinID, geneID = proteinID[0], geneID[0]
                    gene_dic[geneID].add(proteinID)
    return gene_dic


# Run twice through the protein fasta file and pull out longest isoforms
# first time through is in order to figure out the longest isoform, second time through is to write out the longest
# isoform to a file
# each isoform will be written to a different file in order to allow for phobius to work most efficiently
def find_longest_isoforms(gff_file, fasta_file):
    refseq_dic = return_gff_protein_dic(gff_file=gff_file)
    if fasta_file[-3:] == ".gz":
        openfn = gzip.open
    else:
        openfn = open
    with openfn(fasta_file) as fh:
        current = ""
        protein_len = defaultdict(int)
        for line in fh:
            if line.startswith(">NP"):
                line = line[1:].strip().split()
                current = line[0]
            else:
                protein_len[current] += len(line.strip())
    longest_isoform = set()
    for geneID in refseq_dic:
        longest = ''
        length = 0
        for refseq in refseq_dic[geneID]:
            assert refseq in protein_len
            if protein_len[refseq] > length:
                longest = refseq
                length = protein_len[refseq]
        longest_isoform.add(longest)
    with openfn(fasta_file) as fh:
        write = False
        for line in fh:
            if line.startswith('>'):
                write = False
            if line.startswith('>NP'):
                write = True
                current = line[1:].strip().split()[0]
                if current in longest_isoform:
                    wfile = open(current+".fa",'w')
                    wfile.write(">"+current+'\n')
                    continue
                else:
                    write = False
            elif current in longest_isoform and write:
                wfile.write(line)
    wfile.close()
    return


if __name__ == "__main__":
    files = sys.argv[1:]
    find_longest_isoforms(*files)



