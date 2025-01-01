# This will take a fasta file and output the longest coding sequence (ATG->STOP)
#

from collections import defaultdict
import sys

code = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
        'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
        'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
        }

def help():
    return '\nUSE AS FOLLOWS:\npython longestcds_nonredundant.py <input_file> > <output_filename>\n\n'

def find_longest_coding(sequence):
    seq1 = [sequence[_:_+3].upper() for _ in range(0,3*(len(sequence)/3),3)]
    seq2 = [sequence[_:_+3].upper() for _ in range(1,3*(len(sequence)-1)/3,3)]
    seq3 = [sequence[_:_+3].upper() for _ in range(2,3*(len(sequence)-2)/3,3)]
    longest = ''
    for seq in [seq1,seq2,seq3]:
        for pos,codon in enumerate(seq):
            if codon == 'ATG':
                newCDS = 'ATG'
                for pos2,codon2 in enumerate(seq[pos+1:]):
                    newCDS += codon2
                    if code[codon2] == '*': break
                if len(newCDS) > len(longest):
                    longest = newCDS
    return longest

def find_longest_nonredundant_isoforms(input_file):
    seqdic = defaultdict(str)
    with open(input_file) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line.strip()[1:]
            elif line.strip():
                seqdic[name] += line.strip()
    new_seqdic = {}
    seqlist = set()
    for seq in seqdic:
        newseq = find_longest_coding(seqdic[seq])
        if newseq not in seqlist:
            new_seqdic[seq] = newseq
            seqlist.add(newseq)
    for seq in new_seqdic:
        sys.stdout.write('>{0}\n'.format(seq))
        sys.stdout.write('{0}\n'.format('\n'.join(new_seqdic[seq][_:_+50] for _ in range(0,len(new_seqdic[seq]),50))))
    return

if __name__ == '__main__':
    sys.stderr.write(help())
    file = sys.argv[1:][0]
    find_longest_nonredundant_isoforms(input_file=file)


