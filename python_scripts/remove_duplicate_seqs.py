# remove duplicate sequences from the msa used to make a pssm
# use as cat <input_file> | python remove_duplicate_seqs.py > <output_file_dedup>

import sys
def remove_dups():
    seq_set = set()
    for line in sys.stdin:
        seq_set.add(line.strip())
    seq_set = sorted(list(seq_set))
    for element in seq_set:
        sys.stdout.write(element+"\n")

if __name__ == '__main__':
    remove_dups()
