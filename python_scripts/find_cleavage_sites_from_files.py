
import glob
import os
import sys
from collections import defaultdict
def find_proteins_with_cleavage_sites():
    sys.stdout.write("#\n")
    file_basenames = set()
    for i in glob.glob('./results/*'):
        file_basenames.add(os.path.splitext(i)[0].split('.')[1]+'.')
    for file in file_basenames:
        cleavage_sites = defaultdict(int)
        for posterior_file in glob.glob(file+"*"):
            with open(posterior_file) as fh:
                for line in fh:
                    cleavage_sites[line.strip().split('\t')[0]] += 1
        written = False
        for site in cleavage_sites:
            if cleavage_sites[site] == 4:
                if not written:
                    sys.stdout.write(file.split('/')[-1]+'\n')
                    written = True
                sys.stdout.write(site+'\n')
        if written:
            sys.stdout.write("#\n")
            written = False
if __name__ == "__main__":
    find_proteins_with_cleavage_sites()
