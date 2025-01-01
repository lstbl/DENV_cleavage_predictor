# this will take the directory that the phobius output files are in and output the sequence of the protein

import glob
def get_seq(phobius_directory,genbank_accession):
    if not phobius_directory[-1] == '/':
        phobius_directory += '/'
    file = glob.glob(phobius_directory+genbank_accession+'*')
    assert file
    file = file[0]
    seq = ''
    with open(file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line:
                seq += line[1]
    return seq

