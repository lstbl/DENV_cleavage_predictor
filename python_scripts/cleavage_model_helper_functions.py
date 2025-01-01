from collections import defaultdict
import glob
import numpy as np
import sys
from sklearn import svm
from sklearn.linear_model import LogisticRegression
import os


# parse the posterior file output from the HMM "Phobius" and return a dictionary of overlapping amino acid sequences of
# length N (default 30), with the properties of the cleavage site (offset M from the right) in a tuple of the form
# (non-cytoplasmic, cytoplasmic, membrane). The dictionary will look something like this:
# {"ARDCGGEQRGGDAERQF":(0.455,0.545,0.003) ... }
def parse_phobius_posteriors(posterior_filename="filename",length=30,cleavage_offset=(13,)):
    with open(posterior_filename) as ofile:
        posterior_info = []
        for line in ofile:
            if line.startswith("#"): continue
            line = line.strip().split('\t')
            if line:
                posterior_info.append(tuple(line[1:5]))
    motif_info_dic = dict()
    avoid = set()
    for pos in range(len(posterior_info)-length):
        seq = ''.join(_[0] for _ in posterior_info[pos:pos+length])
        if seq in motif_info_dic:
            del motif_info_dic[seq]
            avoid.add(seq)
        if seq in avoid:
            continue
        result = []
        for offset in cleavage_offset:
            result += [float(_) for _ in posterior_info[pos+offset-1][1:4]]
        motif_info_dic[seq] = tuple(result)
    return motif_info_dic



# Parse a file with a MSA of known cleavage sites to determine the fraction of any specific amino acid at a given
# poistion. You will have to modify this input based on how much info you plan on including in your model. Currently
# I am putting in a 30 amino acid sequence, with the cleavage site between positions 13 and 14.
# Use a pseudocount of 1 for every amino acid at every position
# The ouput will be a dicionary that will take a tuple input and output the percent identity of that AA at a given
# position. The positions will be zero-based to keep with the python format
def calculate_percentages(true_positives="positive_control_seqs_for_training.txt"):
    AAs = "GPAVLIMCFYWHKRQNEDST"
    seqs = []
    with open(true_positives) as ofile:
        first =True
        seqlen=0
        for line in ofile:
            line = line.strip()
            seqs.append([_ for _ in line])
            if not first:
                assert len(line) == seqlen
                first = False
            seqlen = len(line)
    pos_dic = defaultdict(int)
    for AA in AAs:
        for pos in range(seqlen):
            pos_dic[pos,AA] += 1
    for colpos in range(seqlen):
        for rowpos in range(len(seqs)):
            pos_dic[colpos,seqs[rowpos][colpos]] += 1
    # make pos_dic into percentages
    for entry in pos_dic:
        pos_dic[entry] = pos_dic[entry]/float(len(seqs))
    return pos_dic

# function to merge dictionaries:
def merge_two_dicts(x,y):
    z = {}
    for entry in x:
        if entry not in y:
            z[entry] = x[entry]
    for entry in y:
        if entry not in x:
            z[entry] = y[entry]
    return z

# Make dictionary with all possible AA combinations from posterior files.
def collect_posterior_info(folder = "folder_with_phobius_files",length=30,offset=(13,)):
    big_dic = {}
    for file in glob.glob(folder+'/*'):
        big_dic = merge_two_dicts(big_dic,parse_phobius_posteriors(posterior_filename=file,length=length,cleavage_offset=offset))
    return big_dic

# Make positive and negative training data
def create_training_data(posterior_file ="./posteriors", positives ="positive_control_seqs_for_training.txt", length=30,
                         offset=(13,), chemistry=False, posteriors=True):
    percentage_dic = calculate_percentages(true_positives=positives)
    if posteriors:
        sequence_info = collect_posterior_info(folder=posterior_file, length=length, offset=offset)
    positive_seqs = set()
    if chemistry:
        with open(chemistry) as fh:
            AA_chemistry = {}
            for line in fh:
                if line.startswith("#"): continue
                else:
                    line = line.strip().split('\t')
                    AA_chemistry[line[0]] = line[1:]
    with open(positives) as ofile:
        for line in ofile:
            positive_seqs.add(line.strip())
    with open('positive_controls.training','w') as wfile:
        wfile.write(
            "#sequence\t{0}\tcytoplasmic_posterior_prob\tnon-cytoplasmic_posterior_prob\ttransmembrane_posterior_prob\n".format(
                '\t'.join("pos" + str(_) for _ in range(1, len(list(positive_seqs)[0]) + 1))))
        for seq in positive_seqs:
            if not posteriors:
                line = seq + "\t"
                line += "\t".join(str(percentage_dic[pos, seq[pos]]) for pos in range(len(seq))) + "\t"
                if chemistry:
                    line += "\t".join("\t".join(_ for _ in AA_chemistry[AA]) for AA in seq) + "\t"
                wfile.write(line)
            elif posteriors:
                if seq in sequence_info:
                    line = seq+"\t"
                    line += "\t".join(str(percentage_dic[pos,seq[pos]]) for pos in range(len(seq)))+"\t"
                    if chemistry:
                        line += "\t".join("\t".join(_ for _ in AA_chemistry[AA]) for AA in seq)+"\t"
                    line += "\t".join(str(_) for _ in sequence_info[seq])+"\n"
                    wfile.write(line)
    with open("negative_controls.training",'w') as wfile:
        wfile.write(
            "#sequence\t{0}\tcytoplasmic_posterior_prob\tnon-cytoplasmic_posterior_prob\ttransmembrane_posterior_prob\n".format(
                '\t'.join("pos" + str(_) for _ in range(1, len(list(positive_seqs)[0]) + 1))))
        for seq in sequence_info:
            if seq not in positive_seqs:
                line = seq+"\t"
                line += "\t".join(str(percentage_dic[pos,seq[pos]]) for pos in range(len(seq)))+"\t"
                if chemistry:
                    line += "\t".join("\t".join(_ for _ in AA_chemistry[AA]) for AA in seq)+"\t"
                line += "\t".join(str(_) for _ in sequence_info[seq])+"\n"
                wfile.write(line)


# prepare training data produced by the "create_training_data" function
def prepare_training_data(positives = "positive_controls.training", negatives = "negative_controls.training"):
    XOR_list = []
    data_list = []
    with open(positives) as ofile:
        for line in ofile:
            if line.startswith("#"): continue
            XOR_list.append(True)
            line = line.strip().split('\t')
            data_list.append([float(_) for _ in line[1:]])
    with open(negatives) as ofile:
        for line in ofile:
            if line.startswith("#"): continue
            XOR_list.append(False)
            line = line.strip().split('\t')
            data_list.append([float(_) for _ in line[1:]])
    return np.array(XOR_list), np.array(data_list)



def test_sequence_for_cleavage(model = "trained_model",phobius_file = "<phobius_posterior_file>",seq_len=30,
cleavage_offset=(13,), chemistry = False,probability=False,suffix='result'):
    percentage_dic = calculate_percentages()
    posterior_info = []
    sequence_info = parse_phobius_posteriors(posterior_filename=phobius_file,length=seq_len,cleavage_offset=cleavage_offset)
    if chemistry:
        with open(chemistry) as fh:
            AA_chemistry = {}
            for line in fh:
                if line.startswith("#"): continue
                else:
                    line = line.strip().split('\t')
                    AA_chemistry[line[0]] = line[1:]
    sequence = ''
    with open(phobius_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            line = line.strip().split('\t')
            if line:
                sequence += line[1]
    cleavage_sites = []
    seqs_to_test = []
    os.system('mkdir -p results')
    for start in range(len(sequence)-seq_len):
        seq_data = []
        seq = sequence[start:start+seq_len]
        cleavage_sites.append(seq)
        if seq not in sequence_info: continue
        seq_data += [float(percentage_dic[pos,seq[pos]]) for pos in range(len(seq))]
        if chemistry:
            for AA in seq:
                seq_data += [float(_) for _ in AA_chemistry[AA]]
        seq_data += [float(_) for _ in sequence_info[seq]]
        seqs_to_test.append(seq_data)
    test_results = model.predict(seqs_to_test)
    if probability:
        probs = model.predict_proba(seqs_to_test)
    with open('./results/'+phobius_file.split('/')[-1]+'cleavagesites.{0}'.format(suffix),'w') as wfile:
        for pos,result in enumerate(test_results):
            if result and probability:
                wfile.write('{0} {1} {2}\t{3}\t{4}\n'.format(pos, cleavage_sites[pos], pos + len(cleavage_sites[pos]),
                                                             probs[pos], test_results[pos]))
            elif result:
                wfile.write('{0} {1} {2}\n'.format(pos, cleavage_sites[pos], pos+len(cleavage_sites[pos])))

if __name__ == "__main__":
    models = [svm.SVC(kernel='linear', C=1),svm.SVC(kernel='linear', C=10),svm.SVC(kernel='rbf', C=10),LogisticRegression(C=10)]
    phobius_file = sys.argv[1:][0]
#    create_training_data(offset=tuple(range(1, 31)), chemistry="AA_properties.txt")
    XOR_list, data_list = prepare_training_data()
    test_sequence_for_cleavage(model=svm.SVC(kernel='linear', C=1).fit(data_list,XOR_list), phobius_file=phobius_file, seq_len=30,
                               cleavage_offset=tuple(range(1,31)), chemistry='AA_properties.txt', probability=False, suffix='linearC1')

    test_sequence_for_cleavage(model=svm.SVC(kernel='linear', C=10).fit(data_list, XOR_list), phobius_file=phobius_file, seq_len=30,
                               cleavage_offset=tuple(range(1,31)), chemistry='AA_properties.txt', probability=False, suffix='linearC10')

    test_sequence_for_cleavage(model=svm.SVC(kernel='rbf', C=10).fit(data_list, XOR_list), phobius_file=phobius_file, seq_len=30,
                               cleavage_offset=tuple(range(1,31)), chemistry='AA_properties.txt', probability=False, suffix='rbfC10')

    test_sequence_for_cleavage(model=LogisticRegression(C=10).fit(data_list, XOR_list), phobius_file=phobius_file, seq_len=30,
                               cleavage_offset=tuple(range(1,31)), chemistry='AA_properties.txt', probability=True, suffix='LogRegC10')
