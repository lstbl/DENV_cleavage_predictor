
from collections import defaultdict
import numpy as np
import glob
import sys
from sklearn.linear_model import LogisticRegression
from sklearn import svm


def calculate_percentages(true_positives="positive_control_seqs_for_training_8AA.txt"):
    AAs = "GPAVLIMCFYWHKRQNEDST"
    seqs = []
    with open(true_positives) as ofile:
        first =True
        seqlen=0
        for line in ofile:
            line = line.strip()
            if line:
                seqs.append([_ for _ in line])
                if not first:
                    assert len(line) == seqlen
                    first = False
                seqlen = len(line)
    pos_dic = defaultdict(int)
    # give a pseudocount of 1
    for AA in AAs:
        for pos in range(seqlen):
            pos_dic[pos,AA] += 1
    for colpos in range(seqlen):
        for rowpos in range(len(seqs)):
            pos_dic[colpos,seqs[rowpos][colpos]] += 1
    # make pos_dic into percentages
    for entry in pos_dic:
        pos_dic[entry] = pos_dic[entry]/float(len(seqs))
    return dict(pos_dic)

def chemistry_info(chemistry_file="AA_properties.txt"):
    with open(chemistry_file) as fh:
        AA_chemistry = {}
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split('\t')
                AA_chemistry[line[0]] = [float(_) for _ in line[1:]]
    return AA_chemistry

def sequence_info(sequence="SEQUENCE", positional_percentages = {("percentages_dict",) : 0.1}, chemistry_dic = {("chemistry_info",) : [0,1,1,0]}):
    seq_info = []
    for pos,AA in enumerate(sequence):
        seq_info.append(positional_percentages[pos,AA])
        seq_info += chemistry_dic[AA]
    return seq_info


def targets_from_phobius(phobius_file='./posteriors/human_STING_posteriors.txt', test_len=8):
    with open(phobius_file) as ofile:
        sequence = ''
        membrane = []
        cytoplasmic = []
        for line in ofile:
            if line.startswith("#"): continue
            line = line.strip().split('\t')
            if line:
                sequence += line[1]
                membrane.append(float(line[4]))
                cytoplasmic.append(float(line[2]))
    end = 0
    for pos in range(len(sequence)-1,-1,-1):
        if membrane[pos] >= 0.5:
            end = pos
            break
    if end == 0:
        return [],[]
    if end + 20 < len(sequence):
        end += 21
    else:
        end = len(sequence)
    boarders = []
    current = []
    for pos in range(end):
        if current or pos == end-1:
            if cytoplasmic[pos] < 0.5 or pos == end-1:
                current.append(pos)
                boarders.append(tuple(current))
                current = []
        elif not current:
            if cytoplasmic[pos] > 0.5:
                current.append(pos)
    if len(boarders[-1]) == 1:
        boarders = boarders[:-1]
    sequences = []
    for boarder in boarders:
        if boarder[1]-boarder[0] < test_len: continue
        sequences.append(sequence[boarder[0]:boarder[1]+1])
    return sequences, boarders

def split_sequence(sequence = "", split_seq = "", length = 8):
    if len(sequence) < len(split_seq):
        return [],[]
    split_positions = [0]
    for pos in range(len(sequence)-len(split_seq)):
        if sequence[pos:pos+len(split_seq)] == split_seq and pos > split_positions[-1]:
            split_positions += [pos, pos+len(split_seq)]
    if split_positions[-1] == 0:
        return [sequence], [(0,len(sequence)-1)]
    split_positions.append(len(sequence))
    assert len(split_positions) % 2 == 0
    split_sequences = [sequence[split_positions[_]:split_positions[_+1]] for _ in range(0, len(split_positions),2)]
    split_seq_positions = [(split_positions[_],split_positions[_+1]) for _ in range(0, len(split_positions),2)]
    return_sequences = []
    return_positions = []
    for pos,seq in enumerate(split_sequences):
        if len(seq) >= length:
            return_sequences.append(split_sequences[pos])
            return_positions.append(split_seq_positions[pos])
    return return_sequences, return_positions


def create_sequences_from_targets(sequence_list = [] , positions_list = [],length=8, avoid = []):
    if avoid:
        for avoided_seq in avoid:
            new_sequence_list = []
            new_positions_list = []
            for pos,sequence in enumerate(sequence_list):
                start = positions_list[pos][0]
                seqs, positions = split_sequence(sequence = sequence, split_seq = avoided_seq, length = length)
                new_sequence_list += seqs
                new_positions_list += [(positions[_][0]+start,positions[_][1]+start) for _ in range(len(positions))]
        sequence_list, positions_list = new_sequence_list, new_positions_list
    new_sequence_list = []
    new_positions_list = []
    for pos, sequence in enumerate(sequence_list):
        start = positions_list[pos][0]
        for startpos in range(len(sequence)-length+1):
            new_sequence_list.append(sequence[startpos:startpos+length])
            new_positions_list.append(( start + startpos, start + startpos + length))
    return new_sequence_list, new_positions_list

def create_training_data(positive_training_example = 'positive_control_seqs_for_training_8AA_dedup.txt',
                         phobius_files = './posteriors', percentage_dic = calculate_percentages(),
                         chemistry_dic=chemistry_info(), training_output_prefix = "training_examples"):
    # Create positive training examples
    positive_seqs = []
    with open(positive_training_example) as fh:
        for line in fh:
            line = line.strip()
            if line:
                positive_seqs.append(line)
    with open(training_output_prefix+".positives",'w') as wfile:
        wfile.write("#Sequence\t{0}\n".format('\t'.join('pos'+str(_)+' percentage\tpos'+str(_)+' size\tpos'+str(_)+
                                                      ' charge\tpos'+str(_)+' hydrophobicity\tpos'+str(_)+' aromatic'
                                                        for _ in range(len(positive_seqs[0])))))
        for seq in positive_seqs:
            wfile.write(seq +'\t'+'\t'.join(str(_) for _ in
                                  sequence_info(sequence=seq, positional_percentages=percentage_dic,
                                                chemistry_dic=chemistry_dic))+'\n')

    # Create negative training examples
    files = glob.glob(phobius_files+'/*')
    negative_seqs = []
    for file in files:
        sequences, boarders = targets_from_phobius(phobius_file=file)
        sequences, boarders = create_sequences_from_targets(sequence_list=sequences, positions_list=boarders,
                                                            avoid=positive_seqs)
        negative_seqs += sequences
    assert len(negative_seqs[0]) == len(positive_seqs[0])
    with open(training_output_prefix+".negatives",'w') as wfile:
        wfile.write("#Sequence\t{0}\n".format('\t'.join('pos'+str(_)+' percentage\tpos'+str(_)+' size\tpos'+str(_)+
                                                      ' charge\tpos'+str(_)+' hydrophobicity\tpos'+str(_)+' aromatic'
                                                        for _ in range(len(negative_seqs[0])))))
        negative_seqs = set(negative_seqs)
        for seq in negative_seqs:
            wfile.write(seq + '\t' + '\t'.join(str(_) for _ in
                                               sequence_info(sequence=seq, positional_percentages=percentage_dic,
                                                             chemistry_dic=chemistry_dic)) + '\n')

def prepare_training_data(training_output_prefix = "training_examples"):
    XOR_list = []
    data_list = []
    with open(training_output_prefix+".positives") as ofile:
        for line in ofile:
            if line.startswith("#"): continue
            XOR_list.append(True)
            line = line.strip().split('\t')
            data_list.append([float(_) for _ in line[1:]])
    with open(training_output_prefix+".negatives") as ofile:
        for line in ofile:
            if line.startswith("#"): continue
            XOR_list.append(False)
            line = line.strip().split('\t')
            data_list.append([float(_) for _ in line[1:]])
    return np.array(XOR_list), np.array(data_list)

def test_sequence_for_cleavage(percentage_dic=calculate_percentages(),
                                   chemistry_dic=chemistry_info(),
                                   phobius_file='./posteriors/human_STING_posteriors.txt',
                                   sequence_list=[],
                                   return_probs=False,
                                   fit_model=LogisticRegression(C=1).fit([[1], [1]], [3, 1]),
                                   return_bool = False,
                                   probability=False):


    if sequence_list:
        sequences = sequence_list
    else:
        sequences, boarders = create_sequences_from_targets(*targets_from_phobius(phobius_file=phobius_file))
    seq_data = []
    for seq in sequences:
        seq_data.append(sequence_info(sequence=seq, positional_percentages=percentage_dic, chemistry_dic=chemistry_dic))
    test_results = fit_model.predict(seq_data)
    if probability:
        probs = fit_model.predict_proba(seq_data)
    if return_probs:
        probs_to_return = []
        for pos, result in enumerate(test_results):
            probs_to_return.append(probs[pos][1])
        return probs_to_return
    if return_bool:
        return test_results
    for pos, result in enumerate(test_results):
        if sequence_list and probability:
            sys.stdout.write('{0} {1} {2}\n'.format(sequences[pos], probs[pos], test_results[pos]))
        elif sequence_list:
            sys.stdout.write('{0} {1}\n'.format(sequences[pos], test_results[pos]))
        elif result and probability:
            sys.stdout.write('{0} {1} {2}\t{3}\t{4}\n'.format(boarders[pos][0], sequences[pos], boarders[pos][1], probs[pos],
                                                   test_results[pos]))
        elif result:
            sys.stdout.write('{0} {1} {2}\n'.format(boarders[pos][0], sequences[pos], boarders[pos][1]))
    return

# this function will be used to test all human genes on the cluster. It will use the
def test_gene_for_cleavage(percentage_dic=calculate_percentages(),
                           chemistry_dic=chemistry_info(),
                           phobius_files='./posteriors'):
    XOR_list, data_list = prepare_training_data()
    SVM = svm.SVC(kernel='rbf',C=10).fit(data_list,XOR_list)
    LR = LogisticRegression(C=10).fit(data_list,XOR_list)
    for phobius_file in glob.glob(phobius_files+'/*'):
        sequences, boarders = create_sequences_from_targets(*targets_from_phobius(phobius_file=phobius_file))
        if not sequences: continue
        seq_data = []
        for seq in sequences:
            seq_data.append(sequence_info(sequence=seq, positional_percentages=percentage_dic, chemistry_dic=chemistry_dic))
        test_results = SVM.predict(seq_data)
        probs = LR.predict_proba(seq_data)
        to_write = ''
        for pos, result in enumerate(test_results):
            if result:
                to_write += '{0} {1} {2}\t{3}\n'.format(boarders[pos][0]+1, sequences[pos], boarders[pos][1]+1,
                                                        probs[pos])
        if to_write:
            outfile = phobius_file.split('/')[-1].split('.')[0]
            with open('./results/' + outfile + '.result', 'w') as wfile:
                wfile.write(to_write)


if __name__ == "__main__":
    test_gene_for_cleavage()



