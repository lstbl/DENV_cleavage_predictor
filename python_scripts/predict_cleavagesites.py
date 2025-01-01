from collections import defaultdict
import numpy as np
import glob
import sys
import getopt
from sklearn.linear_model import LogisticRegression
from sklearn import svm
import pandas as pd

def create_linear_scaling_function(input_data):
    maxval = max(input_data)
    minval = min(input_data)
    # a(maxval) + b = 1
    # a(minval) + b = -1
    a = 2./(maxval-minval)
    b = 1. - a*maxval
    def scalingfn(inp):
        return round(a*inp + b,5)
    return scalingfn

'''
simulations = ['' for _ in range(1000)]
def find_choice(weighted_distance, choicevalue):
    total = sum(weighted_distance[_][1] for _ in range(len(weighted_distance)))
    cumulant = 0
    for pos,entry in enumerate(weighted_distance):
        cumulant += entry[1]/total
        if cumulant >=choicevalue:
            return entry[0]
for pos in range(8):
    weighted_dist = [(AA,-get_weighted_dist(AA,pos,df,percent_dic)) for AA in 'GPAVLIMCFYWHKRQNEDST']
    #weighted_dist = [(weighted_dist[_][0],weighted_dist[_][1]-min(weighted_dist,key=itemgetter(1))[1]) for _ in range(len(weighted_dist))]
    #weighted_dist = [(weighted_dist[_][0],np.exp(weighted_dist[_][1])) for _ in range(len(weighted_dist))]
    weighted_dist = [(weighted_dist[_][0], np.exp(weighted_dist[_][1])) for _ in range(len(weighted_dist))]
    for pos2 in range(len(simulations)):
        simulations[pos2] += find_choice(weighted_dist,np.random.random())
'''





def calculate_percentages(true_positives="positive_control_seqs_for_training_8AA.txt",normalize=True):
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

    # add every AA to every position:
    for AA in AAs:
        for pos in range(seqlen):
            pos_dic[pos,AA] += 0

    for colpos in range(seqlen):
        for rowpos in range(len(seqs)):
            pos_dic[colpos,seqs[rowpos][colpos]] += 1
    # make pos_dic into percentages
    data_points = []
    for entry in pos_dic:
        percentage = pos_dic[entry]/float(len(seqs))
        data_points.append(percentage)
        pos_dic[entry] = percentage
    if normalize:
        # normalize percentages to a [-1,1] scale
        scalefn = create_linear_scaling_function(data_points)
        for entry in pos_dic:
            pos_dic[entry] = scalefn(pos_dic[entry])

    return dict(pos_dic)

def chemistry_info(chemistry_file="AA_properties.txt"):
    with open(chemistry_file) as fh:
        AA_chemistry = {}
        AAs = ''
        for line in fh:
            if line.startswith("#"):
                line = line[1:].strip().split('\t')
                categories = line[1:]
                datas = [[] for _ in range(len(categories))]
            else:
                line = line.strip().split('\t')
                AA = line[0]
                AAs += AA
                values = [float(_) for _ in line[1:]]
                for pos,value in enumerate(values):
                    datas[pos].append(value)
                    AA_chemistry[AA,categories[pos]] = value
    assert len(AAs) == 20
    #normalize AA_chemistry values
    normalizing_fns = [create_linear_scaling_function(_) for _ in datas]
    for pos,category in enumerate(categories):
        scalefn = normalizing_fns[pos]
        for AA in AAs:
            AA_chemistry[AA,category] = scalefn(AA_chemistry[AA,category])
    return AA_chemistry

def AA_Grantham_dist(distfile='AA_Grantham_dists.txt'):
    df = pd.read_table(distfile,header=0,index_col=0)
    normalizing_fn = create_linear_scaling_function(df.get_values().flatten())
    for AA1 in df.columns:
        for AA2 in df.columns:
            df.loc[AA1,AA2] = normalizing_fn(df.loc[AA1,AA2])
    return df

def get_weighted_dist(AA = '', pos = 0, grantham_dist = pd.DataFrame(), AA_percent_sim = dict()):
    return sum(grantham_dist.loc[AA,cur_AA]*AA_percent_sim[pos,cur_AA] for cur_AA in grantham_dist.columns)

def get_training_info(file = './parameter_files/parameters.csv'):
    return pd.read_csv(file,header=0,index_col=0)


# must modify the AA_info dictionary to contain information in the AA_chemistry, grantham_dists, and percent similarity
# input
def sequence_info(sequence="SEQUENCE", granthamdist = pd.DataFrame(), AA_percent_sim = {}, AA_chemistry = {},parameters = pd.DataFrame()):
    seq_info = [sequence]
    for pos,AA in enumerate(sequence):
        for parameter in parameters.index:
            if pd.notnull(parameters.loc[parameter,str(pos)]):
                if parameter == 'grantham_dist':
                    seq_info.append(get_weighted_dist(AA,pos,granthamdist,AA_percent_sim))
                elif parameter == 'positional_percentage':
                    seq_info.append(AA_percent_sim[pos,AA])
                else:
                    assert (AA,parameter) in AA_chemistry
                    seq_info.append(AA_chemistry[AA,parameter])
    return seq_info


# generates target areas from a phobius file (outputs 2 lists, first is a list of strings (Amino Acids); second is a
# list of tuples (inclusive boarders of amino acids)
def targets_from_phobius(phobius_file='./posteriors/human_STING_posteriors.txt', test_len=8,cytoplasmic_only = True):
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
            if cytoplasmic_only:
                if cytoplasmic[pos] < 0.5 or pos == end-1:
                    current.append(pos)
                    boarders.append(tuple(current))
                    current = []
            else:
                if membrane[pos] >= 0.5 or pos == end-1:
                    current.append(pos)
                    boarders.append(tuple(current))
                    current = []
        elif not current:
            if cytoplasmic_only:
                if cytoplasmic[pos] > 0.5:
                    current.append(pos)
            else:
                if membrane[pos] < 0.5:
                    current.append(pos)
    if len(boarders[-1]) == 1:
        boarders = boarders[:-1]
    sequences = []
    new_boarders = []
    for boarder in boarders:
        if boarder[1]-boarder[0] < test_len: continue
        sequences.append(sequence[boarder[0]:boarder[1]+1])
        new_boarders.append(boarder)
    return sequences, new_boarders

# determines if there is a known positive sequence in a domain and cuts it out. If the remaining sequences are less than
# 'length' (should be almost always 8 AA long), then it removes them altogether. This is used to construct the negative
# training examples
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

# takes sequences generated by "targets from phobius" and returns a list of sequences of length 'length' (should be
# 8). This will also remove sequences that overlap with known cut sequences.
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

# takes the known 8AA cleavage sequences and generates data based on the "parameters" file
# note, the positive controls MUST be 8AA long
def create_training_data(positive_training_example = 'positive_training_8AA_dedup.txt',
                         phobiusfiles_from_posexamples ='./posteriors', percentage_dic = calculate_percentages(),
                         granthamDist = AA_Grantham_dist(), chemistry_dic=chemistry_info(),
                         training_output_prefix = "training_examples", length = 8, parameters = get_training_info()):
    # Create positive training examples
    positive_seqs = []
    with open(positive_training_example) as fh:
        for line in fh:
            line = line.strip()
            if line:
                assert len(line) == length
                positive_seqs.append(line)
    with open(training_output_prefix+".positives",'w') as wfile:
        header = '#Sequence\t'
        for pos in range(length):
            for parameter in parameters.index:
                if pd.notnull(parameters.loc[parameter,str(pos)]):
                    header += 'pos{0}_{1}\t'.format(pos,parameter)
        wfile.write(header[:-1]+'\n')
        for seq in positive_seqs:
            wfile.write(seq +'\t'+'\t'.join(str(_) for _ in
                                            sequence_info(sequence=seq, granthamdist=granthamDist,
                                                          AA_percent_sim=percentage_dic, AA_chemistry=chemistry_dic,
                                                          parameters=parameters))+'\n')

    # Create negative training examples
    files = glob.glob(phobiusfiles_from_posexamples + '/*')
    negative_seqs = []
    for file in files:
        sequences, boarders = targets_from_phobius(phobius_file=file)
        sequences, boarders = create_sequences_from_targets(sequence_list=sequences, positions_list=boarders,
                                                            avoid=positive_seqs, length=length)
        negative_seqs += sequences
    assert len(negative_seqs[0]) == len(positive_seqs[0])
    with open(training_output_prefix+".negatives",'w') as wfile:
        wfile.write(header[:-1] + '\n')
        negative_seqs = set(negative_seqs)
        for seq in negative_seqs:
            wfile.write(seq + '\t' + '\t'.join(str(_) for _ in
                                               sequence_info(sequence=seq, granthamdist=granthamDist,
                                                             AA_percent_sim=percentage_dic, AA_chemistry=chemistry_dic,
                                                             parameters=parameters)) + '\n')

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

# This can be used in an iterative fashion to predict cleavage of a number of different phobius files
# or a list of sequences can be passed to it to predict cleavage of those sequences.
def test_sequences_for_cleavage(percentage_dic=calculate_percentages(),
                                chemistry_dic=chemistry_info(),
                                granthamDist = AA_Grantham_dist(),
                                parameters = get_training_info(),
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
        seq_data.append(sequence_info(sequence=seq, granthamdist=granthamDist, AA_percent_sim=percentage_dic,
                                      AA_chemistry=chemistry_dic, parameters=parameters))
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

'''
# this function will be used to test all human genes on the cluster. DEPRECATED, DO NOT USE UNLESS MODIFIED
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
'''

if __name__ == "__main__":
    args = sys.argv[1:]
    # calculate percentages of AA at positions:
    percent_AA_normal = calculate_percentages()
    percent_AA_nonnormal = calculate_percentages(normalize=False)
    # get chemistry info for each AA:
    AA_chemistry = chemistry_info()
    # get grantham distances for each AA:
    AA_Granthamdist = AA_Grantham_dist()
    # get training info from parameter file:
    for file in glob.glob('./parameters/*'):
        parameters = get_training_info(file)
        # create positive and negative training data
        create_training_data(positive_training_example='positive_training_8AA_dedup.txt',
                             phobiusfiles_from_posexamples='./posteriors', percentage_dic=percent_AA_normal,
                             granthamDist=AA_Granthamdist, chemistry_dic=AA_chemistry,
                             training_output_prefix="training_examples", length=8, parameters=parameters)
        XOR, DATA = prepare_training_data(training_output_prefix='training_examples')




