#DENV Cleavage Predictor

Neural network (PyTorchd) and SVM (scikit-learn) models to classify 8-mer peptide windows as cleaved or uncleaved by dengue virus protease.
The repo includes pre-trained models so you can reproduce the paper figures/tables out-of-the-box. You can also re-train from scratch by toggling a the "DoNewRuns" flag in the notebook.

What’s inside
-DengueProtease.ipynb — end-to-end workflow (load data → encode 8-mers → train/evaluate NN & SVM → score/rank test sites → plot).
-training_examples.positives_nopct.tsv — labeled cleaved 8-mers.
-training_examples.negatives_nopct.tsv — labeled uncleaved 8-mers.
-sites_to_test_ER_and_ER_membrane_with_encoding_nopct.tsv — unlabeled candidate 8-mers to score/rank.
-Pretrained model outputs (examples used by the notebook when DoNewRuns=False):
-nn_results_100iter_1.pkl.gz
-svm_results_100iter_1.pkl.gz
NOTE: the training examples contain amino acid biological properties that are located in the "AA_properties.txt" file

File sizes: several data/model files are ~75 MB each (below GitHub’s 100 MB limit; Git LFS not required).

#Quick start
##Set up the environment
(Optional) Create from environment.yml 
conda env create -f environment.yml
conda activate <env-name>

##Or install key deps
conda install -y python>=3.10 numpy pandas scikit-learn matplotlib
pip install torch --index-url https://download.pytorch.org/whl/cpu


###Open the notebook
jupyter notebook DengueProtease.ipynb
Choose pretrained vs. retrain (toggle "DoNewRuns")
To use pretrained models (default): In the first cell, ensure `DoNewRuns = False`


The notebook will load:
nninfile = './nn_results_100iter_1.pkl.gz'
svinfile = './svm_results_100iter_1.pkl.gz'

To retrain from scratch: set `DoNewRuns = True`
and run all cells. This performs multiple stochastic training runs (default Niter=100) for both models and writes fresh result files like:
nn_results_{Niter}iter_xxx.pkl.gz
svm_results_{Niter}iter_xxx.pkl.gz

How it works (high level)
Encoding: 8-mer sequences are one-hot encoded and added to biological properties in the input files

Models:
Neural Net (PyTorch): The neural net (NN) has three-layers. The first two layers have input and output lengths 240, followed by ReLUs. The third layer has 240 input and 2 output, followed, for training, by a sigmoid and with weighted BCE loss function.

SVM (scikit-learn): RBF kernel with C=5.0 (defaults visible in code).

Evaluation: accuracy/precision and ROC-style diagnostics are computed during training/testing--essentially perfect performance on training data.

Scoring & Ranking: the notebook scores all candidate 8-mers (from a filtered human proteome), aggregates across runs, and reports ranks from both NN and SVM. Predefined sets of experimentally tested positives/negatives are highlighted in the final scatter/rank views.


##Reproducibility tips
-Training uses stochastic elements; results and ranks will vary slightly between runs.
-For deterministic reruns, set a seed (where applicable) before training.


#Citation

If you use this code, data, or pretrained models in a publication, please cite the associated manuscript once available.

#License
All rights reserved