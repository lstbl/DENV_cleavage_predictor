# DENV Cleavage Predictor

Neural network (PyTorch) and SVM (scikit-learn) models to classify 8-mer peptide windows as cleaved or uncleaved by dengue virus protease.
The repo includes pre-trained models so you can reproduce the paper figures/tables out-of-the-box. You can also re-train from scratch by toggling a the `DoNewRuns` flag in the notebook. You can also train a completely new set of models on positive and negative training data 8-mers by toggling the `FREEZE` variable and providing the required positive and negative training data in .tsv files

What’s inside
-DengueProtease.ipynb — end-to-end workflow (load data → encode 8-mers → train/evaluate NN & SVM (or load pre-trained models) → score/rank test sites → plot).
-training_examples.positives_nopct.tsv — labeled cleaved 8-mers derived from DENV genomes
-training_examples.negatives_nopct.tsv — labeled uncleaved 8-mers derived from DENV genomes
-sites_to_test_ER_and_ER_membrane_with_encoding_nopct.tsv — unlabeled candidate 8-mers to score/rank derived from human ER proteins
-Pretrained model outputs (examples used by the notebook when DoNewRuns=False):
-nn_results_100iter_{1,2}.pkl.gz
-svm_results_100iter_{1,2}.pkl.gz
NOTE: the training examples contain amino acid biological properties that are located in the `AA_properties.txt` file

File sizes: several data/model files are ~75 MB each (below GitHub’s 100 MB limit; Git LFS not required).

# Quick start
## Suggested workflow
1) install miniconda (e.g. found at https://www.anaconda.com/docs/getting-started/miniconda/install or by using a package manager)
2) create environment with dependency using provided .yml file:
`conda env create -f environment.yml`
3) activate environment
`conda activate <env-name>`
4) Change `dire` to the working directory on your local machine

The jupyter notebook can be run using jupyterlab, ipython notebook, or third party IDE (e.g. VSCode)

### Open the notebook (using jupyter notebook)
`jupyter notebook`
Choose pretrained vs. retrain (toggle `DoNewRuns`)
To use pretrained models (default): In the first cell, ensure `DoNewRuns = False`


The notebook will load:
nninfile = './nn_results_100iter_1.pkl.gz'
svinfile = './svm_results_100iter_1.pkl.gz'

To retrain from scratch: set `DoNewRuns = True`
and run all cells. This performs multiple stochastic training runs (default Niter=100) for both models and writes fresh result files like:
```
nn_results_{Niter}iter_xxx.pkl.gz
svm_results_{Niter}iter_xxx.pkl.gz
```

to re-train a model with new data, you will need to provide a new files to replace
```
training_examples.positives_nopct.tsv
training_examples.negatives_nopct.tsv
```
This can be done by toggling the `FREEZE` variable to True and entring the appropriate filenames
 
files and a new 8mer test file (currently labeled as  `sites_to_test_ER_and_ER_membrane_with_encoding_nopct.tsv` for this study)

# Overview of models
Encoding: 8-mer sequences are one-hot encoded and added to biological properties in the input files

Models:
Neural Net (PyTorch): The neural net (NN) has three-layers. The first two layers have input and output lengths 240, followed by ReLUs. The third layer has 240 input and 2 output, followed, for training, by a sigmoid and with weighted BCE loss function.

SVM (scikit-learn): RBF kernel with C=5.0 (defaults visible in code).

Evaluation: accuracy/precision and ROC-style diagnostics are computed during training/testing--essentially perfect performance on training data for the current study (discussed in manuscript).

Scoring & Ranking: the notebook scores all candidate 8-mers (from a filtered human proteome), aggregates across runs, and reports ranks from both NN and SVM. Predefined sets of experimentally tested positives/negatives are highlighted in the final scatter/rank views.

## Runtime
Neural net runtime:
- on a consumer-grade NVIDIA GeForce RTX 3090 NN training takes ~10 sec per iteration (~16m to run 100 iterations for NN model)
- on CPU (Intel i7 2.2 GHz) takes ~180s per iteration (~5h to complete 100 iterations)

SVM runtime:
- on CPU (Intel i7 2.2 GHz) takes ~10s per iteration (~16m to complete 100 iterations)

For the current study, the results are not significantly impacted by reducing the number of iterations. This can be tried if the runtime is too burdensome for an exploratory analysis (change the value for the `Niter` parameter)

## Reproducibility tips
-Training uses stochastic elements; results and ranks will vary slightly between runs (note pre-loaded models differ slightly from the model used to produce TableS1)
-For deterministic reruns, set a seed (e.g. 
```
    seed_value = 42
    np.random.seed(seed_value)  # Seed NumPy's RNG
    random.seed(seed_value)    # Seed Python's built-in random module
    torch.manual_seed(seed_value) # Seed PyTorch's CPU RNG
```
where applicable) before training.


# Citation

If you use this code, data, or pretrained models in a publication, please cite the associated manuscript once available.

# License
free to use, modify, and distribute with attribution; provided “as is,” without warranty.