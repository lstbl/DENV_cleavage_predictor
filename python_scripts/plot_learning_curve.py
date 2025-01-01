import numpy as np
import matplotlib.pyplot as plt
from sklearn.learning_curve import learning_curve
from sklearn import svm



def plot_learning_curve(estimator, title, DATA, DATA_LABELS, ylim=None, cv=None,
                        n_jobs=4, train_sizes=np.linspace(.05, 1.0, 10)):
    plt.figure()
    plt.title(title)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, DATA, DATA_LABELS, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    plt.legend(loc="best")
    return plt

def plot_SenSpec_curve(only_positives = [], only_negatives = [], model = svm.SVC(kernel='rbf', C=10),title="rbf Kernel",iterations=100):
    TPR_means = []
    TPR_stdev = []
    SPC_means = []
    SPC_stdev = []
    PPV_means = []
    PPV_stdev = []
    F1_means = []
    F1_stdev = []
    percent_sampled = []
    for N in np.linspace(0.05,0.8,10):
        TPR_results = []
        SPC_results = []
        PPV_results = []
        F1_results = []
        for sampleN in range(iterations):
            pos_samples = set(np.random.choice(range(len(only_positives)),size=int(len(only_positives)*N),replace=False))
            pos_test = set(range(len(only_positives))) - pos_samples
            positives, positives_to_test = np.array([only_positives[_] for _ in pos_samples]), \
                                           np.array([only_positives[_] for _ in pos_test])
            neg_samples = set(np.random.choice(range(len(only_negatives)),
                                               size=int(len(only_negatives)*N),replace=False))
            neg_test = set(range(len(only_negatives))) - neg_samples
            negatives, negatives_to_test = np.array([only_negatives[_] for _ in neg_samples]), \
                                           np.array([only_negatives[_] for _ in neg_test])
            XOR_positives = np.array([True for _ in positives])
            XOR_negatives = np.array([False for _ in negatives])
            training_data = np.concatenate((positives,negatives))
            training_XOR = np.append(XOR_positives,XOR_negatives)
            training_model = model
            training_model.fit(training_data,training_XOR)
            TPR_results.append(sum(training_model.predict(positives_to_test))/float(len(positives_to_test)))
            SPC_results.append((len(negatives_to_test)-sum(training_model.predict(negatives_to_test)))/float(len(negatives_to_test)))
            PPV_results.append(sum(training_model.predict(positives_to_test))/float(sum(training_model.predict(negatives_to_test))+sum(training_model.predict(positives_to_test))))
            F1_results.append((2*TPR_results[-1]*PPV_results[-1])/(TPR_results[-1]+PPV_results[-1]))
        percent_sampled.append(N)
        TPR_means.append(np.mean(TPR_results))
        TPR_stdev.append(np.std(TPR_results))
        SPC_means.append(np.mean(SPC_results))
        SPC_stdev.append(np.std(SPC_results))
        PPV_means.append(np.mean(PPV_results))
        PPV_stdev.append(np.std(PPV_results))
        F1_means.append(np.mean(F1_results))
        F1_stdev.append(np.std(F1_results))
    TPR_means, TPR_stdev = np.array(TPR_means), np.array(TPR_stdev)
    SPC_means, SPC_stdev = np.array(SPC_means), np.array(SPC_stdev)
    PPV_means, PPV_stdev = np.array(PPV_means), np.array(PPV_stdev)
    F1_means, F1_stdev = np.array(F1_means), np.array(F1_stdev)
    plt.figure()
    plt.title(title)
    plt.xlabel("percent sampled from positive and negative training samples")
    plt.ylabel("Sensitivity/Specificity")
    plt.grid()
    plt.fill_between(percent_sampled,TPR_means-TPR_stdev,TPR_means+TPR_stdev,alpha=0.1,color='r')
    plt.plot(percent_sampled,TPR_means,'o-',color='r',label="Sensitivity")
    plt.fill_between(percent_sampled,SPC_means-SPC_stdev,SPC_means+SPC_stdev,alpha=0.1,color='g')
    plt.plot(percent_sampled,SPC_means,'o-',color='g',label="Specificity")
    plt.fill_between(percent_sampled,PPV_means-PPV_stdev,PPV_means+PPV_stdev,alpha=0.1,color='b')
    plt.plot(percent_sampled,PPV_means,'o-',color='b',label="Precision")
    plt.fill_between(percent_sampled, F1_means - F1_stdev, F1_means + F1_stdev, alpha=0.1, color='k')
    plt.plot(percent_sampled, F1_means, 'o-', color='k', label="F1-score")
    plt.legend(loc="best")
    return plt

