import os, sys, click, numpy as np, pandas as pd
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.model_selection import GridSearchCV, ShuffleSplit, KFold, cross_val_score
try :
    from .configure import readFasta, uopen
except :
    from configure import readFasta, uopen

@click.command()
@click.option('-g', '--gtog', help='outfile', required=True)
def main(gtog) :
    param_vec = []
    group_vec = []
    group_tag = {}
    with uopen(gtog, 'rt') as fin :
        for line in fin :
            if line.startswith('@') :
                part = line.strip().split()
                group_tag[part[1]] = len(group_tag)
            else :
                part = line.strip().split()
                grp = group_tag[part[3]]
                param = np.array(part[5:]).astype(np.float32)
                param_vec.append(param)
                group_vec.append(grp)
    group_vec = np.array(group_vec, dtype=int)
    param_vec = np.array(param_vec, dtype=np.float32)
    non_nested_scores = np.zeros(10)
    nested_scores = np.zeros(10)
    for i in np.arange(10) :
        inner_cv = ShuffleSplit(n_splits=10, train_size=0.2, test_size=0.2, random_state=i)
        outer_cv = ShuffleSplit(n_splits=10, train_size=0.2, test_size=0.2, random_state=i+1)
        #inner_cv = KFold(n_splits=4, shuffle=True, random_state=i)
        #outer_cv = KFold(n_splits=4, shuffle=True, random_state=i+1)
        #param_grid = {'C':[1, 10, 100, 1000]}
        #clf = GridSearchCV(LinearSVC(dual=False, max_iter=10000), param_grid=param_grid, cv=inner_cv, n_jobs=-1, verbose=1)

        param_grid = {'average': [True, False],
                      #'l1_ratio': np.linspace(0, 1, num=3),
                      'alpha': np.power(10, np.arange(-4, 1, dtype=float))}
        est = SGDClassifier(loss='log', penalty='l2', fit_intercept=True)
        clf = GridSearchCV(est, param_grid=param_grid, cv=inner_cv, n_jobs=-1, verbose=1)

        clf.fit(param_vec, group_vec)
        non_nested_scores[i] = clf.best_score_
        nested_score = cross_val_score(clf, X=param_vec, y=group_vec, cv=outer_cv, n_jobs=-1, verbose=1)
        nested_scores[i] = nested_score.mean()
        print(group_vec)

if __name__ == '__main__' :
    main()
