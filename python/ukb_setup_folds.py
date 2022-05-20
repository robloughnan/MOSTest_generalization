import argparse
import numpy as np
import pandas as pd
from os import path, makedirs

def main(id_file, out, split_start, split_end, split_frac):
    clust = pd.read_csv(id_file, sep='\t')
    if clust.columns.isin(['cluster']).sum() > 0:
        clust = clust.dropna(subset=['cluster']) # Get rid of hubby individuals
        clusters = clust.loc[:, 'cluster'].unique()
    else:
        clusters = clust.iloc[:, 0].unique()
        clust['cluster'] = clust.iloc[:, 0].values
    for i in np.arange(split_start, split_end+1):
        np.random.seed(i) # Makes sure folds should be the same each time it is called
        out_dir = path.join(out, f'''fold_{i}''')
        makedirs(out_dir, exist_ok=True)
        train_out = path.join(out_dir, 'training_sample.txt')
        test_out = path.join(out_dir, 'testing_sample.txt')
        sub_clust = np.random.choice(clusters, int(len(clust)*split_frac), replace=False)
        ind = clust.loc[:, 'cluster'].isin(sub_clust).values
        train_indiv  = clust.loc[ind, 'IID']
        test_indiv  = clust.loc[~ind, 'IID']
        pd.Series(train_indiv).to_csv(train_out, index=0, header=0)
        pd.Series(test_indiv).to_csv(test_out, index=0, header=0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Function creates cross-validation folds from id file, for conducting discovery and replication')
    parser.add_argument('-id_file', help='Filepath to IDs of individuals to be included in cross validation folds', type=str, default=None)
    # /path/to/include_indiv/UKB_king_cutoff.in.id - for maximally unrelated sample
    parser.add_argument('-out', help='Out directory to create file', type=str, default=None)
    parser.add_argument('-split_start', help='Index for start fold', type=int, default=1)
    parser.add_argument('-split_end', help='Index for last fold', type=int, default=10)
    parser.add_argument('-split_frac', help='Split fraction of full sample', type=float, default=0.67)

    opt = parser.parse_args()

    main(opt.id_file, opt.out, opt.split_start, opt.split_end, opt.split_frac)