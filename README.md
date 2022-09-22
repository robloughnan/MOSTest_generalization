# Introduction

MOSTest is a tool for join genetic analysis of multiple traits, using multivariate analysis to boost the power of discovering associated loci. Once discovery is performed within one dataset (e.g. UK Biobank) for vertex level imaging data the scripts in this repository can be used to validate within another dataset (e.g. ABCD).

# Setup Training and Testing Folds 

This step is only required if cross validation is performed within a single dataset.

Run `./python/ukb_setup_folds.py` as:
``python ./python/ukb_setup_folds.py -id_file /path/to/clust/UKB_king_cutoff.in.id -out /path/to/results`` 

Where `/path/to/clust/UKB_king_cutoff.in.id` is a file listing the IDs to be included analysis, here a maximally unrelated set of individuals created using `plink –king-cutoff 0.0625`. This will create 10 folders in `/path/to/results` of the format `fold_#`, each containing two files `training.txt` and `testing.txt` listing individuals for training and testing of min-P and MOSTest.

# min-P and MOSTest Discovery

This step performs discovery, for a single fold, on pre-residualized multivariate phenotypes contained in a `.tsv` file with row names a IIDs in associated plink files and a header indicating the phenotype name. Update filepaths of `bfile`, `pheno` and `out_dir` in `minp_mostest_discovery.m` and then run:
``matlab -nojvm -nodisplay -r minp_mostest_discovery``

This step can take very long, if performed on a GPU enabled machine computations will be performed using GPU. Additionally, if performing this step across multiple training folds it may be suitable to split the folds across a computing cluster. 

# min-P and PVS/MOSTest Replication

`minp_mostest_replication.m` demonstrates how to call `MOSTest_PVS_validation.m` for a single fold. `MOSTest_PVS_validation.m` takes outputs of `mostest_light.m` (`disc_bfile.m`, `disc_zstats.m` and `disc_pfile.m`) as inputs. This function calls `pre_weights.m` to select lead SNPs and extract z statistics, as well as whitening matricies `C0s_inv` to decorrelate z stats for regularization. After these weights have been produced, seperate PolyVertex Scores (`pvs`) are calculated for each discovered loci and at different regularizations. These pvs are then correlated with each respective loci to perform replication. 
