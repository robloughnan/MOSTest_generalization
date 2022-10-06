% This script runs the jobs
addpath('../mostest/');

% Need to update to where own 1k data is
LD_ref='~/data/1000kGenome/phase3_EUR/chr@';
% Update to path of own python2 kernel (used in pyprune)
py2_kern = '~/anaconda3/envs/py2/bin/python';

% Perform validation on dummy data
results_dir =  '../../dummy_data'; % Discovery UKB GWAS results, also where replciation results will be saved to
pheno_dir = results_dir;
num_eigval_to_keep = 10;  % Regularization parameter for covariance matrix

% Discovery files
disc_fold_dir = results_dir;
disc_bfile = [results_dir, '/genetics_rndm'];
disc_zstats = sprintf('%s/zmat_mostest_disc_chrm%s.mat', disc_fold_dir, '%i');
disc_pfile = sprintf('%s/zmat_mostest_disc_pvals.mat', disc_fold_dir);

% Validation files
pheno = sprintf('%s/pheno.tsv', pheno_dir);
valid_bfile = disc_bfile;        
include_indiv = sprintf('%s/testing_sample.txt', disc_fold_dir);
outdir = disc_fold_dir;
if ~exist(outdir, 'dir'), mkdir(outdir); end

MOSTest_PVS_validation;

% ABCD validation
% measures={'area', 'thickness', 'sulc'};
% SVD_pheno = ''; % no SVD pheno
% pheno_dir = [results_dir, 'ABCD_phenotypes'];
% for m = 1:3
%     measuretype = measures{m};                                                                            
%     fprintf('Processing %s\n', measuretype);
%     if strcmp(measuretype, 'area')                                                                        
%         num_eigval_to_keep = 10;                                                                   
%     else                                                                                           
%         num_eigval_to_keep = 20;                                                                   
%     end 

%     % Validation files
%     pheno = sprintf('%s/ABCD_%s_ic3_sm0_resid.tsv',pheno_dir, measuretype);
%     valid_bfile = '/path/to/ABCD_bfile';

%     % Discovery files
%     disc_fold_dir = sprintf('%s/fold_%i/', results_dir, fold);
%     disc_bfile = '/path/to/UKB_bfile';
%     disc_zstats = sprintf('%s/UKB_%s_ic3_sm0_reg%i_most_disc%s_chrm%s_zmat.mat', disc_fold_dir, measuretype, num_eigval_to_keep, SVD_pheno, '%i');
%     disc_pfile = sprintf('%s/UKB_%s_ic3_sm0_reg%i_most_disc%s_pvals.mat.mat', disc_fold_dir, measuretype, num_eigval_to_keep, SVD_pheno);
%     remap_vert = false;

%     exclude_indiv = NaN;  
%     include_indiv = '/path/to/ABCD_EUR_unrel.txt';  
%     outdir = sprintf('/path/to/ABCD/replication/fold_%d/%s', fold, measuretype);                                                                                                                                                                                                               
%     if ~exist(outdir, 'dir'), mkdir(outdir); end
%     save_pvs = false;
%     MOSTest_PRS_validation;
% end
% exit()