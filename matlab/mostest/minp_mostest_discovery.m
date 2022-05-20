%% Script demonstrates how to call mostest_light and mostest_light_pvals for minp and 
%% MOSTest discovery for a single fold (fold = 1)

% This script runs the jobs
disp('Running mmatlab_job.m script...')
filpth = fileparts(mfilename('fullpath'));
addpath([filpth, '/../mostest']);


fold = 1;

bfile = '/path/to/UKB_bfile'; 
results_dir = '/path/to/results/';
pheno_dir = [results_dir, 'UKB_phenotypes'];

% Loop over measures
for m = 1:3
    measures = {'thickness', 'area', 'sulc'}; 
    meas = measures{m};      
    % Different values of regularization                                                                      
    if strcmp(meas, 'area')                                                                        
        num_eigval_to_keep = 10;                                                              
    else                                                                                           
        num_eigval_to_keep = 20;                                                                  
    end      
    pval_file = dir(sprintf('%s/fold_%d/UKB_%s_ic3_sm0_reg%d_most_disc_pvals.mat.mat', results_dir, job_id, meas, num_eigval_to_keep));                                                  
    pheno = sprintf('%s/phenotypes/UKB_%s_ic3_sm0_resid.tsv', pheno_dir, meas);

    out = sprintf('%s/fold_%d/UKB_%s_ic3_sm0_reg%d_most_disc', results_dir, job_id, meas, num_eigval_to_keep);   
    subselect_indiv = sprintf('%s/fold_%d/training_sample.txt', results_dir, job_id);

    if ~exist(fileparts(out), 'dir'), mkdir(fileparts(out)); end
    SVD_pheno=false;
    nsubj=nan;
    snps=nan;
    save_chroms=true;
    zscore_geno = true;
    mostest_light;  
    stats_file = [out, '.mat'];
    out = [out, '_pvals.mat'];
    mostest_light_pvals;  
end 
exit()