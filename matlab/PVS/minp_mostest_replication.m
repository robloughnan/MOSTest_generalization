% This script runs the jobs
disp('Running mmatlab_job.m script...')
addpath('../mostest/');

fold = 1;
% UKB validation
measures={'area', 'thickness', 'sulc'};
results_dir =  '/path/to/results/'; % Discovery UKB GWAS results, also where replciation results will be saved to
pheno_dir = [results_dir, 'UKB_phenotypes'];

SVD_pheno = ''; % no SVD pheno
for m = 1:3
    measuretype = measures{m};     
    fprintf('Processing %s:\n',measuretype);                                                                      
    if strcmp(measuretype, 'area')                                                                        
        num_eigval_to_keep = 10;                                                                   
    else                                                                                           
        num_eigval_to_keep = 20;                                                                   
    end 
    fprintf('\tProcessing fold: %d\n', fold);
    % Discovery files
    disc_fold_dir = sprintf('%s/fold_%i/', results_dir, fold);
    disc_bfile = '/path/to/UKB_bfile';
    disc_zstats = sprintf('%s/UKB_%s_ic3_sm0_reg%i_most_disc%s_chrm%s_zmat.mat', disc_fold_dir, measuretype, num_eigval_to_keep, SVD_pheno, '%i');
    disc_pfile = sprintf('%s/UKB_%s_ic3_sm0_reg%i_most_disc%s_pvals.mat.mat', disc_fold_dir, measuretype, num_eigval_to_keep, SVD_pheno);

    % Validation files
    pheno = sprintf('%s/UKB_%s_ic3_sm0_resid.tsv', pheno_dir, measuretype);
    % pheno = sprintf('/space/syn48/1/data/UKB/MOSTest/ymat_%s_ic3_sm0.mat', measuretype);
    valid_bfile = disc_bfile;        
    include_indiv = sprintf('%s/testing_sample.txt', disc_fold_dir);
    outdir = sprintf('%s/%s/', disc_fold_dir, measuretype);
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    MOSTest_PRS_validation;
end

% ABCD validation
pheno_dir = [results_dir, 'ABCD_phenotypes'];
for m = 1:3
    measuretype = measures{m};                                                                            
    fprintf('Processing %s\n', measuretype);
    if strcmp(measuretype, 'area')                                                                        
        num_eigval_to_keep = 10;                                                                   
    else                                                                                           
        num_eigval_to_keep = 20;                                                                   
    end 

    % Validation files
    pheno = sprintf('%s/ABCD_%s_ic3_sm0_resid.tsv',pheno_dir, measuretype);
    valid_bfile = '/path/to/ABCD_bfile';

    % Discovery files
    disc_fold_dir = sprintf('%s/fold_%i/', results_dir, fold);
    disc_bfile = '/path/to/UKB_bfile';
    disc_zstats = sprintf('%s/UKB_%s_ic3_sm0_reg%i_most_disc%s_chrm%s_zmat.mat', disc_fold_dir, measuretype, num_eigval_to_keep, SVD_pheno, '%i');
    disc_pfile = sprintf('%s/UKB_%s_ic3_sm0_reg%i_most_disc%s_pvals.mat.mat', disc_fold_dir, measuretype, num_eigval_to_keep, SVD_pheno);
    remap_vert = false;

    exclude_indiv = NaN;  
    include_indiv = '/path/to/ABCD_EUR_unrel.txt';  
    outdir = sprintf('/path/to/ABCD/replication/fold_%d/%s', fold, measuretype);                                                                                                                                                                                                               
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    save_pvs = false;
    MOSTest_PRS_validation;
end
exit()