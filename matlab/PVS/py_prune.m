function [bim_survive] = py_prune(bim, pval_file, stat, LD_ref, py2_kern)
    % bim_survive is a vector the same length as bim with Trues where it survives pruning from plink
    sumstat_py = '../python/sumstats.py';
    % LD_ref = '/space/syn03/1/data/cfan/1000Genome/phase3/build37_released/plink_binary/EUR.chr@.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_clean';
    
    % Renmae columns for saving out in standardized form
    bim.Properties.VariableNames(1:6) = {'CHR', 'SNP', 'na', 'BP', 'A1', 'A2'};

    bim{:, 'PVAL'} = 10.^-bim.([stat, '_pval']);

    % Drop NAN
    bim_out = bim(~isnan(bim{:, 'PVAL'}), :);

    % Save out sumstats format
    sumstat_out = strrep(pval_file, '.mat', ['_', stat, '.sumstats']);
    save_cols = {'SNP', 'CHR', 'BP', 'PVAL', 'A1', 'A2', 'N'};
    writetable(bim_out(:, save_cols), sumstat_out, 'FileType', 'text', 'Delimiter', '\t');

    clump_out = [sumstat_out, '.clump'];
    py_clump_cmd = sprintf('%s %s clump --clump-field PVAL --plink plink --sumstats %s --bfile-chr %s --exclude-ranges 6:25119106-33854733 --clump-p1 5E-8 --out %s',...
     py2_kern, sumstat_py, sumstat_out, LD_ref, clump_out);
    
    status = system(py_clump_cmd);

    % Read in ressults
    clumped = readtable([clump_out, '.loci.csv'], 'FileType', 'text', 'Delimiter', '\t');
    bim_survive = ismember(bim{:, 'SNP'}, clumped{:, 'LEAD_SNP'});
end

