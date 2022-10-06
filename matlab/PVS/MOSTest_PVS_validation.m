filpth = fileparts(mfilename('fullpath'));
addpath([filpth, '/../mostest/']);
load('SurfView_surfs.mat');
icnum = 4;  icnvert = size(icsurfs{icnum}.vertices,1); % Icosahedral order to use for stats
smfnum = ''; % Need to change this to zero somoothing

% Discovery files
if ~exist('weights_file', 'var')
    if ~exist('disc_bfile', 'var'), error('bfile is required'); end
    if ~exist('disc_zstats', 'var'), error('discovery zstats file is required'); end
    if ~exist('disc_pfile', 'var'), error('discovery pfile file is required'); end
end


% Validation files
if ~exist('pheno', 'var'), error('pheno file is required'); end
if ~exist('valid_bfile', 'var'), error('Vaidation bfile is needed'); end
if ~exist('exclude_indiv', 'var'), exclude_indiv = NaN; end; 
if ~exist('include_indiv', 'var'), include_indiv = NaN; end; 
if ~exist('clusters_file', 'var'), clusters_file = NaN; end; % File containing cluster assignments for random effects of associations

if ~exist('outdir', 'var'),   error('out file prefix is required'); end

if ~exist('miss_snp', 'var'), miss_snp = 0.05; end; 
if ~exist('miss_indiv', 'var'), miss_indiv = 0.05; end; 
if ~exist('remap_vert', 'var'), remap_vert = false; end; % Maps verticies from new2orig (sometimes needed for ABCD)
if ~exist('create_plots', 'var'), create_plots = false; end;
if ~exist('save_pvs', 'var'), save_pvs = false; end;
if ~exist('force_compute_weights', 'var'), force_compute_weights=true; end % if set to false will read in pre-computed weights


tic;
% Prepare weights
save_weights_file = strrep(disc_pfile, 'pvals', 'weights');
if (~exist('weights_file', 'var') && ~isfile(save_weights_file)) || force_compute_weights
    fprintf('Preping weights\n');
    [bim, zmat_orig, C0s_inv, SVD_proj_weights] = prep_weights(disc_bfile, disc_zstats, disc_pfile, LD_ref, py2_kern);
    survive = bim{:, 'minp_survive'} | bim{:, 'most_survive'};
    bim_survive = bim(survive, :);
    save(save_weights_file, '-v7.3', 'bim_survive', 'zmat_orig', 'C0s_inv', 'SVD_proj_weights');
elseif isfile(save_weights_file)
    fprintf('Loading weights from %s\n', save_weights_file);
    load(save_weights_file);
else
    fprintf('Loading weights from %s\n', weights_file);
    load(weights_file);
end

if ~endsWith(pheno, '.mat')
    brain_dat = readtable(pheno, 'ReadVariableNames', true, 'ReadRowNames', true, 'FileType', 'text');
else
    brain_dat = load(pheno, 'ymat');
    valid_fam = readtable([valid_bfile, '.fam'], 'FileType', 'text', 'delimiter', ' ', 'ReadVariableNames', false);
    if strcmp(class(valid_fam{:, 'Var1'}), 'double')
        brain_dat = array2table(brain_dat.ymat, 'RowNames', cellstr(num2str(valid_fam{:, 1})));
    else
        brain_dat = array2table(brain_dat.ymat, 'RowNames', valid_fam{:, 1});
    end
end
% If SVD was performed phenotype for discovery
if ~isnan(SVD_proj_weights) 
    brain_dat{:, :} = brain_dat{:, :}*SVD_proj_weights; 
    % Only take top components
    Npcs = size(zmat_orig.most, 2);
    brain_dat = brain_dat(:, 1:Npcs);
end
% May need to remap verticies - different vertex ordering for freesurfer vs MMIL internal
if remap_vert
    brain_dat = remap_verticies(icsurfs, icsurfs_orig, icnvert, icnum, brain_dat, 'new2orig');
end

% Read in relevant validation SNPs (subject x SNPs)
valid_genotypes = read_genetics_chrms(valid_bfile, bim_survive, 1, 1, true, false);

% Align genetics and brain imaging
[ia, ib] = ismember(valid_genotypes.Properties.RowNames, brain_dat.Properties.RowNames);   
valid_genotypes = valid_genotypes(ia, :);
brain_dat = brain_dat(ib(ia), :);
if ~isnan(exclude_indiv)
    % Exclude individuals
    exclude = readtable(exclude_indiv, 'Delimiter', 'tab', 'ReadVariableNames', false, 'Filetype', 'text');
    if isnumeric(exclude{:, 1})
        individ_ind = ~ismember(brain_dat.Properties.RowNames, cellstr(num2str(exclude{:, 1})));
    else
        individ_ind = ~ismember(brain_dat.Properties.RowNames, exclude{:, 1});
    end
elseif ~isnan(include_indiv)
    % Include individuals
    include = readtable(include_indiv, 'Delimiter', 'tab', 'ReadVariableNames', false, 'Filetype', 'text');
    if isnumeric(include{:, 1})
        individ_ind = ismember(brain_dat.Properties.RowNames, cellstr(num2str(include{:, 1})));
    else
        individ_ind = ismember(brain_dat.Properties.RowNames, include{:, 1});
    end
end

if (~isnan(exclude_indiv)) | (~isnan(include_indiv))
    brain_dat = brain_dat(individ_ind, :);
    valid_genotypes = valid_genotypes(individ_ind, :);
end

if (~isnan(clusters_file))
    fprintf('Reading in clusters from %s\n', clusters_file);
    clusters = readtable(clusters_file, 'Delimiter', 'tab', 'ReadVariableNames', true, 'Filetype', 'text');
    if isnumeric(clusters{:, 'IID'})
        [ia, ib] = ismember(cellstr(num2str(clusters{:, 'IID'})), valid_genotypes.Properties.RowNames);
    else
        [ia, ib] = ismember(clusters{:, 'IID'}, valid_genotypes.Properties.RowNames);
    end
    assert(sum(ia)==size(valid_genotypes, 1));
    clusters = clusters{ia, 'cluster'};
end

% Change IDs so that they can be used as variables names in table
numeric_ids = startsWith(bim_survive{:, 'id'}, arrayfun(@(x) num2str(x), [0:9], 'UniformOutput', false));
bim_survive{:, 'id_clean'} = bim_survive{:, 'id'};
bim_survive{numeric_ids, 'id_clean'} = arrayfun(@(x) ['X' strrep(x{:}, ':', '_')], bim_survive{numeric_ids, 'id'}, 'UniformOutput', false);
% Min p
minp_snps = bim_survive(bim_survive{:, 'minp_survive'}, :);
[ia, ib] = ismember(valid_genotypes.Properties.VariableNames, minp_snps{:, 'id_clean'});
minp_bim = minp_snps(ib(ia), :);
% Cycle through different regularization schemes for MOSTest
mostest_snps = bim_survive(bim_survive{:, 'most_survive'}, :);

if ~isnan(SVD_proj_weights) 
    % Polygenic score to predict brain data
    pgs= struct();
    pgs.minp = double(valid_genotypes{:, ia}) * zmat_orig.minp(ib(ia), :);
    pgs.minp = array2table(pgs.minp, 'RowNames', valid_genotypes.Properties.RowNames);


    regs = fieldnames(C0s_inv);
    for reg_i = 1:numel(regs)
        reg_name = regs{reg_i};
        wvec = C0s_inv.(reg_name) * zmat_orig.most';
        [ia, ib] = ismember(valid_genotypes.Properties.VariableNames, mostest_snps{:, 'id_clean'});
        pgs.(['most_', reg_name]) = double(valid_genotypes{:, ia}) * wvec(:, ib(ia))';
        pgs.(['most_', reg_name]) = array2table(pgs.(['most_', reg_name]), 'RowNames', brain_dat.Properties.RowNames);    
        mostest_bim = mostest_snps(ib(ia), :);
    end
    % Test association with brain
    fields = fieldnames(pgs);
    ts = struct();
    coefs = struct();
    ps = struct();
    for f = 1:numel(fields)
        field = fields{f};
        [corr_coef, t, n] = nancorr_col(double(brain_dat{:, :}), pgs.(field){:, :});
        ts.(field) = t;
        coefs.(field) = corr_coef;
        ps.(field) = 2*normcdf(-abs(ts.(field)));
    end
    save(sprintf('%s/PGS_assoc_stats.mat', outdir), 'ts', 'coefs', 'ps', 'minp_bim', 'mostest_bim');
end



% Perform PVS in Validation set to get (subject x SNP matrix)
pvs = struct();
% Min p
pvs.minp = brain_dat{:, :} * zmat_orig.minp';
pvs.minp = array2table(pvs.minp, 'VariableNames', minp_snps{:, 'id_clean'},...
'RowNames', brain_dat.Properties.RowNames);
% Cycle through different regularization schemes for MOSTest
regs = fieldnames(C0s_inv);
for reg_i = 1:numel(regs)
    reg_name = regs{reg_i};
    wvec = C0s_inv.(reg_name) * zmat_orig.most';
    pvs.(['most_', reg_name]) = brain_dat{:,:} * wvec;
    pvs.(['most_', reg_name]) = array2table(pvs.(['most_', reg_name]), 'VariableNames', mostest_snps{:, 'id_clean'}, 'RowNames', brain_dat.Properties.RowNames);
end
if save_pvs
    fprintf('Saving %s/PVS.tsv', outdir);
    writetable(pvs.(['most_', reg_name]), sprintf('%s/PVS.tsv', outdir), 'Delimiter', '\t', 'FileType', 'text', 'WriteRowNames', true);
end

% Test dosage
fields = fieldnames(pvs);
ts = struct();
coefs = struct();
ps = struct();
dfs = struct();
for f = 1:numel(fields)
    field = fields{f};
    [ia, ib] = ismember(valid_genotypes.Properties.VariableNames, pvs.(field).Properties.VariableNames);
    if (isnan(clusters_file))
        if f==1, fprintf('Performing association without random effect\n'); end;
        [corr_coef, t, n] = nancorr_col(double(valid_genotypes{:, ia}), double(pvs.(field){:, ib(ia)}));
        df = n -2;
    else
        if f==1, fprintf('Performing association with random effect\n'); end;
        [corr_coef, t, n, df] = lme_col(double(valid_genotypes{:, ia}), double(pvs.(field){:, ib(ia)}), clusters);
    end
    ts.(field) = t;
    coefs.(field) = corr_coef;
    ps.(field) = 2*normcdf(-abs(ts.(field)));
    dfs.(field) = df;
    if strcmp(field, 'minp')
        minp_bim = minp_snps(ib(ia), :);
    else
        mostest_bim = mostest_snps(ib(ia), :);
    end
end
% MOSTest discovery performed in valdiation set
if exist('rediscovery_pval_file', 'var')
    pvals = load(rediscovery_pval_file);
    % Read in bim file as not saved in weights file
    bim = readtable([disc_bfile, '.bim'], 'FileType', 'Text', 'Delimiter', '\t');
    bim.Properties.VariableNames(1:6) = {'chr', 'id', 'na', 'bp', 'a1', 'a2'};
    bim{pvals.ivec_snp_good, 'resdisc_most_pval'} = 0;
    bim{pvals.ivec_snp_good, 'resdisc_most_pval'} = pvals.most_log10pval_orig;

    [ia, ib] = ismember(mostest_bim{:, 'id'}, bim{:, 'id'});
    mostest_bim{ia, 'resdisc_most_pval'} = bim{ib(ia), 'resdisc_most_pval'};
end
fprintf('Saving %s/PVS_assoc_stats.mat', outdir);
save(sprintf('%s/PVS_assoc_stats.mat', outdir), 'ts', 'coefs', 'ps', 'dfs', 'minp_bim', 'mostest_bim');

% Create bar plots of number of SNPs that validate
if create_plots
    figure(1); clf;
    validate_snps = [];
    fields = fieldnames(ps);
    for f = 1:numel(fields)
    % for f = [1, 10, 11]
        field = fields{f};
        validate_snps(f, 1) = sum(ps.(field)<(0.05/length(ps.(field))));
        FDR = mafdr(ps.(field));
        validate_snps(f, 2) = sum(FDR<0.05);
        hold on;
    end
    bar(validate_snps(:, :));
    title(measuretype, 'FontSize', 16);
    xticks(1:numel(fields));
    xticklabels(strrep(fields, '_', ' '));
    xtickangle(90);
    legend({'Bonferonni', 'FDR'}, 'Location','NorthEastOutside', 'FontSize', 16);
    ylabel('Number of Validated Loci');
    % saveas(gcf, './MOSTest_validation/matlab/myfigure.pdf');
    ax = gca;
    ax.XAxis.FontSize = 18;
    ax.YAxis.FontSize = 18;
    set(gcf,'Position',[100 100 500 300])
    saveas(gcf, sprintf('%s/number_of_loc_methods.pdf', outdir));

    % Significance in UKB vs ABCD
    figure(2); clf;
    disc_pval = mostest_bim{mostest_bim{:, 'most_survive'}, 'most_pval'};
    valid_pval = -log10(double(ps.most_N_NaN));
    valid_pval(isinf(valid_pval)) = max(disc_pval);
    scatter(disc_pval, valid_pval, 5);
    rho_pear = corr(double(disc_pval), valid_pval', 'Rows', 'complete');
    rho_spear = corr(double(disc_pval), valid_pval', 'Type', 'Spearman', 'Rows', 'complete');
    xlabel('Discovery -log_{10}(P)');
    ylabel('Validation -log_{10}(P)');
    title(sprintf('%s: \nSpearman r=%.2f \n Pearson r=%.2f', measuretype, rho_spear, rho_pear));
    saveas(gcf, sprintf('%s/Disc_Valid_scatter_sig.pdf', outdir));

    % Plot out t distribution
    figure(3); clf;
    hist(ts.most_N_NaN, 100);
    title(sprintf('%s: T-stat in Validation', measuretype));
    xlim([-max(abs(ts.most_N_NaN)), max(abs(ts.most_N_NaN))]);
    hold on;
    yl = ylim;
    plot([0, 0], [0, yl(2)],'black'); 
    saveas(gcf, sprintf('%s/t_stats.pdf', outdir));

    % Create qq plot
    figure(4); clf;
    for f = 1:2:numel(fields)
        field = fields{f};
        plot(-log10(1/length(ps.(field)):1/length(ps.(field)):1), sort(-log10(ps.(field)), 'descend')); 
        hold on;
    end
    title(measuretype);
    legend(strrep(fields(1:2:numel(fields)), '_', ' '));
    ylabel('Observed -log_{10}(P)')
    xlabel('Theoretical -log_{10}(P)')
    plot([0, 8], [0, 8], '--'); 
    saveas(gcf, sprintf('%s/qq_plot.pdf', outdir));
end
time_taken = toc();
fprintf('Finshed association, time taken %.2f\n seconds\n', time_taken);
