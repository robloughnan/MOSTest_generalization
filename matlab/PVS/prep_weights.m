function [bim, zmat_orig, C0s_inv, SVD_proj_weights] = prep_weights(bfile, z_stats_file, pval_file, bim_target)
    filpth = fileparts(mfilename('fullpath'));
    addpath([filpth '/../mostest']);
    r2_thresh=0.2;
    

    bim = readtable([bfile, '.bim'], 'FileType', 'Text', 'Delimiter', '\t');
    bim.Properties.VariableNames(1:6) = {'chr', 'id', 'na', 'bp', 'a1', 'a2'};

    
    pvals = load(pval_file);
    % Add N column
    bim{:, 'N'} = NaN;
    bim{pvals.ivec_snp_good, 'N'} = pvals.nvec;
    % Prune p values
    for stat = {'most', 'minp'}
        fprintf('Pruning %s snps\n', stat{:});
        pval_vec = zeros(length(pvals.ivec_snp_good), 1);
        pval_vec(pvals.ivec_snp_good) = -pvals.([stat{:}, '_log10pval_orig']);
        bim{:, [stat{:}, '_pval']} = -1*(pval_vec);
        bim{~pvals.ivec_snp_good, [stat{:}, '_pval']} = NaN;
        if ~exist('bim_target', 'var')
            bim{:, [stat{:}, '_survive']} = py_prune(bim, pval_file, stat{:}); % plink prune
            % bim{:, [stat{:}, '_survive']} = prune(bfile, pval_vec, log10(5E-8), r2_thresh);
        else
            % Set bim survive to target - to enable those snps to be saved
            bim_target_surv  = bim_target(bim_target{:, [stat{:}, '_survive']}, :);
            ia = ismember(bim{:, 'id'}, bim_target_surv{:, 'id'});
            bim{:, [stat{:}, '_survive']} = ia;
        end
    end

    % Save z-stats for most significant hits
    % Read in from zstats from smallest chromosome
    zstats = load(sprintf(z_stats_file, max(bim{:, 'chr'}))); 
    if ismember('keep', fieldnames(zstats))
        keep = zstats.keep;
    else
        keep = logical(ones(1, size(zstats.zmat_orig, 2)));
    end

    if ismember('SVD_proj_weights', fieldnames(zstats)) 
        if ~isnan(zstats.SVD_proj_weights)
            SVD_proj_weights = zeros(length(keep));
            SVD_proj_weights(keep, keep) = zstats.SVD_proj_weights;
            keep = logical(ones(1, size(zstats.zmat_orig, 2))); % Keep now refers to the input vector for saving zmat_orig
        else
            SVD_proj_weights = NaN;
        end
    else
        SVD_proj_weights = NaN;
    end

    zmat_orig = struct('most', single(zeros(sum(bim{:, 'most_survive'}), length(keep))),...
    'minp', single(zeros(sum(bim{:, 'minp_survive'}), length(keep))));

    i = struct('minp', 1, 'most', 1);
    fprintf('Extracting significant zstats\n');
    for chr = [unique(bim{:, 'chr'})]'
        fprintf('Processing chromosome %i\n', chr);
        fprintf('Loading %s ...', sprintf(z_stats_file, chr));
        chr_zmat = load(sprintf(z_stats_file, chr), 'zmat_orig');
        fprintf('\t Complete\n')
        chr_bim = bim(bim{:, 'chr'} == chr, :);
        for stat = {'most', 'minp'}
            snp_ind = logical(chr_bim{:, [stat{:}, '_survive']});
            zmat_ind = [i.(stat{:}): i.(stat{:}) + sum(snp_ind) - 1];
            zmat_orig.(stat{:})(zmat_ind, keep) = chr_zmat.zmat_orig(snp_ind, :);
            if stat{:} =='minp'
                % Only preseverve most significant snp
                z_orig_tmp = abs(zmat_orig.(stat{:})(zmat_ind, :));
                [~, max_ind] = max(z_orig_tmp, [], 2);
                mask = logical(zeros(length(max_ind), size(z_orig_tmp, 2))); % Mask to preserve only top zstat
                for j = 1:length(max_ind), mask(j, max_ind(j)) = true; end
                z_orig_tmp(~mask) = 0;
                zmat_orig.(stat{:})(zmat_ind, :) = z_orig_tmp;
            end
            i.(stat{:}) = i.(stat{:}) + sum(snp_ind);
        end
    end

    % Create weightings
    if ismember('ymat_corr', fieldnames(zstats))
        C0 = zstats.ymat_corr; % Compute overall phenotypic covariance matrix -- should residualize phenotypes first, then inverse rank normalize, then normalize variance (if needed)
        kmaxlist = [1 5 10 20 50 100 200 500 inf, nan]; % Nan is idenety matrix (no regularizing)
        [U S] = svd(C0); s = diag(S);
        C0s_reg = struct(); C0s_inv = struct;
        for kmaxi = 1:length(kmaxlist)
            kmax = kmaxlist(kmaxi);
            if ~isnan(kmax)
                s_reg = max(s(min(kmax,length(s))),s);
                % Set verticies outside of 'keep' to zero
                C0_reg = zeros(size(zmat_orig.most, 2));
                C0_inv = zeros(size(zmat_orig.most, 2));
                C0_reg(keep, keep) = U*diag(s_reg)*U';
                C0_inv(keep, keep) = U*diag(s_reg.^-1)*U';
                C0s_reg.(['N_', num2str(kmax)]) = C0_reg;
                C0s_inv.(['N_', num2str(kmax)]) = C0_inv;
            else
                C0s_inv.(['N_', num2str(kmax)]) = eye(size(C0_inv));
            end
        end
    else
        C0s_inv = [];
    end
end

% Get keep index for Brain imaging
% ymat = readtable(pheno_file, 'FileType', 'Text', 'Delimiter', '\t');
% keep = (min(ymat{:, :})~=max(ymat{:, :}));
% assoc = table(bim{:, 1}, bim{:, 2}, bim{:, 4}, bim{:, 5}, bim{:, 6},...
% zeros(size(bim, 1), 1), zeros(size(bim, 1), 1), zeros(size(bim, 1), 1), pval_vec,...
% 'VariableNames', {'CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'BETA', 'STAT', 'P'});
