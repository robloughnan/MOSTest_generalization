function [genotypes, bim_tot, bim_target] = read_genetics_chrms(genetics_file, bim_target, miss_snp, miss_indiv, flip_snps, fill_missing)
    if ~exist('flip_snps', 'var'), flip_snps=true; end
    if ~exist('fill_missing', 'var'), fill_missing=true; end
    % Create another column for filter SNPs by position
    bim_target{:, 'chrbp'} = arrayfun(@(x, y) [num2str(x), ':',  num2str(y)], bim_target{:, 'chr'}, bim_target{:, 'bp'}, 'UniformOutput', false);
    fam_file = [strrep(genetics_file, '#', '1'), '.fam'];
    fam = readtable(fam_file, 'FileType', 'Text');
    nsubj = size(fam, 1);
    if contains(genetics_file, '#')
        for chrom = 1:22
            fprintf('Reading chromosome %i\n', chrom);
            bim_targ_ind = bim_target{:, 'chr'} == chrom;
            chr_bp_vec_filt = bim_target{bim_targ_ind, 'chrbp'};
            fileprefix = strrep(genetics_file, '#', num2str(chrom));
            bim = readtable([fileprefix, '.bim'], 'FileType', 'Text');
            bim.Properties.VariableNames = {'chr', 'id', 'pos', 'bp', 'a1', 'a2'};
            bim{:, 'chrbp'} = arrayfun(@(x, y) [num2str(x), ':',  num2str(y)], bim{:, 'chr'}, bim{:, 'bp'}, 'UniformOutput', false);
            % Find matching SNPS
            ia = ismember(bim{:, 'chrbp'}, bim_target{bim_targ_ind, 'chrbp'});
            bim = bim(ia, :);
            snps = find(ia);
            genomat = PlinkRead_binary2(nsubj, snps, fileprefix);
            bim_target_subset = bim_target(bim_targ_ind, :);
            [ia, ib] = ismember(bim{:, 'chrbp'}, bim_target_subset{:, 'chrbp'});
            assert(sum(ia)==length(ia));
            if flip_snps
                mismatch_snps = strcmp(bim{:, 'a2'}, bim_target_subset{ib(ia), 'a1'});
                genomat(:, mismatch_snps) = (genomat(:, mismatch_snps)*-1) + 2;
                % Reset missing SNPs to -1
                genomat(genomat==3) = -1;
            end
            % Concatonate into total matrix
            if ~exist('genotypes', 'var')
                genotypes = genomat;
                bim_tot = bim;
            else
                genotypes = [genotypes, genomat];
                bim_tot = [bim_tot; bim];
            end
        end
    % If files are not split across chromosomes
    else
        bim_tot = readtable([genetics_file, '.bim'], 'FileType', 'Text');
        bim_tot.Properties.VariableNames = {'chr', 'id', 'pos', 'bp', 'a1', 'a2'};
        bim_tot{:, 'chrbp'} = arrayfun(@(x, y) [num2str(x), ':',  num2str(y)], bim_tot{:, 'chr'}, bim_tot{:, 'bp'}, 'UniformOutput', false);
        % Find matching SNPS
        ia = ismember(bim_tot{:, 'chrbp'}, bim_target{:, 'chrbp'});
        bim_tot = bim_tot(ia, :);
        snps = find(ia);
        genotypes = PlinkRead_binary2(nsubj, snps, genetics_file);
        bim_target_subset = bim_target(:, :);
        [ia, ib] = ismember(bim_tot{:, 'chrbp'}, bim_target_subset{:, 'chrbp'});
        assert(sum(ia)==length(ia));
        if flip_snps
            mismatch_snps = strcmp(bim_tot{:, 'a2'}, bim_target_subset{ib(ia), 'a1'});
            genotypes(:, mismatch_snps) = (genotypes(:, mismatch_snps)*-1) + 2;
            % Reset missing SNPs to -1
            genotypes(genotypes==3) = -1;
        end
    end
    
    missing_entries = (genotypes==-1);
    missing_indiv_frac = sum(missing_entries)/size(genotypes, 1);
    missing_snp_frac = sum(missing_entries, 2)/size(genotypes, 2);

    % Set remainder to 0
    if fill_missing
        genotypes(genotypes==-1) = 0; 
    else
        genotypes(genotypes==-1) = NaN; 
    end

    keep_snps = missing_indiv_frac<miss_indiv;
    keep_indiv = missing_snp_frac<miss_snp;
    genotypes = genotypes(:, keep_snps);
    genotypes = genotypes(keep_indiv, :);
    bim_tot = bim_tot(keep_snps, :);

    % rename bims with bad names
    bad_names = arrayfun(@(x) contains(x, ':') | contains(x, ';'), bim_tot{:, 'id'});
    [ia, ib] = ismember(bim_tot{bad_names, 'chrbp'}, bim_target{:, 'chrbp'});
    assert (sum(ia) == sum(bad_names)); % Should be the case from ealier on

    % Clean names for putting in table
    bim_tot{bad_names, 'id'} = bim_target{ib(ia), 'id'};
    % For any remaining ':' replace with '_' for valid variable names
    numeric_ids = startsWith(bim_tot{:, 'id'}, arrayfun(@(x) num2str(x), [0:9], 'UniformOutput', false));
    bim_tot{numeric_ids, 'id'} = arrayfun(@(x) ['X' strrep(x{:}, ':', '_')], bim_tot{numeric_ids, 'id'}, 'UniformOutput', false);
    bim_target{ib(ia), 'id'} = arrayfun(@(x) ['X' strrep(x{:}, ':', '_')], bim_target{ib(ia), 'id'}, 'UniformOutput', false);
    % Still bad names
    still_bad_names = arrayfun(@(x) contains(x, ';'), bim_tot{:, 'id'});
    [ia, ib] = ismember(bim_tot{still_bad_names, 'chrbp'}, bim_target{:, 'chrbp'});
    bim_tot{still_bad_names, 'id'} = arrayfun(@(x) [strrep(x{:}, ';', '_')], bim_tot{still_bad_names, 'id'}, 'UniformOutput', false);
    bim_target{ib(ia), 'id'} = arrayfun(@(x) [strrep(x{:}, ';', '_')], bim_target{ib(ia), 'id'}, 'UniformOutput', false);
    if isnumeric(fam{:, 1})
        genotypes = array2table(genotypes, 'VariableNames', bim_tot{:, 'id'}, 'RowNames', cellstr(num2str(fam{keep_indiv, 1})));
    else
        genotypes = array2table(genotypes, 'VariableNames', bim_tot{:, 'id'}, 'RowNames', fam{keep_indiv, 1});
    end
end
