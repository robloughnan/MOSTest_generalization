
% Discovered Loci
n_loci_mostest_disc = array2table(zeros(10, 1));
n_loci_minp_disc = array2table(zeros(10, 1));
n_loci_mostest_valid = array2table(zeros(10, 1));
n_loci_minp_valid = array2table(zeros(10, 1));
r2_minp = array2table(zeros(10, 1));
r2_most = array2table(zeros(10, 1));
results_dir = '/path/to/results/';
for fold = 1:10
    disc_dir = sprintf('%s/fold_%i', results_dir, fold);
    assoc_stats = sprintf('%s/validation/PVS_assoc_stats.mat', disc_dir, );
    load(assoc_stats);
    n_loci_mostest_disc{fold, 1} = size(mostest_bim, 1);
    n_loci_minp_disc{fold, 1} = size(minp_bim, 1);
    r2_minp{fold, 1} = max((coefs.minp).^2);

    if strcmp(measuretype, 'area'), field = 'most_N_10'; else field = 'most_N_20'; end
    n_loci_mostest_valid{fold, 1} = sum(ps.(field)<(0.05));
    n_loci_minp_valid{fold, 1} = sum(ps.minp<(0.05));
    % n_loci_mostest_valid{fold, measuretype} = sum(ps.(field)<(0.05/length(ps.(field))));
    % n_loci_minp_valid{fold, measuretype} = sum(ps.minp<(0.05/length(ps.minp)));
    r2_most{fold, measuretype} = max((coefs.(field)).^2);
end

mean(r2_minp{:, :})
std(r2_minp{:, :})

mean(r2_most{:, :})
std(r2_most{:, :})

n_loci_tot = [mean(n_loci_minp_disc{:, :});
mean(n_loci_minp_valid{:, :});
mean(n_loci_mostest_disc{:, :});
mean(n_loci_mostest_valid{:, :})]';

sd_loci_tot = [std(n_loci_minp_disc{:, :});
std(n_loci_minp_valid{:, :});
std(n_loci_mostest_disc{:, :});
std(n_loci_mostest_valid{:, :})]';

n_loci_table = array2table(n_loci_tot, 'RowNames', meas, 'VariableNames', {'minP_disc', 'minP_valid', 'MOSTest_disc', 'MOSTest_valid'});
sd_loci_table = array2table(sd_loci_tot, 'RowNames', meas, 'VariableNames', {'minP_disc', 'minP_valid', 'MOSTest_disc', 'MOSTest_valid'});
out_fig1 = sprintf('%s/figure_1', results_dir);
mkdir(out_fig1);
writetable(n_loci_table, sprintf('%s/UKB_crossval_n_loci.tsv', out_fig1), 'FileType', 'text', 'Delimiter', '\t', 'WriteRowNames', true, 'WriteVariableNames', true);
writetable(sd_loci_table, sprintf('%s/UKB_crossval_sd_loci.tsv', out_fig1), 'FileType', 'text', 'Delimiter', '\t', 'WriteRowNames', true, 'WriteVariableNames', true);