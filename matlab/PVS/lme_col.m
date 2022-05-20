function [corr_coef, t, n, df] = lme_col(genotypes, pvs, grouping)
    n_snps = size(genotypes, 2);
    corr_coef = nan(1, n_snps);
    t = nan(1, n_snps);
    n = nan(1, n_snps);
    df = nan(1, n_snps);
    X0 = ones(size(genotypes, 1), 1);
    Z = ones(size(genotypes, 1), 1);
    for i = 1:n_snps
        X = [X0, genotypes(:, i)];
        y = pvs(:, i);
        if length(unique(X(:, 2)))~=1
            lme = fitlmematrix(X,y,Z,grouping);
            corr_coef(i) = lme.Coefficients{2, 'Estimate'};
            t(i) = lme.Coefficients{2, 'tStat'};
            n(i) = lme.NumObservations;
            df(i) = lme.Coefficients{2, 'DF'};
        end
    end
end