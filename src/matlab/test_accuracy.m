% test_accuracy.m
%
% Generate relatively small random matrices and test the accuracy of
% Clarkson-Woodruff against normal A \ b.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)
close all;

warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:singularMatrix');
warning('on', 'MATLAB:nearlySingularMatrix');

% Parameters.
n = 2;
k = 5;
p = 3;
density = 0.9;
epsilon = 0.7;
N_iter = 10; % number of calls to Clarkson-Woodruff per matrix
N_mats = 100; % number of matrices to test

t = round((n/epsilon)^2 * log(n/epsilon)^6);

for m = round([0.1*t, 0.5*t])
    errs = zeros(N_iter * N_mats, 1);
    for ii = 1 : N_mats
        A = sprandn(m, n, density) / sqrt(m);
        b = sprandn(m, 1, density) / sqrt(m);

        % Get true solution.
        true_x = A \ b;

        for jj = 1 : N_iter
            approx_x = clarkson_woodruff_ls(A, b, t, k, p);
            %approx_x = randomized_low_rank_row_extraction_ls(A, b, k, p);
            %approx_x = randomized_low_rank_factorization_ls(A, b, k, p);
            errs(jj + (ii-1) * N_iter) = norm(approx_x - true_x) / norm(true_x);
        end
    end

    figure;
    histogram(errs, 20, 'Normalization', 'probability');
    xlabel('Relative error');
    ylabel('Normalized counts');
    title(sprintf('Accuracy of Clarkson-Woodruff Algorithm\n m = %d, n = %d',...
        m, n));

    n_correct = sum(errs < 1e-8);
    fprintf('Accuracy of Clarkson-Woodruff Algorithm\n m = %d, n = %d\n',...
        m, n)
    fprintf('Accuracy: %3.1f%% (%d / %d)\n',...
        100 * n_correct / length(errs), n_correct, length(errs));
end