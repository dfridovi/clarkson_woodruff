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
m = 10; n = 5; r = 1;  % size of matrices
N_iter = 20; % number of calls to Clarkson-Woodruff per matrix
N_mats = 100; % number of matrices to test

t = 1000;
p = 6;

for k = 5:5:20
    errs = zeros(N_iter * N_mats, 1);
    for ii = 1 : N_mats
        A = randn(m, n) / sqrt(m);
        %A(:, r+1:n) = 0;
        b = randn(m, 1) / sqrt(m);
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
    title(sprintf('Accuracy of Clarkson-Woodruff Algorithm\n t = %d, k = %d, p = %d',...
        t, k, p));

    n_correct = sum(errs < 1e-8);
    fprintf('Accuracy of Clarkson-Woodruff Algorithm\n t = %d, k = %d, p = %d\n',...
        t, k, p)
    fprintf('Accuracy: %3.1f%% (%d / %d)\n',...
        100 * n_correct / length(errs), n_correct, length(errs));
end