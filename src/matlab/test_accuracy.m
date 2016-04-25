% test_accuracy.m
%
% Generate relatively small random matrices and test the accuracy of
% Clarkson-Woodruff against normal A \ b.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)

warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:singularMatrix');

% Parameters.
n = 10;  % size of matrices
N_iter = 20; % number of calls to Clarkson-Woodruff per matrix
N_mats = 100; % number of matrices to test

t = 20;
k = 5;
p = 10;

errs = zeros(N_iter * N_mats, 1);
for ii = 1 : N_mats
    A = randn(n) / sqrt(n);
    b = randn(n, 1) / sqrt(n);
    true_x = A \ b;
    
    for jj = 1 : N_iter
        approx_x = clarkson_woodruff_ls(A, b, t, k, p);
        errs(jj + (ii-1) * N_mats) = norm(approx_x - true_x) / norm(true_x);
    end
end

figure;
histogram(errs, 20, 'Normalization', 'probability');
xlabel('Log of error');
ylabel('Normalized counts');
title(sprintf('Accuracy of Clarkson-Woodruff Algorithm\n t = %d, k = %d, p = %d',...
    t, k, p));

n_correct = sum(errs < 1e-8);
fprintf('Accuracy of Clarkson-Woodruff Algorithm\n t = %d, k = %d, p = %d\n',...
    t, k, p)
fprintf('Accuracy: %3.1f%% (%d / %d)\n',...
    100 * n_correct / length(errs), n_correct, length(errs));