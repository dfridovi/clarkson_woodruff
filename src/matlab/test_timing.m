% test_accuracy.m
%
% Generate relatively small random matrices and test the accuracy of
% Clarkson-Woodruff against normal A \ b.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)
close all; clear; clc;

warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:singularMatrix');
warning('on', 'MATLAB:nearlySingularMatrix');

% Parameters.
m = 128;
n = 5;
k = n;
e = 1;
p = 6;
t = ceil(n^2 / e^2 * log(n / e)^6);

N_iter = 200; % number of calls to Clarkson-Woodruff per matrix
N_mats = 200; % number of matrices to test

density = 0.01:0.02:0.37+eps;

cw_times = zeros(numel(density), N_iter * N_mats);
rlrs_times = zeros(numel(density), N_iter * N_mats);

for d = 1 : numel(density)
    d
    for ii = 1 : N_mats
        ii
        A = sprandn(m, n, density(d)) / sqrt(m);
        b = sprandn(m, 1, density(d)) / sqrt(m);
        true_x = A \ b;
        
        for jj = 1 : N_iter
            tic
            clarkson_woodruff_ls(A, b, t, k, p);
            cw_times(d, jj + (ii-1) * N_iter) = toc;
            
            tic
            randomized_low_rank_factorization_ls(A, b, k, p);
            rlrs_times(d, jj + (ii-1) * N_iter) = toc;
        end
    end
end