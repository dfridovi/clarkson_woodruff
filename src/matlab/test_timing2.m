% test_timing2.m
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
m = 1024;
n = 5;
k = n;
e = 1.0;
p = 6;
t = ceil(n^2 / e^2 * log(n / e)^6);

N_iter = 1 : 15;
% N_iter = 200; % number of calls to Clarkson-Woodruff per matrix
N_mats = 200; % number of matrices to test

density = 0.15;

cw_times = zeros(numel(N_iter), N_mats);
rlrs_times = zeros(numel(N_iter), N_mats);
cw_errors = zeros(numel(N_iter), N_mats);
rlrs_errors = zeros(numel(N_iter), N_mats);
cw_resids = zeros(numel(N_iter), N_mats);
rlrs_resids = zeros(numel(N_iter), N_mats);

for hh = 1 : numel(N_iter)
    hh
    for ii = 1 : N_mats
        
        A = sprandn(m, n, density) / sqrt(m);
        b = sprandn(m, 1, density) / sqrt(m);
        true_x = A \ b;
        true_x_norm = norm(true_x);
        
        tic
        best_resid = inf;
        best_error = inf;
        for jj = 1 : N_iter(hh)
            xhat = clarkson_woodruff_ls(A, b, t, k, p);
            this_resid = norm(A*xhat-b)/norm(A*true_x-b) - 1;
            if this_resid < best_resid
                best_resid = this_resid;
            end
            this_error = norm(true_x - xhat) / norm(true_x);
            if this_error < best_error
                best_error = this_error;
            end
        end
        cw_times(hh, ii) = toc;
        cw_resids(hh, ii) = best_resid;
        cw_errors(hh, ii) = best_error;
        
        tic
        xhat = randomized_low_rank_factorization_ls(A, b, k, p);
        rlrs_times(hh, ii) = toc;
        rlrs_resids(hh, ii) = norm(A*xhat-b)/norm(A*true_x-b) - 1;
        rlrs_errors(hh, ii) = norm(true_x-xhat) / norm(true_x);
    end
end

cw_errors_average = mean(cw_errors, 2);
rlrs_errors_average = mean(rlrs_errors, 2);
cw_times_average = mean(cw_times, 2);
rlrs_times_average = mean(rlrs_times, 2);
cw_resids_average = mean(cw_resids, 2);
rlrs_resids_average = mean(rlrs_resids, 2);

subplot(1,2,1);
hold on; box on; grid on;
plot(N_iter, cw_resids_average, '-b', 'linewidth', 2);
plot([4 4], [0 2], '--r', 'linewidth', 2);
axis([1 15 0 0.01]);
xlabel('Iterations', 'fontsize', 20, 'fontweight', 'bold');
ylabel('$$\vert\vert A\hat{x} - b \vert\vert / \vert\vert Ax - b\vert\vert - 1$$', 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold');
lh = legend('CW with \epsilon = 1', 'Iterations before more time than RLRF');
set(lh, 'fontsize', 14);

subplot(1,2,2);
hold on; box on; grid on;
plot(N_iter, cw_errors_average, '-b', 'linewidth', 2);
plot([4 4], [0 2], '--r', 'linewidth', 2);
xlabel('Iterations', 'fontsize', 20, 'fontweight', 'bold');
axis([1 15 0 2]);
ylabel('$$\vert\vert \hat{x} - x\vert\vert / \vert\vert x \vert\vert$$', 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold');
lh = legend('CW with \epsilon = 1', 'Iterations before more time than RLRF');
set(lh, 'fontsize', 14);