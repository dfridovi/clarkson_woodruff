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
m = 10;
k = 2;
p = 2;

epsilon = 0.8;
N_iter = 10; % number of calls to Clarkson-Woodruff per matrix
N_mats = 100; % number of matrices to test

t_0 = ceil((n/epsilon)^2 * log(n/epsilon)^6);
ts = ceil([t_0, 10*t_0, 200*t_0]);

figure;
for tt = [1, 2, 3]
    t = ts(tt);
    errs = zeros(N_iter * N_mats, 1);
    for ii = 1 : N_mats
        A = randn(m, n, density) / sqrt(m);
        b = randn(m, 1, density) / sqrt(m);

        % Get true solution.
        true_x = A \ b;

        for jj = 1 : N_iter
            approx_x = clarkson_woodruff_ls(A, b, t, k, p);
            %approx_x = randomized_low_rank_row_extraction_ls(A, b, k, p);
            %approx_x = randomized_low_rank_factorization_ls(A, b, k, p);
            errs(jj + (ii-1) * N_iter) = norm(approx_x - true_x) / norm(true_x);
        end
    end
    
    % Compute success percentage.
    n_correct = sum(errs < 1e-8);
    percent_correct = 100 * n_correct / length(errs);
    fprintf('Accuracy of Clarkson-Woodruff Algorithm\n m = %d, n = %d, t = %d\n',...
        m, n, t)
    fprintf('Accuracy: %3.1f%% (%d / %d)\n',...
        percent_correct, n_correct, length(errs));
    
    f = subplot(3, 1, tt);
    s = sprintf('Empirical Errors for CW for Random %d-by-%d Matrices', m, n);
    hold on;
    if (tt == 1)
        title(f, strcat(s, ', $$\varepsilon = 0.8$$'),...
             'interpreter', 'latex');
    end
    histogram(errs(errs < 2), 20, 'Normalization', 'probability');
    ylim([0 1]);
    xlim([0 2]);
    set(gca,'fontsize',14)
    xlabel('Relative error', 'fontsize', 14);
    ylabel('Normalized counts', 'fontsize', 14);
    legend(sprintf('t = %d', t));
    text(0.8, 0.8, sprintf('Percent correct: %3.1f%%', percent_correct),...
        'fontsize', 14);


end

%legend(sprintf('t = %d', t_0), sprintf('t = %d', 10*t_0), ...
%    sprintf('t = %d', 200*t_0), 'fontsize', 14);
