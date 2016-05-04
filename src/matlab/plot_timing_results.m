% Plot timing results.
close all; clear; clc
load cw_timing2;

cw_times_average = mean(cw_times, 2);
rlrs_times_average = mean(rlrs_times, 2);
cw_times_average(20:end) = [];
rlrs_times_average(20:end) = [];
density(20:end) = [];

figure; hold on; grid on; box on;

% Plot vertical line where O(nnz(A)) dominates O(n*t);
% x = [n*t n*t]; 
% y = [0 1];
% plot(x, y, '-k', 'linewidth', 2);

plot(density*n*m, cw_times_average, '-r', 'linewidth', 2);
plot(density*n*m, rlrs_times_average, '--b', 'linewidth', 2);
xlabel('nnz($$A$$)', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Time (s)', 'fontsize', 20);
lh = legend('CW($$A$$, $$b$$)', 'RLRF($$A$$, $$b$$)');
set(lh, 'interpreter', 'latex', 'fontsize', 18);

title('Timing comparison for $$A \in R^{1024 x 5}$$', 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold');
