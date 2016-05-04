% Demonstrates how absurdly impractical this algorithm is by displaying the
% absolute minimum number of rows necessary for a given number of columns
% before the CW algorithm has the slimmest hope of ever beating RLRS on
% time.

make_image = false;

n = 2.^(1.0:.1:18);
epsilons = 0.1 : 0.3 : 1; % Worst possible algorithmic accuracy.

% Log log plot.
figure;
for e = 1 : numel(epsilons)
    y = (n./epsilons(e)).^2 .* log(n./epsilons(e)).^6;
    loglog(n, y, 'linewidth', 2); 
    hold on;
end

% Also plot the maximum number of rows that will fit in 16 Gb RAM for each
% number of columns.
m = 2e9 ./ n;
loglog(n, m, 'linewidth', 2, 'linestyle', '-', 'color', 'k');

% Amount of data able to fit in largest data center in the world.
m = (12e18 / 8) ./ n;
loglog(n, m, 'linewidth', 2, 'linestyle', '-.', 'color', 'k');

% Amount of data in the world.
m = (2.9e20 / 8) ./ n;
loglog(n, m, 'linewidth', 2, 'linestyle', ':', 'color', 'k');

grid on; box on;
title('Minimum feasible number of rows for Clarkson-Woodruff');
xlabel('Columns, $$n$$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Rows, $$m$$', 'interpreter', 'latex', 'fontsize', 20);
lh = legend('\epsilon = 0.1',...
            '\epsilon = 0.4',...
            '\epsilon = 0.7',...
            '\epsilon = 1.0',...
            '16 Gb RAM',...
            '12 Eb World`s Largest Data Center',...
            '0.2 Zb Memory on Earth in 2011');
set(lh, 'fontsize', 16);

if make_image
    addpath ~/Documents/MATLAB/export_fig;
    set(gcf, 'color', 'none');
    export_fig -pdf 'MinimumFeasibleRows';
    set(gcf, 'color', 'white');
end