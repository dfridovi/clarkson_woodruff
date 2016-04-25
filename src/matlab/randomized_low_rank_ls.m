% randomized_low_rank_ls.m
%
% An implementation of randomized low rank factorization via row extraction,
% based on the description in the class notes.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)

function [x] = randomized_low_rank_ls(A, b, k, p)

% Returns an approximate solution to the problem
%
%           arg min ||Ax - b||
%                x            2
%
% when A is a sparse mxn matrix. The algorithm is essentially expressed as
% 6 steps:
% 1. Construct the matrix F = D * FFT * C, where D is a diagonal matrix
%    with complex diagonal entries of unit size, FFT is the FFT basis, and
%    C extracts a random subset of the columns of FFT (k + p of them, where
%    k is the rank of A and p is some small positive integer).
% 2. Construct Y = A * F
% 3. Compute the QR factorization of Y: [Q, R] = qr(Y)
% 4. Find the k + p most independent rows of Q and reorder so that they are
%    at the top: eg. for some permutation P, P * Q = [Q1; Q2].
% 5. Let X = P * Q * inv(Q1) = [I; Q2 * inv(Q1)].
% 6. Let P * A = [A1; A2], where A1 has k + p rows, and A ~ P' * X * A1.
% 7. Solve the least squares problem as usual, given this low-rank
%    QR decomposition.

% Constants.
% TODO: HOW TO SET p??
[m, n] = size(A);

% Y = A * F is just setting the columns of Y to be sub-sampled, phase-
% shifted FFTs of columns of A.
Y = zeros(m, k+p); 
col_inds = randi(n, k+p, 1);
phase = 2.0 * pi * rand(n, 1);
D = exp(1i * phase);

for ii = 1 : numel(col_inds)
    Y(:, ii) = fft(A(:,col_inds(ii)) * D(col_inds(ii)));
end

% Compute QR factorization of Y.
[Q, ~] = qr(Y);

% Get the k + p most independent rows of Q... Estimate as the first k + p
% indices in the permutation vector returned by a QR decomposition.
[~, ~, perm] = qr(Q, 'vector');
perm_t(perm) = 1:length(perm);
Q1 = Q(perm(1:k + p), :);
Q2 = Q(perm(k + p + 1:length(perm)), :);

% Create X matrix.
X = [eye(k + p); Q2 * Q1'];

% Create low-rank approximation to A.
A1 = A(perm(1:k + p), :);
temp = X * A1;
A_approx = temp(perm_t(:), :);

% Solve.
x = A_approx \ b;