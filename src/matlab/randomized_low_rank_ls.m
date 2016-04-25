% randomized_low_rank_ls.m
%
% An implementation of randomized low rank factorization via row extraction,
% based on the description in the class notes.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)

function [x] = randomized_low_rank_ls(A, b, k)

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
%    at the top: eg. for some permutation P, P * Q = [Q1, Q2].T
% 5. The QR factorization of A is approximately Q1 * R.
% 6. Solve the least squares problem as usual, given this low-rank
%    QR decomposition.

% Constants.
% TODO: HOW TO SET p??
[m, n] = size(A);
p = round(0.25 * k);

% Y = A * F is just setting the columns of Y to be sub-sampled, phase-
% shifted FFTs of columns of A.
Y = zeros(m, k+p); 
col_inds = randi(n, k+p, 1);
phase = 2.0 * pi * randn(n, 1);
D = exp(1i * phase);

for ii = 1 : numel(col_inds)
    Y(:, ii) = fft(A(:,ii) * D(ii, ii));
end

%col_indices = np.random.randint(n, size=k + p)
%phase = 2.0 * np.pi * np.random.rand(n)
%D = np.exp(1j * phase)

%for jj, ii in enumerate(col_indices):
%    Y[:, jj] = np.fft.fft(A[:, ii] * D[ii])

% Compute QR factorization of Y.
[Q, R, P] = qr(Y);

% Get the k + p most independent rows of Q.
% TODO: WHAT DOES THIS EVEN MEAN? Q IS ORTHOGONAL!
