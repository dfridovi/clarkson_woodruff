% randomized_low_rank_factorization_ls.m
%
% An implementation of randomized low rank factorization,
% based on the description in the class notes.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)

function [x] = randomized_low_rank_factorization_ls(A, b, k, p)

% Constants.
[m, n] = size(A);

% Y = A * F. F is just a bunch of random Gaussians.
F = randn(n, k + p) / sqrt(m);
Y = A * F;

% Compute QR factorization of Y.
[Q, ~] = qr(Y);

% Form B = Q' * A so that A ~ Q * B = Q * Q' * A
B = Q' * A;

% Solve.
x = B \ (Q' * b);