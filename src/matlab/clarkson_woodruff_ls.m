% clarkson_woodruff.m
%
% An implementation of the Clarkson-Woodruff algorithm, based on class 
% notes and the original paper, which can be found at 
% http://arxiv.org/pdf/1207.6365v4.pdf.
%
% Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%          Erik Nelson (eanelson@eecs.berkeley.edu)

function [x] = clarkson_woodruff_ls(A, b)

% Returns an approximate solution to the problem
%
%                arg min ||Ax - b||
%                     x            2
%
% when A is a sparse mxn matrix. The algorithm is essentially expressed as
% four steps (where we define some error bound e and set
% t = (n/e)^2 log(n/e)^6):
% 1. Construct S = PD a subspace embedding matrix, where P (txm) is random
%    columns of the identity matrix and D is diagonal mxm matrix where the
%    entries are IID Bernoulli variables (either 1 or -1 with equal
%    probability).
% 2. Construct A' = SA.
% 3. Construct b' = Sb.
% 4. Solve arg min ||A'x - b'||   using randomized low rank factorization
%               x              2
%    via row extraction.

% Set constants.
[n, m] = size(A);
e = 0.98;
t = round((n/e)^2 * (log(n/e)^6));

% Construct S matrix by doing implicit matrix multiplication.
D = rand(m, 1);
D(D>0.5) = 1.0;
D(D<=0.5) = -1.0;

col_inds = randi(t, m, 1);
S = zeros(t, m);
for ii = 1:numel(col_inds)
    S(col_inds(ii), ii) = D(ii);
end

% Construct A' matrix.
A_prime = S * A;

% Construct the b' matrix.
b_prime = S * b;

% Solve the 'prime' problem with randomized low rank factorization via
% row extraction.
% TODO: HOW TO SEST RANK k??
k = ceil(0.5 * t);
p = t - k - 1;
x = randomized_low_rank_ls(A_prime, b_prime, k, p);