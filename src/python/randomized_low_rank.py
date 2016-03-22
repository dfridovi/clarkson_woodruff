"""
randomized_low_rank.py

An implementation of randomized low rank factorization via row extraction,
based on the description in the class notes.

Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
         Erik Nelson (eanelson@eecs.berkeley.edu)
"""

import numpy as np

def RandomizedLowRankLS(A, b):
    """
    Returns an approximate solution to the problem

                arg min ||Ax - b||
                     x            2

    when A is a sparse mxn matrix. The algorithm is essentially expressed as
    6 steps:
    1. Construct the matrix F = D * FFT * C, where D is a diagonal matrix
       with complex diagonal entries of unit size, FFT is the FFT basis, and
       C extracts a random subset of the columns of FFT (k + p of them, where
       k is the rank of A and p is some small positive integer).
    2. Construct Y = A * F
    3. Compute the QR factorization of Y: [Q, R] = qr(Y)
    4. Find the k + p most independent rows of Q and reorder so that they are
       at the top: eg. for some permutation P, P * Q = [Q1, Q2].T
    5. The QR factorization of A is approximately Q1 * R.
    6. Solve the least squares problem as usual, given this low-rank
       QR decomposition.
    """

    # TODO!
