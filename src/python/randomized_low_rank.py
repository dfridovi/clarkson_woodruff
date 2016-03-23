"""
randomized_low_rank.py

An implementation of randomized low rank factorization via row extraction,
based on the description in the class notes.

Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
         Erik Nelson (eanelson@eecs.berkeley.edu)
"""

import numpy as np

def RandomizedLowRankLS(A, b, k):
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

    # Constants.
    # QUESTION: HOW TO SET p??
    m, n = A.shape
    p = int(0.25 * k)

    # Y = A * F is just setting the columns of Y to be sub-sampled, phase-
    # shifted FFTs of columns of A.
    Y = np.matlib.zeros((m, k + p))
    col_indices = np.random.randint(n, size=k + p)
    phase = 2.0 * np.pi * np.random.rand(n)
    D = np.exp(1j * phase)

    for jj, ii in enumerate(col_indices):
        Y[:, jj] = np.fft.fft(A[:, ii] * D[ii])

    # Compute QR factorization of Y.
    Q, R = np.linalg.qr(Y)

    # Get the k + p most independent rows of Q.
    # QUESTION: WHAT DOES THIS EVEN MEAN? Q IS ORTHOGONAL!
