"""
clarkson_woodruff.py

An implementation of the Clarkson-Woodruff algorithm, based on class notes and
the original paper, which can be found at http://arxiv.org/pdf/1207.6365v4.pdf.

Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
         Erik Nelson (eanelson@eecs.berkeley.edu)
"""

import numpy as np
from randomized_low_rank import RandomizedLowRankLS

def ClarksonWoodruffLS(A, b):
    """
    Returns an approximate solution to the problem

                arg min ||Ax - b||
                     x            2

    when A is a sparse mxn matrix. The algorithm is essentially expressed as
    four steps (where we define some error bound e and set
    t = (n/e)^2 log(n/e)^6):
    1. Construct S = PD a subspace embedding matrix, where P (txm) is random
       columns of the identity matrix and D is diagonal mxm matrix where the
       entries are IID Bernoulli variables (either 1 or -1 with equal
       probability).
    2. Construct A' = SA.
    3. Construct b' = Sb.
    4. Solve arg min ||A'x - b'||   using randomized low rank factorization
                  x              2
        via row extraction.
    """

    # Set constants.
    m = A.shape[0]
    n = A.shape[1]
    e = 1e-4
    t = round((n/e) * (n/e) * (np.log(n/e)**6))

    # Construct S matrix by doing implicit matrix multiplication.
    D = np.random.rand(m)
    D[D > 0.5] = 1.0
    D[D <= 0.5] = -1.0

    col_indices = np.random.randint(t, size=m)
    S = np.matlib.zeros((t, m))
    for jj, ii in enumerate(col_indices):
        S[ii, jj] = D[jj]

    # Construct A' matrix.
    A_prime = S * A

    # Construct the b' matrix.
    b_prime = S * b

    # Solve the 'prime' problem with randomized low rank factorization via
    # row extraction.
    return RandomizedLowRankLS(A_prime, b_prime)
