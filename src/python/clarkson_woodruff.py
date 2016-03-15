"""
clarkson_woodruff.py

An implementation of the Clarkson-Woodruff algorithm, based on class notes and the original
paper, which can be found at http://arxiv.org/pdf/1207.6365v4.pdf.

Authors: David Fridovich-Keil (dfk@eecs.berkeley.edu)
         Erik Nelson (eanelson@eecs.berkeley.edu)
"""

def ClarksonWoodruffLS(A, b):
    """
    Returns an approximate solution to the problem

                arg min ||Ax - b||
                     x            2

    when A is a sparse mxn matrix. The algorithm is essentially expressed as four steps (where
    we define some error bound e and set t = (n/e)^2 log(1/e)^6):
    1. Construct S = PD a subspace embedding matrix, where P is t random columns of the
       identity matrix and D is diagonal mxm matrix where the entries are IID Bernoulli
       variables (either 1 or -1 with equal probability).
    2. Construct A' = SA.
    3. Construct b' = Sb.
    4. Solve arg min ||A'x - b'||   using randomized low rank factorization via row extraction.
                  x              2
    """

    # TODO!!
