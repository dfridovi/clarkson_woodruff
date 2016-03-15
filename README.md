# UC Berkeley, Math 221 Final Project
A timing study of low-rank factorization and regression algorithms. The authors are **David Fridovich-Keil** and **Erik Nelson**.

## Summary
We propose to study the computational complexity of the Clarkson-Woodruff algorithm for low-rank factorization and regression on sparse matrices, and compare its performance with other methods. Our project will include the following components:

1. A study of the effects of number of non-zero entries on algorithmic cost
2. A timing comparison against other randomized and non-randomized algorithms for solving the least-squares problem
3. What's the constant in front of the big O?
4. A Monte Carlo analysis of the distribution of error, and a trade-off curve demonstrating the mean error versus number of repetitions of the algorithm
5. Does this error distribution change with matrix size/sparsity/choice of entries in A and b?
