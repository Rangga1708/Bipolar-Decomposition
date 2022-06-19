# Bipolar-Decomposition
It is one of matrix factorization method for nonsingular complex matrix. Based on [(Bhatia, 2013)](https://www.sciencedirect.com/science/article/pii/S0024379513005612), every nonsingular complex matrix $Z$ can be factorized as
$$Z = e^Le^{iT}e^{iK}e^{S}$$
where $S$ and $T$ are real symmetric matrices, and $K$ and $L$ are real skew-symmetric matrices.

For further details and proofs, I suggest you to read the paper. This repository will implement the bipolar decomposition on Python.
