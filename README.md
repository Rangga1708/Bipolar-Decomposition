# Bipolar-Decomposition
It is one of matrix factorization method for nonsingular complex matrix. Based on [(Bhatia, 2013)](https://www.sciencedirect.com/science/article/pii/S0024379513005612), every nonsingular complex matrix $Z$ can be factorized as
$$Z = e^Le^{iT}e^{iK}e^{S}$$
where $S$ and $T$ are real symmetric matrices, and $K$ and $L$ are real skew-symmetric matrices.

For further details and proofs, I suggest you to read the paper. Although it's rare to find complex matrix in our daily life, maybe it can be used on quantum mechanics field since it uses complex matrix (as far as I know). 

This repository will implement the bipolar decomposition on Python.

## How to Use?
Since it's a local library, you don't need to do `pip install` or anything. Just make sure the `bipolar_decomposition.py` file is in the same folder as the code who's calling it.

Let's import some libraries we need first.
```
import numpy as np
import scipy.linalg as la
import bipolar_decomposition as bd
```

Define the matrix we want to decompose.
```
Z = np.matrix([[1j,0,-1j],
               [0,1,-1-4j],
               [2-1j,1j,3]])
```

Then, let's decompose it.
```
L,T,K,S = bd.bipolar(Z)
```

If you want to check the results, just multiply them all (don't forget to exponentiate them first).
```
la.expm(L) @ la.expm(1j * T) @ la.expm(1j * K) @ la.expm(S)
```

Further details on how to use it can be checked on `Documentation.ipynb` notebook.

Last but not least, hope it will be useful. :)
