# Understanding the Lanczos Algorithm and Tridiagonalization

The Lanczos algorithm is an iterative method used to approximate the eigenvalues and eigenvectors of a real, symmetric matrix, A. It's particularly useful for large sparse matrices, which can be computationally expensive to diagonalize directly. The algorithm iteratively constructs an orthogonal basis set and a tridiagonal matrix T, which can then be used to find the eigenvalues of A efficiently. This document explains the Lanczos algorithm and highlights how the resulting matrix T is tridiagonal and can be employed in a divide-and-conquer approach to determine the eigenvalues of A.

# Lanczos Algorithm Overview

## Initialization:

* We start with a random initial vector, $\mathbf{v}$.
* Normalize $\mathbf{v}$ to have a unit length.
* Compute $\mathbf{u} = A\mathbf{v}$.
* Set $\mathbf{v}$ as the first column vector of the matrix $V$
* Compute $ a = \mathbf{u}^T\mathbf{v}$
* Set the first diagonal entry of $T$ to $a$
* Compute $w$ as the residual between $u$ and its projection onto $v$ :

```math
\mathbf{w} = \mathbf{u} - a\mathbf{v}
```

## Iteration:

* At each iteration $i$, we perform the following steps:
  
  1. Compute the norm $b$ of the vector $\mathbf{w}_i$.
  2. Confront the norm $b$ with respect to a tolerance (or machine precision) $\epsilon$:
     - If it's bigger, then compute $\mathbf{v}_i = \frac{\mathbf{w}_i}{ b} $.
     - If it's lower, then select a new random vector $\mathbf{v}_i$ orthogonal to all previous $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_i$ using the *Gram-Schmidt* process.
  3. Set $\mathbf{v}_i$ as the $i$-th column of V.
  4. Compute $\mathbf{u}_{i+1}$ from $\mathbf{v}_i$: $\mathbf{u}_{i+1} = A\mathbf{v}_i$
  5. Compute $\mathbf{a}_{i+1} = \mathbf{u}_{i+1}^T\mathbf{v}_i$
  6. Compute $\mathbf{w}_{i+1} = \mathbf{u} - a\mathbf{v_i} - b\mathbf{v}_{i-1}$
  7. Update the tridiagonal matrix $T$ with the computed values $a$ as diagonal entry and the norm $b$ as off-by-1 diagonal entries.

## Repeat:

Repeat the iteration process until convergence, where convergence can be defined based on a tolerance value or a predefined number of iterations.

## Result:

The Lanczos algorithm produces two main results:
An orthogonal basis set, V, whose columns are the vectors generated during the iterations.
A tridiagonal matrix, T, which is symmetric and tridiagonal.

# Tridiagonal Matrix T

The matrix $T$ produced by the Lanczos algorithm is a symmetric tridiagonal matrix with the following properties:

* The diagonal entries ($T[i, i]$) contain approximations of the eigenvalues of the original matrix $A$.
* The off-diagonal entries ($T[i, i+1]$ and $T[i+1, i]$) contain information about the off-diagonal elements of $A$.

The key advantage of having a tridiagonal matrix $T$ is that it simplifies the problem of finding the eigenvalues of A. Once we have T, we can use various numerical techniques to compute the eigenvalues efficiently. One common approach is the divide-and-conquer method.

# Divide-and-Conquer Approach

The divide-and-conquer method is a technique used to find the eigenvalues of a tridiagonal matrix. The code and method presented closely follows https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf

It works as follows:

## Partition:

Given the tridiagonal matrix $T$, we split it into two smaller tridiagonal matrices, $T_1$ and $T_2$, by removing a row and a column at a certain position (usually the middle).
Apply a *rank-1* correction to both $T_1$ and $T_2$: subtract the off-diagonal term $\beta$ that was cut outside the division of $T$ into the two smaller sub-matrices from the last diagonal (bottom-right) element of $T_1$ and the first (top-left) diagonal element of $T_2$.

## Eigenvalue Computation:

Recursively calling the method on the two smaller submatrices $T_1$ and $T_2$ will yield $Q_1, Q_2$ and $D_1, D_2$ respectively the eigenvectors and eigenvalues of submatrices $T_1, T_2$.

## Conquer:

Combine the eigenvalues of $T_1$ and $T_2$ to obtain the eigenvalues of $T$. In particular, let $v^T$ be the real vector composed by concatenating the last row of $Q_1$ with the first row of $Q_2$ and

```math
D = \begin{bmatrix} D_1 &\\ & D_2 \end{bmatrix}
```

then $Q_0\Lambda Q_0^T = D + \beta vv^T$, where $\Lambda$ is the matrix of eigenvalues for $T$ and columns of $Q_0$ are set of orthonormal vectors. Finding the eigenvalues is then a matter of solving the *Secular Equation*, formally:

$$
f(\lambda) := 1 - \beta \sum_{k=1}^{n}{\frac{v_k^2}{\lambda - d_k}} = 0
$$

While the columns of $Q_0$ are easily obtainable as:

$$
q = \frac{(D-\lambda I)^{-1}v}{||(D-\lambda I)^{-1}v||}
$$

Then, since the spectral decomposition of $T$ is:

```math
T = \begin{bmatrix} Q_1 &\\ & Q_2 \end{bmatrix} Q_0\Lambda Q_0^T \begin{bmatrix} Q_1^T &\\ & Q_2^T \end{bmatrix}
```

finding the eigenvectors is a matter of finding

```math
Q = \begin{bmatrix} Q_1 &\\ & Q_2 \end{bmatrix} Q_0
```

## Repeat:

Repeat the partitioning and eigenvalue computation process recursively until you have the eigenvalues of the original matrix $A$ to the desired precision.

## Conclusion

In summary, the Lanczos algorithm is an iterative method for approximating the eigenvalues and eigenvectors of a real, symmetric matrix $A$. It constructs an orthogonal basis set $V$ and a tridiagonal matrix $T$. The tridiagonal matrix $T$ can be used in a divide-and-conquer approach to efficiently compute the eigenvalues of $A$. This technique is particularly useful when dealing with large sparse matrices, where direct diagonalization may be impractical.

