import numpy as np
import random
from math import abs

epsilon = 1e-10

def psi(v, d, l, start, stop):
    x = 0
    for k in range(start, stop):
        x += (v[k]**2) / (d[start] - l)
    return x

# derivative of psi wrt l (the eigenvalue)
def dpsi_dl(v, d, l, start, stop):
    x = 0
    for k in range(start, stop):
        x += (v[k]**2) / ((d[start] - l)**2)
    return x

# Reference: https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf

def dnc_eigenvalues(T: np.matrix):
    [n, m] = T.shape

    assert n == m, "Not a square matrix"
    
    # Base case
    if n <= 2:
        return solve_base(T)

    r = n >> 1
    
    assert T[r, r-1] == T[r-1, r], "Not a tridiagonal Matrix!"
    b = T[r, r-1]

    # Divide
    T1 = T[:r, :r]
    T2 = T[r:, r:]

    # rank-1 correction:
    T1[r-1, r-1] -= b
    T2[0, 0] -= b

    # Recursion
    [Q1, d1] = dnc_eigenvalues(T1)
    [Q2, d2] = dnc_eigenvalues(T2)

    # Conquer
    d = np.hstack((d1, d2))
    v = np.hstack((Q1[-1,:], Q2[0,:]))

    [d, v] = sorted(zip(d, v), lambda x: x[0]) # Sort the diagonal values

    # (D + bzz^T)v = 0 ==> 1 + bz^T(D - vI)^-1 = 0 ==> 1 - b SUM_k=1^n (z_k^2 / (v - d_k)) = 0
    eig = np.empty(n)
    X = np.empty((n,n))

    for i in range(n):

        # We know the eigenvalues reside in between intervals d[0] < l_0 < d[1] < l_1 < ... < l_(n-1)
        # if b is positive, otherwise: l_0 < d[0] < l_1 < d[1] < ... < d[n-1]
        if b > 0:
            if i < n-1:
                l = random.uniform(d[i], d[i+1])
            else:
                l = random.uniform(d[n-1], d[n-1]+1)
        else:
            if i > 0:
                l = random.uniform(d[i-1], d[i])
            else:
                l = random.uniform(d[0]-1, d[0])

        
        err = float("inf")          # Initialize error to infinity

        while abs(err) >= epsilon:

            tmp = b*dpsi_dl(v, d, l, 0, i)
            c_1 = tmp*((d[i] - l)**2)
            c_1_hat = b*psi(v, d, l, 0, i) - tmp*(d[i] - l)         # (5.28)
            
            tmp = b*dpsi_dl(v, d, l, i+1, n)
            c_2 = tmp*((d[i+1] - l)**2)
            c_2_hat = b*psi(v, d, l, i+1, n) - tmp*(d[i+1] - l)     # (5.29)

            c_3 = 1 + c_1_hat + c_2_hat
            err = c_3 + c_1 / (d[i] - l) + c_2 / (d[i+1] - l)       # (5.30) see also (5.25) (5.26)

            # Newton iteration l_(n+1) = l_n - h(l_n) / h'(l_n)
            l = l - err / ( c_1 / ((d[i] - l)**2) + c_2 / ((d[i+1] - l)**2) )
        
        # We now have eigenvalue l
        eig[i] = l

        # Calculate eigenvector for l
        x = np.array([(1/(l - d[k]))*v[k] for k in range(n)])       # (5.15)
        x /= np.linalg.norm(x)

        # Record eigenvector into Q
        X[:,i] = x.T

    Q = np.zeros(n,n)
    Q[:r, :r] = Q1
    Q[r:, r:] = Q2
    Q = np.matmul(Q, X)

    return Q, eig