import numpy as np
import random

epsilon = 1e-10

def psi(v, d, l, start, stop):
    x = 0
    for k in range(start, stop):
        x += (v[k]**2) / (d[start] - l)
    return x

def dpsi_dl(v, d, l, start, stop):
    x = 0
    for k in range(start, stop):
        x += (v[k]**2) / ((d[start] - l)**2)
    return x

# https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf

def dnc_eigenvalues(T: np.matrix):
    [n, m] = T.shape

    assert n == m, "Not a square matrix"
    
    # Base case
    if n <= 2:
        [Q1, T1] = solve_base(T)


    r1 = n >> 1
    r2 = m >> 1

    assert T[r1, r2-1] == T[r1-1, r2], "Not a tridiagonal Matrix!"
    b = T[r1, r2-1]

    # Divide
    T1 = T[:r1, :r2]
    T2 = T[r1:, r2:]

    # rank-1 correction:
    T1[r1-1, r2-1] -= b
    T2[0, 0] -= b

    # Recursion
    [Q1, d1] = dnc_eigenvalues(T1)
    [Q2, d2] = dnc_eigenvalues(T2)

    # Conquer
    d = np.hstack((d1, d2))
    v = np.hstack((Q1[-1,:], Q2[0,:]))

    [d, v] = sorted(zip(d, v), lambda x: x[0]) # Sort the diagonal values

    # (D + bvv^T)x = 0 ==> 1 + bv^T(D - lI)^-1 = 0 ==> 1 - b SUM_k=1^n (v_k^2 / (l - d_k)) = 0
    # Search in between the intervals of values of d
    eigenvalues = np.empty(d.size)

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

        
        err = float(0x7FFFFFFF)

        tmp = b*dpsi_dl(v, d, l, 0, i)
        c_1 = tmp*((d[i] - l)**2)
        c_1_hat = b*psi(v, d, l, 0, i) - tmp*(d[i] - l)
        
        tmp = b*dpsi_dl(v, d, l, i+1, n)
        c_2 = tmp*((d[i+1] - l)**2)
        c_2_hat = b*psi(v, d, l, i+1, n) - tmp*(d[i+1] - l)

        c_3 = 1 + c_1_hat + c_2_hat
        
        while err > epsilon:
            # ... TBC        
    

    
    return Q, eigenvalues
