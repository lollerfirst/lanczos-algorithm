import numpy as np
import random
from math import abs, sqrt


# machine precision
epsilon = 1e-10

def solve_base(A):
    # det(A-lI) = 0 => (A[0,0]-l)*(A[1,1]-l)-A[1,0]*A[0,1] = 0
    a = A[0,0], b = A[0,1], c = A[1,0], d = A[1,1]

    eig = np.array([((a+d) + sqrt((a+d)**2 - 4 * (a*d-b*c))) / 2,
                    ((a+d) - sqrt((a+d)**2 - 4 * (a*d-b*c))) / 2])
    Q = np.empty((2,2))

    for i in range(2):
        F = np.copy(A)
        F[0,0] -= eig[i]
        F[1,1] -= eig[i]

        # Mini Gaussian Elimination
        k = F[1,0] / F[0,0]
        F[1,:] = -k * F[0,:] + F[1,:]
        
        if abs(F[1,1]) > epsilon:
        
            Q[1,i] = 0
            Q[0,i] = 0
        
        else:
            Q[1,i] = 1
            Q[0,i] = -F[0,(i+1)&1] / F[0, i]

    return Q, eig
    

def psi(v, d, l, start, stop):
    x = 0
    for k in range(start, stop+1):
        x += (v[k]**2) / (d[start] - l)
    return x

# derivative of psi wrt l (the eigenvalue)
def dpsi_dl(v, d, l, start, stop):
    x = 0
    for k in range(start, stop+1):
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

    # (D + bvv^T)l = 0 ==> 1 + bv^T(D - lI)^-1 = 0 ==> 1 - b SUM_k=1^n (v_k^2 / (l - d_k)) = 0
    eig = np.empty(n)
    X = np.empty((n,n))

    # TODO introduce deflation:
    #    - Checking of zeros in v^T
    #    - Apply Givens' rotations to v^T and D to form zeros

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
            
            tmp = b*dpsi_dl(v, d, l, i+1, n-1)
            c_2 = tmp*((d[i+1] - l)**2)
            c_2_hat = b*psi(v, d, l, i+1, n-1) - tmp*(d[i+1] - l)     # (5.29)

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