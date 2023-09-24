import numpy as np

def block_lanczos(A, tolerance=1e-10):
    assert A.shape[0] == A.shape[1]
    
    n = A.shape[1]
    V = np.random.randn(n, n)
    T = np.zeros((n,n))

    v = np.random.randn(n)
    v = v.T

    v /= np.linalg.norm(v)
    
    u = np.matmul(A, v)

    V[:,0] = v

    a = np.dot(u.T, v)
    T[0,0] = a

    w = u - a*v

    print(f"v={v}, Av = {u}, V={V}, a={a}, w={w}\n")

    for i in range(1, n):
        b = np.linalg.norm(w)
        
        if b >= tolerance:
            v = w / b
        else:
            v = np.random.randn(n)
            v = v.T
            for j in range(i):
                v -= (np.dot(V[:, j], v) / np.linalg.norm(V[:, j], ord=2)) * V[:, j]  # Gram-Schmidt Process
        
        V[:,i] = v
        u = np.matmul(A, v)
        a = np.dot(u.T, v)
        w = u - a*v - b*V[:,i-1]
        
        T[i,i] = a
        T[i,i-1] = b
        T[i-1,i] = b

        print(f"i={i}: v={v}, u = {u}, V={V}, a={a}\n")

    return V,T

# Example usage:
# Construct a matrix A (you should replace this with your own matrix)
A = np.random.randint(low=1, high=10, size=(6,6))
print(f"A={A}")

# Call the block Lanczos algorithm
V,T = block_lanczos(A)

print(f"V={V}\n T={T}")
