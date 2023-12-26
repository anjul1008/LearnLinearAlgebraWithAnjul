
import random

def print_matrix(matrix):
    for row in matrix:
        print([round(val, 2) for val in row])

def generate_random_matrix(rows, cols):
    matrix = [[random.uniform(1, 10) for _ in range(cols)] for _ in range(rows)]
    
    # Ensure linear independence of columns
    for i in range(cols):
        for j in range(i):
            # Subtract a multiple of the j-th column from the i-th column
            scalar = random.uniform(0, 2)
            matrix = [matrix[k][:j] + [(matrix[k][i] - scalar * matrix[k][j])] + matrix[k][j+1:] for k in range(rows)]

    return matrix

def dot_product(v1, v2):
    return sum(x * y for x, y in zip(v1, v2))

def vector_subtraction(v1, v2):
    return [x - y for x, y in zip(v1, v2)]

def vector_scaling(v, scalar):
    return [x * scalar for x in v]

def matrix_multiply(A, B):
    rows_A, cols_A = len(A), len(A[0])
    rows_B, cols_B = len(B), len(B[0])

    result = [[0] * cols_B for _ in range(rows_A)]

    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                result[i][j] += A[i][k] * B[k][j]

    return result








def construct_elementary_matrix(n, row_op_type, i, j, scalar=1):
    """
    Construct an elementary matrix for a given elementary row operation.

    Parameters:
    - n: Size of the square matrix.
    - row_op_type: Type of elementary row operation (e.g., 'swap', 'scale', 'add').
    - i, j: Indices for row operations (e.g., for 'swap' or 'add').
    - scalar: Scalar multiplier for 'scale' operation (default is 1).

    Returns:
    - Elementary matrix.
    """
    E = [[0] * n for _ in range(n)]
    for k in range(n):
        E[k][k] = 1

    if row_op_type == 'swap':
        E[i][i], E[j][j] = 0, 0
        E[i][j], E[j][i] = 1, 1
    elif row_op_type == 'scale':
        E[i][i], E[j][j] = scalar, scalar
    elif row_op_type == 'add':
        E[j][i] = scalar

    return E

def lu_decomposition(A):
    """
    Perform LU decomposition using Crout's algorithm.

    Parameters:
    - A: Symmetric positive definite matrix.

    Returns:
    - L: Lower triangular matrix.
    - U: Upper triangular matrix.
    """
    n = len(A)
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]

    for k in range(n):
        L[k][k] = 1
        U[k][k] = (A[k][k] - sum(L[k][p] * U[p][k] for p in range(k))) / L[k][k]

        for j in range(k + 1, n):
            U[k][j] = (A[k][j] - sum(L[k][p] * U[p][j] for p in range(k))) / L[k][k]

        for i in range(k + 1, n):
            L[i][k] = (A[i][k] - sum(L[i][p] * U[p][k] for p in range(k))) / U[k][k]

    return L, U

# Example usage:
n = 4
A = [
    [4, -2, 4, -2],
    [-2, 5, -2, 5],
    [4, -2, 8, -2],
    [-2, 5, -2, 6]
]

# Initialize the identity matrix as P (permutation matrix)
P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

# Construct elementary matrices for LU decomposition
for i in range(n - 1):
    pivot_index = max(range(i, n), key=lambda k: abs(A[k][i]))
    if pivot_index != i:
        # Swap rows if necessary
        P_elem = construct_elementary_matrix(n, 'swap', i, pivot_index)
        A = matrix_multiply(P_elem, A)
        P = matrix_multiply(P_elem, P)

    for j in range(i + 1, n):
        # Scale factor for the elimination
        scalar = A[j][i] / A[i][i]
        # Construct elementary matrix for the elimination
        E_elem = construct_elementary_matrix(n, 'add', j, i, -scalar)
        A = matrix_multiply(E_elem, A)

L, U = lu_decomposition(A)

print("Elementary matrices:")
print("P =", P)
print("L =", [[round(L[i][j], 2) for j in range(n)] for i in range(n)])
for i in range(n - 1):
    E = construct_elementary_matrix(n, 'add', i + 1, i, -L[i + 1][i])
    print(f"E{i + 1} =", [[round(E[i][j], 2) for j in range(n)] for i in range(n)])

print("\nLU decomposition:")
print("L =", [[round(L[i][j], 2) for j in range(n)] for i in range(n)])




# ==========
def cholesky_decomposition(A):
    n = len(A)
    L = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                sum_val = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = (A[i][j] - sum_val) ** 0.5
            else:
                sum_val = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = (A[i][j] - sum_val) / L[j][j]

    return L


# Example usage:
n = 4
A = [
    [4, -2, 4, -2],
    [-2, 5, -2, 5],
    [4, -2, 8, -2],
    [-2, 5, -2, 6]
]

L = cholesky_decomposition(A)

# Verification: A = LL^T
LT = [[L[j][i] for j in range(n)] for i in range(n)]  # Transpose of L
A_reconstructed = matrix_multiply(L, LT)

print("Original Matrix A:")
for row in A:
    print(row)

print("\nCholesky Decomposition (L):")
for row in L:
    print([round(val, 2) for val in row])

print("\nVerification A = LL^T:")
for row in A_reconstructed:
    print([round(val, 2) for val in row])




# ================

def gram_schmidt(A):
    m, n = len(A), len(A[0])
    Q = [[0] * n for _ in range(m)]
    R = [[0] * n for _ in range(n)]

    for j in range(n):
        v = A[j][:]  # Extract the j-th column as a vector
        for i in range(j):
            R[i][j] = dot_product(Q[i][:], A[j][:])
            v = vector_subtraction(v, vector_scaling(Q[i][:], R[i][j]))

        R[j][j] = sum(x ** 2 for x in v) ** 0.5  # Diagonal element of R
        Q[j][:] = vector_scaling(v, 1 / R[j][j])  # Normalize the vector and assign to Q

    return Q, R



# Example usage:
# import numpy as np

# Generate a matrix A with n linearly independent vectors in m dimensions
m, n = 5, 3
# A = np.random.rand(m, n)
A = generate_random_matrix(m, n)

# Perform QR decomposition
Q, R = gram_schmidt(A)

print("Matrix A:")
print(A)

print("\nMatrix Q:")
print(Q)

print("\nMatrix R:")
print(R)


# ==================


# Example usage:
rows, cols = 5, 4

# Generate a random matrix with linearly independent columns
A = generate_random_matrix(rows, cols)

# Perform QR decomposition
Q, R = gram_schmidt(A)

# Print the matrices
print("Matrix A:")
print_matrix(A)

print("\nMatrix Q:")
print_matrix(Q)

print("\nMatrix R:")
print_matrix(R)

# Observation on the diagonal elements of R
diagonal_elements_R = [R[i][i] for i in range(min(rows, cols))]
print("\nObservation on the diagonal elements of R:")
print(diagonal_elements_R)
