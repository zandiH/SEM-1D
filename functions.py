import numpy as np

def connect_mat(ne, N):
    '''
    Connectivity Matrix, in which rows are based on
    (ne - 1)*N + Var, where Var is changing from 1:N+1
    '''
    C = np.zeros((N+1, ne))
    for i in range(ne):
        for j in range(N+1):
            C[j, i] = i * N + (j+1)
    return C
