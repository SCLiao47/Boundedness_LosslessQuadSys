import numpy as np
from .class_Model_LosslessQuad import Model_LosslessQuad

def model_Lorenz():
    """Model of Lorenz Attractor.
    dx/dt = c + L*x + [x'*Q1*x, ..., x'*Qn*x]'
    """
    name = "Lorenz_Chaotic"
    nx = 3
    c = np.zeros((nx, 1))
    
    sig = 10
    rho = 28
    bet = 8/3
    L = np.array([
        [-sig, sig, 0],
        [rho, -1, 0],
        [0, 0, -bet]
    ])
    
    Q = np.zeros((nx, nx, nx))
    Q[:,:,0] = 0
    Q[:,:,1] = np.array([
        [0, 0, -1/2],
        [0, 0, 0],
        [-1/2, 0, 0]
    ])
    Q[:,:,2] = np.array([
        [0, 1/2, 0],
        [1/2, 0, 0],
        [0, 0, 0]
    ])
    
    return Model_LosslessQuad(name, c, L, Q) 