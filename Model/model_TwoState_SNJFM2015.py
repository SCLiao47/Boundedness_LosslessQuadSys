"""
Model of Two State System from Appendix A of Schlegel and Noack, JFM2015.

ODE:
x1dot = x1 - x1*x2
x2dot = -x2 + x1*x1
"""

import numpy as np
from .class_Model_LosslessQuad import Model_LosslessQuad

def model_TwoState_SNJFM2015():
    name = "TwoState_SNJFM2015"
    nx = 2
    c = np.zeros((nx, 1))
    
    L = np.array([
        [1, 0],
        [0, -1]
    ])
    
    Q = np.zeros((nx, nx, nx))
    Q[:,:,0] = np.array([[0, -0.5], [-0.5, 0]])
    Q[:,:,1] = np.array([[1, 0], [0, 0]])

    return Model_LosslessQuad(name, c, L, Q)