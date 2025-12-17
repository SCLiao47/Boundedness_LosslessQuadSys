import numpy as np
from .class_Model_LosslessQuad import Model_LosslessQuad

def model_TwoState_Alignment(ifAlign: bool = False):
    """Model of a two-state system to illustrate the conservatism of the Schlegel and Noack method.

    Args:
        ifAlign: If True, the constant input c is aligned with the eigenvector corresponding to the larger eigenvalue of L.
    """
    nx = 2
    if ifAlign:
        name = "TwoStateBD_Aligned"
        c = np.zeros((nx, 1))
    else:
        name = "TwoStateBD_NotAligned"
        c = np.array([[0], [1]])

    L = np.diag([-1, -4])

    Q1 = np.array([[0, -0.5], [-0.5, 0]])
    Q2 = np.array([[1, 0], [0, 0]])
    Q = np.stack([Q1, Q2], axis=-1)

    return Model_LosslessQuad(name, c, L, Q)
