import numpy as np
from typing import Callable, Optional
from dataclasses import dataclass, field
import copy

def ode_quadraticDyn(c: np.ndarray, L: np.ndarray, Q: np.ndarray, t: float, x: np.ndarray) -> np.ndarray:
    """Helper function to compute quadratic dynamics."""
    nx = c.shape[0]

    # Constant part
    dx = c.copy().flatten()
    
    # Linear part
    dx += L @ x
    
    # Quadratic part
    for i in range(nx):
        dx[i] += x.T @ Q[:,:,i] @ x
    
    return dx

@dataclass
class Model_LosslessQuad:
    """Class for Lossless Quadratic dynamical systems."""
    name: str
    nx: int
    
    # Dynamics at the origin: xdot = c + Lx + phi(x)
    c: np.ndarray
    L: np.ndarray
    Ls: np.ndarray = field(init=False)  # Symmetric part of linear dynamics
    La: np.ndarray = field(init=False)  # Asymmetric part of linear dynamics
    Q: np.ndarray
    
    # Dynamics at the coordinate shift ydot = d + Ay + phi(y)
    m: np.ndarray = field(init=False)   # Shifted coordinate, default to be 0
    d: np.ndarray = field(init=False)   # Constant dynamics
    A: np.ndarray = field(init=False)   # Linear dynamics
    As: np.ndarray = field(init=False)  # Symmetric part of linear dynamics
    Aa: np.ndarray = field(init=False)  # Asymmetric part of linear dynamics
    
    # Function handles
    ode: Callable = field(init=False)
    dK0: Callable = field(init=False)
    ode_shifted: Callable = field(init=False)
    dKm: Callable = field(init=False)
    
    def __init__(self, name: str, c: np.ndarray, L: np.ndarray, Q: np.ndarray, 
                 m: Optional[np.ndarray] = None):
        """Initialize the model with basic parameters."""
        self.name = name
        self.nx = c.shape[0]
        self.c = c
        self.L = L
        self.Q = Q
        
        # Compute symmetric and asymmetric parts
        self.Ls = 0.5 * (L + L.T)
        self.La = L - self.Ls
        
        # Setup dynamics
        self.setup_dyn()
        
        # Handle coordinate shift
        if m is None:
            m = np.zeros((self.nx, 1))
        self.func_ShiftSystem(m)
    
    def setup_dyn(self):
        """Setup the dynamic functions."""
        self.ode = lambda t, x: ode_quadraticDyn(self.c, self.L, self.Q, t, x)
        self.dK0 = lambda x: self.c.T @ x + x.T @ self.Ls @ x
    
    def func_ShiftSystem(self, m: np.ndarray):
        """Update the system for a coordinate shift."""
        # Update constant part
        self.m = m
        self.d = np.zeros((self.nx, 1))
        for i in range(self.nx):
            self.d[i] = (self.c[i] + 
                        self.L[i,:] @ m + 
                        m.T @ self.Q[:,:,i] @ m)
        
        # Update linear part
        self.A = np.zeros((self.nx, self.nx))
        for i in range(self.nx):
            self.A[i,:] = self.L[i,:] + 2 * m.T @ self.Q[:,:,i]
        
        # Update symmetric parts
        self.As = self.Ls.copy()
        for i in range(self.nx):
            self.As = self.As - m[i] * self.Q[:,:,i]
        self.Aa = self.A - self.As
        
        # Update ode and power functions
        self.ode_shifted = lambda t, x: ode_quadraticDyn(self.d, self.A, self.Q, t, x)
        self.dKm = lambda x: self.d.T @ x + x.T @ self.As @ x
        
        # Create a copy of self
        return copy.deepcopy(self)
    
    def get_inverseTimeModel(self):
        """Get the inverse-time model."""
        cInv = -self.c
        LInv = -self.L
        QInv = -self.Q
        
        return Model_LosslessQuad(self.name, cInv, LInv, QInv, self.m) 