import numpy as np
import cvxpy as cp
from typing import Dict, Tuple, Union, Optional
from scipy.linalg import svd

def check_RelTol(a: float, b: float, option: Dict) -> Tuple[bool, float]:
    """Check relative tolerance between two values."""
    RelErr = abs(a-b) / max(abs(a), abs(b))
    ifPass = RelErr < option['tol']
    return ifPass, RelErr

def setup_sdp_problem(As: np.ndarray, d: np.ndarray, nx: int) -> Tuple[cp.Problem, cp.Variable, cp.Variable]:
    """
    Set up the SDP optimization problem
    
    Parameters:
    -----------
    As : ndarray
        System matrix after shifting
    d : ndarray
        Linear term vector
    nx : int
        System dimension
    
    Returns:
    --------
    prob : cp.Problem
        CVXPY optimization problem
    gam : cp.Variable
        Variable to be minimized
    lam : cp.Variable
        Lagrange multiplier
    """
    # Utility matrices
    Inx = np.eye(nx)
    Znx = np.zeros((nx, nx))
    Znx1 = np.zeros((nx, 1))
    
    # Setup variables
    gam = cp.Variable(1)
    lam = cp.Variable(1, nonneg=True)

    """
    TODO
    - [ ] Set up m as the parameter of CVXpy
    - [ ] compute d and As accordingly
    - [ ] extract the gradient of the SDP solution w.r.t. m
    """
    
    # Construct the LMI constraint
    lmi_left = cp.bmat([
        [cp.reshape(gam, (1,1)), Znx1.T],
        [Znx1, -Inx]
    ])
    
    lmi_right = cp.bmat([
        [cp.Constant([[0]]), d.T/2],
        [d/2, -As]
    ])
    
    constraints = [lmi_left + lam * lmi_right >> 0]
    
    # Create the problem
    prob = cp.Problem(cp.Minimize(gam), constraints)
    
    return prob, gam, lam

def compute_critical_points(Astar: np.ndarray, d: np.ndarray, lam: float, tol: float) -> Tuple[np.ndarray, int]:
    """
    Compute the critical points (ystar) of the system
    
    Parameters:
    -----------
    Astar : ndarray
        Modified system matrix
    d : ndarray
        Linear term vector
    lam : float
        Optimal Lagrange multiplier
    tol : float
        Tolerance for rank computation
    
    Returns:
    --------
    ystar : ndarray or str
        Critical points or description if points form a sphere
    r : int
        Rank of Astar
    """
    nx = Astar.shape[0]
    
    # Compute base solution
    y0 = -lam/2 * np.linalg.solve(Astar, d)
    
    # SVD decomposition for rank
    _, S, V = svd(Astar)
    r = np.sum(np.abs(S) > tol)
    v = V[:, r:].T
    
    if r == nx:
        ystar = y0.reshape(1, -1, 1)
    elif r == nx-1:
        # Solve quadratic equation for additional solutions
        coef1 = v @ Astar @ v.T
        coef2 = 2 * v @ Astar @ y0 + d.T @ v.T
        coef3 = y0.T @ Astar @ y0 + d.T @ y0
        
        s = np.sqrt(coef2**2 - 4*coef1*coef3)
        c = (-coef2 + np.array([s, -s])) / (2*coef1)
        
        ystar = y0 + v.T @ c
        ystar = ystar.reshape(2, -1, 1)
    else:
        ystar = f'A {nx-r}-dimensional sphere.'
    
    return ystar, r

def verify_solutions(ystar: np.ndarray, gam: float, lam: float, As: np.ndarray, d: np.ndarray, option: Dict):
    """
    Verify the computed solutions satisfy optimality conditions
    
    Parameters:
    -----------
    ystar : ndarray
        Computed critical points
    gam : float
        Optimal objective value
    lam : float
        Optimal Lagrange multiplier
    As : ndarray
        System matrix
    d : ndarray
        Linear term vector
    option : dict
        Options dictionary containing tolerance settings
    """
    def Lag(y, lam):
        y = y.reshape(-1, 1) if len(y.shape) == 1 else y
        return float(y.T @ y + lam * (y.T @ As @ y + d.T @ y))
    
    if not isinstance(ystar, str):
        y0 = -lam/2 * np.linalg.solve(np.eye(As.shape[0]) + lam * As, d)
        for i in range(ystar.shape[0]):
            ys = ystar[i]
            
            # Verify optimality conditions
            assert check_RelTol(gam, float(Lag(y0, lam)), option)[0], \
                "Error: gam* == L(ystar, lam*)"
            assert check_RelTol(gam, float(Lag(ys, lam)), option)[0], \
                "Error: gam* == L(ystar, lam*)"
            
            # Check norm condition
            ifPass, RelErr = check_RelTol(gam, float(ys.T @ ys), option)
            if not ifPass:
                print(f"Warning: RelTol not satisfied: gam* == ystar'*ystar with relative error {RelErr}")

def func_TRSize_SDP(model: Dict, option: Optional[Dict] = None) -> Tuple[float, Dict]:
    """
    Solving the size of trapping region by QCQP through dual SDP.
    Corresponding to Section 3.2 and 3.3 of Liao et. al, 2024
    """
    if option is None:
        option = {
            'verbose': False,
            'tol': 1e-6
        }
    
    # Extract parameters
    nx = model.nx
    #  original coordinate
    c = model.c
    Ls = model.Ls
    #  shifted coordinate
    m = model.m
    d = model.d
    As = model.As
    
    # Check if As is negative definite
    try:
        np.linalg.cholesky(-As)
    except np.linalg.LinAlgError:
        print('As is not negative definite!')
        return float('inf'), {}
    
    # Setup and solve SDP
    prob, gam, lam = setup_sdp_problem(As, d, nx)
    
    try:
        prob.solve(solver=cp.MOSEK, verbose=False)
    except:
        prob.solve(verbose=False)
    
    if prob.status == 'optimal':
        # Compute critical points
        Astar = np.eye(nx) + lam.value * As
        ystar, r = compute_critical_points(Astar, d, lam.value, option['tol'])
        
        # Verify solutions
        if r >= nx-1:
            verify_solutions(ystar, gam.value, lam.value, As, d, option)
        
        # Set output
        rTrap = float(np.sqrt(gam.value))
        info = {
            'feasibility': True,
            'cvx_Lag': prob.status,
            'gam': gam.value,
            'lam': lam.value,
            'ystar': ystar,
            'rank': r
        }
        
        if option['verbose']:
            print('Trapping region found!')
            print(f'TR size = {rTrap:.3f}')
            print('y* =')
            print(ystar if isinstance(ystar, str) else ystar.T)
        
        return rTrap, info
    
    else:
        # SDP failed to find a solution, return NaN and empty info
        return float('nan'), {'feasibility': False, 'cvx': None}
    

def func_TRSize_SN(model: Dict) -> Tuple[float, Dict]:
    """
    Estimate the size of trapping region using Schlegel and Noack's method (worst-case spectral analysis).
    """
    # Get largest real eigenvalue
    lam1 = np.max(np.real(np.linalg.eigvals(model.As)))
    
    if lam1 >= 0:
        # TR doesn't exist
        rTrap = float('inf')
    else:
        # Trapping region exists
        rTrap = np.linalg.norm(model.d) / abs(lam1)
    
    return rTrap

def func_TRSize(model: Dict, method: str = 'SDP', option: Optional[Dict] = None) -> Tuple[float, Dict]:
    """
    Compute the size of trapping region using either SDP or SN method.
    
    Args:
        model: Dictionary containing system parameters (As, d, nx)
        method: 'SDP' for QCQP through dual SDP (default) or 'SN' for Schlegel-Noack method
        option: Optional dictionary of solver options
    
    Returns:
        rTrap: Radius of trapping region
        info: Dictionary containing solver information
    """
    if method.upper() == 'SDP':
        return func_TRSize_SDP(model, option)
    elif method.upper() == 'SN':
        return func_TRSize_SN(model)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'SDP' or 'SN'.")