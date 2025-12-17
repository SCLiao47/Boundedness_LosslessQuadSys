import numpy as np
import cvxpy as cp
from typing import Dict, Tuple, Union, Optional
from scipy.linalg import svd

def check_RelTol(a: float, b: float, option: Dict) -> Tuple[bool, float]:
    """Check relative tolerance between two values."""
    RelErr = abs(a-b) / max(abs(a), abs(b))
    ifPass = RelErr < 1e-5
    return ifPass, RelErr

def setup_sdp_problem(As: Union[np.ndarray, cp.Expression], d: Union[np.ndarray, cp.Expression], nx: int) -> Tuple[cp.Problem, cp.Variable, cp.Variable]:
    """
    Set up the SDP optimization problem
    
    Parameters:
    -----------
    As : ndarray or cp.Expression
        System matrix after shifting
    d : ndarray or cp.Expression
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
    Znx1 = np.zeros((nx, 1))
    
    # Setup variables
    gam = cp.Variable(1)
    lam = cp.Variable(1, nonneg=True)

    # Construct the LMI constraint
    # Ensure d is a column vector and d_T is a row vector
    M2 = cp.bmat([[np.array([[0]]), cp.reshape(d, (1, nx), order='C')/2], [cp.reshape(d, (nx, 1), order='C')/2, -As]])
    # Construct the LMI constraint
    M1 = cp.bmat([[cp.reshape(gam, (1,1), order='C'), Znx1.T], [Znx1, -Inx]])
    constraints = [M1 + lam * M2 >> 0]
    
    # Create the problem
    prob = cp.Problem(cp.Minimize(gam), constraints)
    
    return prob, gam, lam

def compute_critical_points(As: np.ndarray, d: np.ndarray, lam: float, tol: float) -> Tuple[np.ndarray, int]:
    """
    Compute the critical points (ystar) of the system
    
    Parameters:
    -----------
    As : ndarray
        System matrix after shifting
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
    nx = As.shape[0]
    Inx = np.eye(nx)
    
    Astar = Inx + lam * As
    
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
        
        print(f"Debug: coef1 = {coef1}")
        print(f"Debug: coef2 = {coef2}")
        print(f"Debug: coef3 = {coef3}")

        s = np.sqrt(coef2**2 - 4*coef1*coef3)
        print(f"Debug: s = {s}")
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
        for i in range(ystar.shape[0]):
            ys = ystar[i]

            # Check norm condition
            ifPass, RelErr = check_RelTol(gam, float(ys.T @ ys), option)
            if not ifPass:
                print(f"Warning: RelTol not satisfied: gam* == ystar'*ystar with relative error {RelErr}")


def func_TRSize_SDP(model: Dict, option: Optional[Dict] = None, m_param: Optional[cp.Parameter] = None) -> Tuple[cp.Problem, cp.Variable, cp.Variable, cp.Expression, cp.Expression]:
    """
    Solving the size of trapping region by QCQP through dual SDP.
    Corresponding to Section 3.2 and 3.3 of Liao et. al, 2024
    
    Args:
        model: Dictionary containing system parameters (nx, c, L, Q, Ls)
        option: Optional dictionary of solver options
        m_param: Optional CVXPY Parameter for the shift vector m. If provided, d and As will be
                 constructed as CVXPY expressions dependent on m_param.
    """
    if option is None:
        option = {
            'verbose': False,
            'tol': 1e-6
        }
    
    # Extract parameters
    nx = model.nx
    c = model.c
    L = model.L
    Ls = model.Ls
    Q = model.Q

    # Determine As and d based on whether m_param is provided
    if m_param is not None:
        # Use func_ShiftSystem to get As and d as CVXPY expressions
        shifted_model = model.func_ShiftSystem(m_param)
        As_to_use = shifted_model.As
        d_to_use = shifted_model.d
    else:
        # Use precomputed As and d from the model (numpy arrays)
        As_to_use = model.As
        d_to_use = model.d
        
        # Check if As is negative definite (only for numpy case)
        try:
            np.linalg.cholesky(-As_to_use)
        except np.linalg.LinAlgError:
            print('As is not negative definite!')
            # This return type needs to be handled carefully, as it's not a problem object
            # For now, return dummy values that will cause an error in bilevel_solver
            return None, None, None, None, None
    
    # Setup SDP
    prob, gam, lam = setup_sdp_problem(As_to_use, d_to_use, nx)
    
    # Return the problem object, gam, and lam for external solving and differentiation
    return prob, gam, lam, As_to_use, d_to_use
    

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