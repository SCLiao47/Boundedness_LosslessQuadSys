import numpy as np
import cvxpy as cp
from typing import Dict, Tuple, Union, Optional
from scipy.linalg import svd

def check_RelTol(a: float, b: float, option: Dict) -> Tuple[bool, float]:
    """Check relative tolerance between two values."""
    RelErr = abs(a-b) / max(abs(a), abs(b))
    ifPass = RelErr < option['tol']
    return ifPass, RelErr

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
    As = model.As
    d = model.d
    nx = model.nx
    
    # Check if As is negative definite
    try:
        np.linalg.cholesky(-As)
    except np.linalg.LinAlgError:
        print('As is not negative definite!')
        return float('inf'), {}
    
    # Utility matrices
    Inx = np.eye(nx)
    Znx = np.zeros((nx, nx))
    Znx1 = np.zeros((nx, 1))
    
    # Lagrangian dual SDP formulation
    flag_Lag = False
    
    # Setup and solve the SDP
    gam = cp.Variable(1)
    lam = cp.Variable(1, nonneg=True)
    
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
    
    # Solve the problem
    prob = cp.Problem(cp.Minimize(gam), constraints)
    try:
        prob.solve(solver=cp.MOSEK, verbose=False)
    except:
        prob.solve(verbose=False)
    
    cvx_status_Lag = prob.status
    
    if cvx_status_Lag == 'optimal':
        flag_Lag = True
        
        # Check Astar <= 0
        Astar = Inx + lam.value * As
        # Check d in range(Astar)
        assert np.linalg.norm(Astar @ d) >= np.finfo(float).eps, "d is not in range of Astar!"
        
        # Compute ystar using KKT condition
        y0 = -lam.value/2 * np.linalg.solve(Astar, d)
        
        # SVD decomposition
        _, S, V = svd(Astar)
        r = np.sum(np.abs(S) > option['tol'])
        v = V[:, r:].T
        
        if r == nx:
            ystar = y0
            ystar = ystar.reshape(1, -1, 1)

            # Verify ystar has shape (1, nx, 1)
            assert ystar.shape == (1, nx, 1), f"ystar shape {ystar.shape} does not match expected shape (1, {nx}, 1)"
        elif r == nx-1:
            # Solve CS, which is a quadratic equation in c
            coef1 = v @ As @ v.T
            coef2 = 2 * v @ As @ y0 + d.T @ v.T
            coef3 = y0.T @ As @ y0 + d.T @ y0
            
            s = np.sqrt(coef2**2 - 4*coef1*coef3)
            c = (-coef2 + np.array([s, -s])) / (2*coef1)
            
            ystar = y0 + v.T @ c
            
            # Verify ystar has shape (2, nx, 1)
            assert ystar.shape == (2, nx, 1), f"ystar shape {ystar.shape} does not match expected shape (2, {nx}, 1)"
        else:
            # rank(Astar) <= nx - 2
            ystar = f'A {nx-r}-dimensional sphere.'

        # Debug
        # print(f"ystar: {ystar}")
        # print(f"r: {r}")
        # print(f"Astar: {Astar}")
        # print(f"d: {d}")
        # print(f"gam: {gam.value}")
        # print(f"lam: {lam.value}")

        # ystar is a 3D array with shape (nsolutions, nx, 1)
        # nsolutions is the number of solutions
        # ystar[i] is the i-th solution with (nx,1) as the initial condition
        
        # Check ystar solutions
        if r >= nx-1:
            def Lag(y, lam):
                # Ensure y is a column vector
                y = y.reshape(-1, 1) if len(y.shape) == 1 else y
                return float(y.T @ y + lam * (y.T @ As @ y + d.T @ y))
            
            if not isinstance(ystar, str):
                for i in range(ystar.shape[0]):
                    ys = ystar[i]
                    
                    # Check gam* = L(ystar, lam*)
                    assert check_RelTol(gam.value, float(Lag(y0, lam.value)), option)[0], \
                        "Error: gam* == L(ystar, lam*)"
                    assert check_RelTol(gam.value, float(Lag(ys, lam.value)), option)[0], \
                        "Error: gam* == L(ystar, lam*)"
                    
                    # Check ystar'*ystar == L(ystar, lam*)
                    ifPass, RelErr = check_RelTol(gam.value, float(ys.T @ ys), option)
                    if not ifPass:
                        print(f"Warning: RelTol not satisfied: gam* == ystar'*ystar with relative error {RelErr}")
        
        # Set output
        rTrap = float(np.sqrt(gam.value))
        info = {
            'feasibility': True,
            'cvx_Lag': cvx_status_Lag,
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