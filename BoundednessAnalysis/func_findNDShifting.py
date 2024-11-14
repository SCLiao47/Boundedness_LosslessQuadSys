import numpy as np
import cvxpy as cp
from typing import Dict, Tuple, Optional

def func_findNDShifting(model: Dict, option: Optional[Dict] = None) -> Tuple[np.ndarray, Dict]:
    """
    Solve trapping region condition using SDP.
    Corresponding to Section 3.1 of Liao et. al, 2024
    """
    # Default options if not provided
    if option is None:
        option = {
            'round_Ndigit': 3,
            'verbose': False,
            'tol': 1e-6
        }
    
    # Initialize output
    m = np.nan
    info = {'feasibility': False}
    
    # System parameters
    nx = model.nx
    Ls = model.Ls
    Q = model.Q
    
    Inx = np.eye(nx)
    
    # Filtering the data matrix
    if not np.isnan(option['round_Ndigit']):
        Ls = np.round(Ls, option['round_Ndigit'])
        Q = np.round(Q, option['round_Ndigit'])
    
    # Solve by SDP
    # This SDP pushes As most further to negative definite
    a = cp.Variable(1)
    m_var = cp.Variable((nx, 1))
    
    As = Ls
    for i in range(nx):
        As = As - m_var[i] * Q[:,:,i]
    
    constraints = [As <= a * Inx]
    prob = cp.Problem(cp.Minimize(a), constraints)
    
    # Try solving with MOSEK
    try:
        prob.solve(solver=cp.MOSEK, verbose=False)
    except:
        # Fallback to default solver if MOSEK not available
        prob.solve(verbose=False)
    
    # Handle unbounded case
    if prob.status == 'unbounded':
        print('Warning: Unregularized TRSDP is unbounded below(a*=-inf). Rerun TRSDP with regularization')
        
        m_bounded = cp.Variable((nx, 1))
        As = Ls
        for i in range(nx):
            As = As - m_bounded[i] * Q[:,:,i]
        
        constraints = [As <= -1e-3 * Inx]
        prob = cp.Problem(cp.Minimize(cp.norm(m_bounded)), constraints)
        
        try:
            prob.solve(solver=cp.MOSEK, verbose=False)
        except:
            prob.solve(verbose=False)
            
        m = m_bounded.value
        a_val = np.max(np.linalg.eigvals(As.value))
    else:
        m = m_var.value
        a_val = a.value[0] if hasattr(a.value, '__len__') else a.value
    
    # Setting output
    if prob.status == 'optimal':
        info['feasibility'] = True
        info['a'] = a_val
        info['As'] = As.value
        info['W'] = constraints[0].dual_value
        
        if a_val < 0:
            info['existTR'] = True
            if option['verbose']:
                print('Coordinated shift s.t. As<0 is found:')
                print(m.T)
                print(f'Most positive eigenvalue a = {a_val:.4f}')
        else:
            info['existTR'] = False
    else:
        if option['verbose']:
            print(f'cvx_status: {prob.status}')
            print('could not find m such that As<0')
    
    info['m'] = m
    return m, info 