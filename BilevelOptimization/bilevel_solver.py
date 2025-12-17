import numpy as np
import cvxpy as cp
from typing import Dict, Tuple, Optional
from BoundednessAnalysis.func_TRSize import func_TRSize_SDP
from Model.class_Model_LosslessQuad import Model_LosslessQuad

def solve_bilevel_optimization(model: Model_LosslessQuad, initial_m: np.ndarray, 
                                 learning_rate: float = 0.01, num_iterations: int = 100,
                                 verbose: bool = True) -> Tuple[np.ndarray, float]:
    """
    Solves the bilevel optimization problem to minimize rTrap by optimizing m.
    
    The upper-level problem minimizes rTrap (output of func_TRSize_SDP).
    The lower-level problem is the func_TRSize_SDP itself.
    
    Args:
        model: The Model_LosslessQuad object.
        initial_m: Initial guess for the shift vector m.
        learning_rate: Step size for gradient descent.
        num_iterations: Maximum number of iterations.
        verbose: Whether to print iteration details.
        
    Returns:
        Tuple of (optimal_m, min_rTrap).
    """
    
    current_m = initial_m.copy()
    
    # Define m as a CVXPY Parameter for differentiation
    m_param = cp.Parameter(model.nx, name='m_param')
    
    for i in range(num_iterations):
        m_param.value = current_m.flatten()
        
        # Get the CVXPY problem object and parameters from func_TRSize_SDP
        prob, gam, lam, As_expr, d_expr = func_TRSize_SDP(model, m_param=m_param, option={'verbose': False})
        
        # If func_TRSize_SDP returned None (due to infeasibility), break
        if prob is None:
            print(f"Iteration {i}: Lower-level problem setup failed. Stopping.")
            break

        # Solve the lower-level problem
        try:
            prob.solve(solver=cp.SCS, verbose=False, requires_grad=True)
        except Exception as e:
            print(f"SCS solve failed: {e}. Stopping.")
            break
            
        print(f"Iteration {i}: Problem status: {prob.status}")
        if prob.status != 'optimal':
            print(f"Iteration {i}: Lower-level problem not optimal. Stopping.")
            break
            
        r_trap = np.sqrt(gam.value.item())

        # Compute the gradient using CVXPY's derivative
        grad_dict = prob.derivative()
        
        # The gradient of gam with respect to m_param is obtained directly
        grad_r_trap = grad_dict[m_param].reshape(-1, 1) # Reshape to column vector

        # Upper-level problem: Update m using gradient descent
        new_m = current_m - learning_rate * grad_r_trap
        
        # Check for convergence (simple L2 norm of change in m)
        if np.linalg.norm(new_m - current_m) < 1e-6:
            if verbose:
                print(f"Iteration {i}: Converged. rTrap = {r_trap:.6f}, m = {current_m.flatten()}")
            current_m = new_m
            break
            
        current_m = new_m
        
        if verbose:
            print(f"Iteration {i}: rTrap = {r_trap:.6f}, m = {current_m.flatten()}")
            
    return current_m, r_trap
