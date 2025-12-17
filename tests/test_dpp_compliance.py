import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
import cvxpy as cp
from Model.class_Model_LosslessQuad import Model_LosslessQuad

def test_dpp_compliance():
    print("\n--- Running DPP Compliance Test ---")

    # 1. Create a dummy Model_LosslessQuad instance
    #    Using arbitrary but valid dimensions and values
    nx = 2
    c = np.array([[1.0], [2.0]])
    L = np.array([[0.1, 0.2], [0.3, 0.4]])
    Q = np.array([[[0.01, 0.02], [0.03, 0.04]], [[0.05, 0.06], [0.07, 0.08]]])
    model = Model_LosslessQuad(name="TestModel", c=c, L=L, Q=Q)

    # 2. Define m_param as a cp.Parameter
    m_param = cp.Parameter(nx, name='m_param')
    m_param.value = np.array([0.5, 0.5]) # Assign an initial numerical value

    # 3. Call model.func_ShiftSystem(m_param) to get As_expr and d_expr
    #    This is where the expressions dependent on m_param are created
    shifted_model = model.func_ShiftSystem(m_param)
    As_expr = shifted_model.As
    d_expr = shifted_model.d

    print(f"Type of As_expr: {type(As_expr)}")
    print(f"Type of d_expr: {type(d_expr)}")

    print(As_expr.is_dpp())
    print(d_expr.is_dpp())

    # 4. Attempt to create a simple CVXPY problem using As_expr and d_expr
    #    We'll try to minimize a simple function of As_expr, for example, its trace
    #    or a sum of its elements, to trigger the DPP check.
    try:
        objective = cp.Minimize(cp.trace(As_expr))
        constraints = [] # No constraints needed for this test
        problem = cp.Problem(objective, constraints)

        print("Attempting to solve a simple problem with requires_grad=True...")
        problem.solve(solver=cp.SCS, requires_grad=True, verbose=True)

        print(f"Problem status: {problem.status}")
        if problem.status == 'optimal':
            print("Simple problem solved successfully with requires_grad=True. DPP compliant.")
        else:
            print("Simple problem did not solve to optimality with requires_grad=True. May indicate DPP issue.")

    except Exception as e:
        print(f"Caught an exception during problem solve: {e}")
        print("This likely confirms the non-DPP compliance of As_expr/d_expr construction.")

    print("\n--- DPP Compliance Test Finished ---")

if __name__ == "__main__":
    test_dpp_compliance()
