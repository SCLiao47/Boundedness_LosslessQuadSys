# Development Plan

## Current Goals:
1.  [DONE] Enable differentiable `func_TRSize_SDP` (initial setup).
2.  [DONE] Develop the initial bilevel optimization framework using finite differences.
3.  [DONE] Initialize `m` for bilevel optimization using `func_findNDShifting.py`.
4.  [IN PROGRESS] Improve gradient calculation in `bilevel_solver.py` by replacing finite differences with `cvxpy`'s native differentiable programming capabilities.
    *   **Current Challenge:** Encountering "Problem is not DPP" errors when using `requires_grad=True` with CVXPY. The issue is likely that the operations constructing `d` and `As` within `Model.func_ShiftSystem` (when `m_param` is a CVXPY object) are not DPP-compliant. Throught the test script `tests/test_dpp_compliance.`, `d` is indeed incompatible. There might be a way to get around with the quadratic form. See https://www.cvxpy.org/tutorial/dpp/index.html.
    *   **Next Steps:** Create a dedicated testing script to isolate and confirm if the construction of `d` and `As` is indeed the source of the non-DPP error. Based on the findings, refactor `Model.func_ShiftSystem` and `func_TRSize_SDP` to ensure DPP compliance.

## Future Tasks:
*   Expand testing to include gradient verification and bilevel convergence for more complex scenarios.
*   Investigate and address persistent MOSEK solver warnings if they impact functionality or accuracy.
*   Ensure the nonlinearity in the model adheres to the lossless properties as described in the TRCVX.pdf paper, especially before utilizing `cvxpy`'s gradient calculation for differentiable optimization.
