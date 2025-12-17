import unittest
import numpy as np
from BoundednessAnalysis.func_findNDShifting import func_findNDShifting
from BoundednessAnalysis.func_TRSize import func_TRSize_SDP, func_TRSize_SN
from Model.class_Model_LosslessQuad import Model_LosslessQuad

class TestBoundednessAnalysis(unittest.TestCase):

    def test_func_findNDShifting(self):
        # Create a simple 2D system for which we know the answer
        name = "TestSystem"
        nx = 2
        c = np.array([[1], [1]])
        L = np.array([[-1, 0], [0, -2]])
        Q = np.zeros((nx, nx, nx))
        model = Model_LosslessQuad(name, c, L, Q)

        # Set options
        option = {
            'round_Ndigit': 3,
            'verbose': False,
            'tol': 1e-6
        }

        # Call the function
        m, info = func_findNDShifting(model, option)

        # Assertions
        self.assertTrue(info['feasibility'])
        self.assertTrue(info['existTR'])
        self.assertLess(info['a'], 0)
        self.assertIsNotNone(m)

    def test_func_TRSize_SDP(self):
        # Create a simple 2D system (similar to the one in the paper's example)
        name = "TestTRSystem"
        nx = 2
        c = np.array([[0], [1]])
        L = np.array([[-1, 0], [0, -4]])
        Q = np.zeros((nx, nx, nx))
        model = Model_LosslessQuad(name, c, L, Q)

        # Shift the model by 0 vector (as done in the script)
        m = np.array([[0], [0]])
        model_shifted = model.func_ShiftSystem(m)

        # Set options
        option = {
            'verbose': False,
            'tol': 1e-6
        }

        # Call the function
        r_sdp, info_sdp = func_TRSize_SDP(model_shifted, option)

        # Assertions (based on the paper's example for the two-state system)
        self.assertTrue(info_sdp['feasibility'])
        self.assertAlmostEqual(r_sdp, 0.288675, places=5) # Value from paper
        self.assertIsNotNone(info_sdp['ystar'])

    def test_func_TRSize_SN(self):
        # Create a simple 2D system (similar to the one in the paper's example)
        name = "TestTRSystem"
        nx = 2
        c = np.array([[0], [1]])
        L = np.array([[-1, 0], [0, -4]])
        Q = np.zeros((nx, nx, nx))
        model = Model_LosslessQuad(name, c, L, Q)

        # Shift the model by 0 vector (as done in the script)
        m = np.array([[0], [0]])
        model_shifted = model.func_ShiftSystem(m)

        # Call the function
        r_sn = func_TRSize_SN(model_shifted)

        # Assertions (based on the paper's example for the two-state system)
        self.assertAlmostEqual(r_sn, 1.0, places=5) # Value from paper

    def test_func_TRSize_SDP_with_m_param(self):
        import cvxpy as cp
        # Create a simple 2D system
        name = "TestTRSystemWithParam"
        nx = 2
        c = np.array([[0], [1]])
        L = np.array([[-1, 0], [0, -4]])
        Q = np.zeros((nx, nx, nx))
        model = Model_LosslessQuad(name, c, L, Q)

        # Define m as a CVXPY Parameter
        m_param = cp.Parameter((nx, 1), name='m_param')
        m_param.value = np.array([[0.1], [0.2]]) # Assign an initial value

        # Set options
        option = {
            'verbose': False,
            'tol': 1e-6
        }

        # Call the function with m_param
        r_sdp, info_sdp = func_TRSize_SDP(model, option, m_param=m_param)

        # Assertions
        self.assertTrue(info_sdp['feasibility'])
        self.assertIsInstance(r_sdp, float)
        self.assertFalse(np.isnan(r_sdp))
        self.assertIsNotNone(info_sdp['gam'])
        self.assertIsNotNone(info_sdp['lam'])
        self.assertIsNone(info_sdp['ystar']) # Should be None when m_param is used
        self.assertIsNone(info_sdp['rank'])  # Should be None when m_param is used

if __name__ == '__main__':
    unittest.main()