import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from Model.class_Model_LosslessQuad import Model_LosslessQuad
from Model.model_TwoState_SNJFM2015 import model_TwoState_SNJFM2015
from BilevelOptimization.bilevel_solver import solve_bilevel_optimization
from BoundednessAnalysis.func_findNDShifting import func_findNDShifting

# Define a simple quadratic system
name = "SimpleBilevelSystem"

model = model_TwoState_SNJFM2015()

# Initial guess for m using func_findNDShifting
initial_m, _ = func_findNDShifting(model)

# If func_findNDShifting returns None or an invalid m, provide a fallback or raise an error
if initial_m is None:
    print("Warning: func_findNDShifting did not return a valid initial m. Using a default initial guess.")
    initial_m = np.array([[0.0], [5.0]])

print("Starting bilevel optimization...")
optimal_m, min_rTrap = solve_bilevel_optimization(model, initial_m,
                                                  learning_rate=0.1,
                                                  num_iterations=50,
                                                  verbose=True)

print("\nBilevel Optimization Finished.")
print(f"Optimal m: {optimal_m.flatten()}")
print(f"Minimum rTrap: {min_rTrap:.6f}")
