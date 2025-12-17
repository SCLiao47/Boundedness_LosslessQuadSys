import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from Model.model_TwoState_Alignment import model_TwoState_Alignment
from BoundednessAnalysis import func_findNDShifting
from BoundednessAnalysis import func_TRSize

def main():
    # Setups
    ifAlign = False
    model_NAl = model_TwoState_Alignment(ifAlign)

    # options
    option = {
        'round_Ndigit': 3,
        'tol': 1e-6,
        'verbose': True
    }

    # verify the model is bounded
    eig_values_NAL = np.linalg.eigvals(model_NAl.Ls)
    assert np.all(eig_values_NAL < 0), 'Some eigenvalue of Ls is not negative. Check model!'

    # compute the radius of TR
    # Shift the model by 0 vector to append elements {d, A, As} to model.
    m = np.array([[0], [0]])
    model_NAl.func_ShiftSystem(m)

    # SDP analysis (proposed)
    r_NAl, info_TR_NAl = func_TRSize(model_NAl, method='SDP', option=option)

    # spectrum analysis (Schlegel and Noack)
    r_SN_NAl = func_TRSize(model_NAl, method='SN')

    print(f"SDP: r = {r_NAl}, SN: r = {r_SN_NAl}")

if __name__ == "__main__":
    main()
