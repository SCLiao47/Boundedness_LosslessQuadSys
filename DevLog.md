This is a log of the development process of the code.

# Setting Up Python Development Environment
To set up the development environment for this project, follow these steps:

1. **Ensure directory structure:**
   - Make sure your directory structure is correct and includes an `__init__.py` file in each package directory to ensure they are recognized as Python packages.

2. **Install the package in development mode:**
   - Run the following command from the root directory of the project:
     ```bash
     pip install -e .
     ```
   - This installs the package in "editable" mode, allowing you to make changes to the source code without needing to reinstall the package.

3. **Update .gitignore:**
   - Add `*.egg-info/` to your `.gitignore` to avoid committing build artifacts to version control.


# Finding the Smallest Trapping Region using differentiability of SDP

In the previous paper, we proposed a SDP-based method to find the tightest trapping region given a coordinate shift. It would be interesting to find the smallest trapping region across all possible coordinate shifts. In this section, we will propose an algoroithm to find the smallest trapping region by leveraging the differentiability of the SDP solution. The algorithm will apply the idea of gradient descent with multiple starting points to find the smallest trapping region.

## Development plan
- [x] implement 2d example
- [ ] get the gradient of the SDP solution w.r.t. the coordinate shift
- [ ] use torch to implement the multiple starting points gradient descent. [Sample](https://mail.google.com/mail/u/2/#inbox/FMfcgzQXKhQSmvVMHQVVStrMFmQJjWhT)
- [ ] Checking the landscape of the SDP solution w.r.t. the coordinate shift. 
  - [ ] Test on the Two state system, validate by gridding the coordinate shift space. 
  - [ ] Test on the Lorenz system, plot the TRSize-vs-iteration curve for multiple starting points.

## Bugs to fix
- [x] the plotting function would create figure when not wanted


## References
1. [Derivative of CVXPy](https://www.cvxpy.org/examples/index.html#derivatives)
   1. [probelm.backward() API](https://www.cvxpy.org/api_reference/cvxpy.problems.html#cvxpy.Problem.backward)
2. [cvxpylayers](https://github.com/cvxgrp/cvxpylayers): [Blog post](https://locuslab.github.io/2019-10-28-cvxpylayers/)

# Log

## 2024-12-17
- add the script `script_TwoState_SNJFM2015.py` to analysis, simulation and visualization of the two-state model
- fix the dimension of `ystar` in the function `func_TRSize`


## 2024-11-14
- add the function `func_TRSize` to compute the size of TR using both SDP and SN method
- edit the script `script_Lorenz.py` and validate the result with matlab implementation

## 2024-11-12 Initialize python branch
- Added script_Lorenz.py to test the functionality of the code.
- The script_Lorenz.py script is a simple script to test the functionality of the code.
- The script_Lorenz.py script is a simple script to test the functionality of the code.