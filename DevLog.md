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

# Todo next step
- [x] implement 2d example
- [ ] differentiable SDP for TR size w.r.t. coordinate shift vector m
- [ ] gradient descent on m to find the optimal coordinate shift for mimizing TR size
- [ ] verify the size on a grid of coordinate shift
- [ ] compare the result from GD and actual grid of size

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