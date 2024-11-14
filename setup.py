from setuptools import setup, find_packages

setup(
    name="boundedness_lossless_quad_sys",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'cvxpy',
        'mosek',
    ],
    python_requires='>=3.6',
) 