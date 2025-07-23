import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from Model.model_TwoState_SNJFM2015 import model_TwoState_SNJFM2015
from BoundednessAnalysis import func_findNDShifting
from BoundednessAnalysis import func_TRSize

@dataclass
class Trajectory:
    """
    Data structure to hold simulation trajectory data
    
    Attributes:
    -----------
    t : ndarray
        Time points array
    solution : ndarray
        Solution array with shape (n_timesteps, n_states)
    """
    t: np.ndarray
    solution: np.ndarray

def simulate(ode, t_max, initial_conditions=None, dt=0.01):
    """
    Simulate the two-state system evolution
    
    Returns:
    --------
    trajectory : Trajectory
        Object containing the time points and solution arrays
    """
    if initial_conditions is None:
        initial_conditions = [1.0, 0.0]  # Start with all population in state 1
        
    t = np.linspace(0, t_max, int(t_max/dt))
    sol = solve_ivp(ode, [t[0], t[-1]], initial_conditions, t_eval=t)
    solution = sol.y.T

    return Trajectory(t=t, solution=solution)

def simulate_multiple(ode, t_max, initial_conditions_list, dt=0.01):
    """
    Simulate the system evolution from multiple initial conditions
    
    Parameters:
    -----------
    ode : callable
        The ODE function to simulate
    t_max : float
        Maximum simulation time
    initial_conditions_list : list
        List of initial condition arrays
    dt : float, optional
        Time step for integration
        
    Returns:
    --------
    list[Trajectory]
        List of trajectory objects for each initial condition
    """
    return [simulate(ode, t_max, ic, dt) for ic in initial_conditions_list]

def plot_phase_portrait_multiple(trajectories, ax=None, colors=None, show=False):
    """
    Plot the phase portrait of the two-state system
    
    Parameters:
    -----------
    trajectories : list[Trajectory]
        List of trajectory objects
    ax : matplotlib.axes.Axes, optional
        The axes object to plot the phase portraits. If None, creates new figure
    colors : list, optional
        List of colors for each trajectory
    show : bool, optional
        Whether to show the plot immediately. Default is False
    
    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axes object containing the phase portraits
    """
    if ax is None:
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
    
    if colors is None:
        colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
    
    # Plot each trajectory with its color and direction arrow
    for traj, color in zip(trajectories, colors):
        # Plot trajectory
        ax.plot(traj.solution[:, 0], traj.solution[:, 1], 
                color=color)
        
        # Add arrow to show direction
        arrow_idx = len(traj.t) // 4  # Add arrow at 1/4 of the trajectory
        ax.arrow(traj.solution[arrow_idx, 0], traj.solution[arrow_idx, 1],
                traj.solution[arrow_idx + 1, 0] - traj.solution[arrow_idx, 0],
                traj.solution[arrow_idx + 1, 1] - traj.solution[arrow_idx, 1],
                head_width=0.02, head_length=0.03, fc=color, ec=color)
    
    ax.set_xlabel('State 1 Population')
    ax.set_ylabel('State 2 Population')
    ax.set_title('Phase Portrait')
    ax.grid(True)
    ax.set_aspect('equal')  # Make the plot square with equal axes
    
    if show:
        plt.tight_layout()
        plt.show()
    
    return ax

def plot_time_evolution(trajectories, ax=None, colors=None, show=False):
    """
    Plot the time evolution of states for multiple trajectories
    
    Parameters:
    -----------
    trajectories : list[Trajectory]
        List of trajectory objects to plot
    ax : matplotlib.axes.Axes, optional
        The axes object to plot on. If None, creates new figure
    colors : list, optional
        List of colors for different trajectories
    show : bool, optional
        Whether to show the plot immediately. Default is False
    
    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axes object containing the time evolution plot
    """
    if ax is None:
        plt.figure(figsize=(10, 6))
        ax = plt.gca()
    
    if colors is None:
        colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
    
    # Plot time evolution for each trajectory
    for traj, color in zip(trajectories, colors):
        ax.plot(traj.t, traj.solution[:, 0], '--', color=color)
        ax.plot(traj.t, traj.solution[:, 1], '-', color=color)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax.set_title('Two-State System Evolution')
    ax.grid(True)
    
    if show:
        plt.tight_layout()
        plt.show()
    
    return ax

def plot_results_multiple(trajectories, show=False):
    """
    Plot the population evolution and phase portraits for multiple trajectories
    
    Parameters:
    -----------
    trajectories : list[Trajectory]
        List of trajectory objects to plot
    show : bool, optional
        Whether to show the plot immediately. Default is False
    
    Returns:
    --------
    tuple
        (fig, (ax1, ax2)) containing the figure and axes objects
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
    
    # Time evolution plot
    plot_time_evolution(trajectories, ax=ax1, colors=colors)
    
    # Phase portrait
    plot_phase_portrait_multiple(trajectories, ax=ax2, colors=colors)
    
    plt.tight_layout()
    
    if show:
        plt.show()
    
    return fig, (ax1, ax2)

def generate_initial_conditions(x_span=(-5, 5), y_span=(-5, 5), n_points=5):
    """
    Generate a grid of initial conditions spanning the specified ranges
    
    Parameters:
    -----------
    x_span : tuple, optional
        (min, max) values for x-axis (State 1 population)
        Default is (-5, 5)
    y_span : tuple, optional
        (min, max) values for y-axis (State 2 population)
        Default is (-5, 5)
    n_points : int, optional
        Number of points along each axis
        Default is 5 (resulting in 25 total points)
    
    Returns:
    --------
    list
        List of initial conditions [x, y] spanning the grid
    """
    # Create grid points
    x = np.linspace(x_span[0], x_span[1], n_points)
    y = np.linspace(y_span[0], y_span[1], n_points)
    X, Y = np.meshgrid(x, y)
    
    # Convert to list of points
    initial_conditions = [[x, y] for x, y in zip(X.flatten(), Y.flatten())]
    
    return initial_conditions

def plot_boundedness_regions(m, r, r_SN, ystar=None, ax=None):
    """
    Plot boundedness regions and critical points
    
    Parameters:
    -----------
    m : ndarray
        Shifting point
    r : float
        Radius from SDP method
    r_SN : float
        Radius from SN method
    ystar : ndarray, optional
        Critical points. If provided, should be shape (n_points, nx, 1)
    ax : matplotlib.axes.Axes, optional
        The axes object to plot on. If None, creates new figure
    show : bool, optional
        Whether to show the plot immediately. Default is False
    
    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axes object containing the plot
    """
    if ax is None:
        plt.figure(figsize=(8, 8))
        ax = plt.gca()

    # Plot coordinate axes only for new figure
    ax.axhline(y=0, color='k', linewidth=2)
    ax.axvline(x=0, color='k', linewidth=2)
    
    # Plot the trapping regions
    theta = np.linspace(0, 2*np.pi, 200)
    X = np.cos(theta)
    Y = np.sin(theta)
    
    # Plot SDP and SN bounds
    ax.plot(X*r + m[0], Y*r + m[1], 'r--', linewidth=3, label='$B(m,R_m^*)$ (SDP)')
    ax.plot(X*r_SN + m[0], Y*r_SN + m[1], 'b--', linewidth=3, label='$B(m,R_m)$ (SN)')
    
    # Plot critical points if available
    if ystar is not None:
        if ystar.ndim == 1:  # If ystar is 1D array
            # Reshape to (1, nx, 1)
            ystar = ystar.reshape(1, -1, 1)
        ax.scatter(m[0]+ystar[:,0,0], m[1]+ystar[:,1,0], s=100, c='purple', marker='D', label='$x^*$')
    
    ax.set_xlabel('State 1 Population')
    ax.set_ylabel('State 2 Population')
    ax.grid(True)
    ax.set_aspect('equal')
    ax.legend()
    
    return ax

def analyze_boundedness(model, option=None, show=False):
    """
    Analyze the boundedness of the two-state system
    
    Parameters:
    -----------
    model : model_TwoState_SNJFM2015
        The model to analyze
    option : dict, optional
        Options for the analysis
    show : bool, optional
        Whether to show plots immediately. Default is False
    """
    if option is None:
        option = {
            'round_Ndigit': 3,
            'verbose': True,
            'tol': 1e-6
        }
    
    # Verify the existence of boundedness region
    m, info_m = func_findNDShifting(model, option)
    model_shifted = model.func_ShiftSystem(m)
    
    # Display the results
    print("\nBoundedness Region Analysis Results:")
    print(f"Shifting point m: {m}")
    print("\nAdditional Information:")
    for key, value in info_m.items():
        print(f"{key}: {value}")
    
    # Solve for the size of boundedness region
    r, info_TR = func_TRSize(model_shifted, method='SDP', option=option)
    r_SN = func_TRSize(model_shifted, method='SN')
    
    print(f"SDP: r = {r}, SN: r = {r_SN}")
    
    # Get critical points if available
    ystar = info_TR.get('ystar', None)
    
    # Plot the boundedness regions
    if show:
        plot_boundedness_regions(m, r, r_SN, ystar, show=show)
    
    return m, r, r_SN, model_shifted, ystar

if __name__ == "__main__":
    '''
    Setup
    '''
    # Model parameters
    t_max = 10.0
    n_grid = 7  # 7x7 grid = 49 trajectories
    span = (-5, 5)  # Same span for both axes
    
    # Create model instance
    model = model_TwoState_SNJFM2015()
    
    '''
    Boundedness Analysis by Liao et al., 2024
    '''
    # Boundedness analysis
    m, r, r_SN, model_shifted, ystar = analyze_boundedness(model)
    

    '''
    Simulation
    '''
    # Generate and simulate trajectories
    initial_conditions = generate_initial_conditions(
        x_span=span,
        y_span=span,
        n_points=n_grid
    )
    trajectories = simulate_multiple(model.ode, t_max, initial_conditions)
    
    '''
    Visualization
    '''
    # Create figure
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    
    # Plot phase portrait
    ax = plot_phase_portrait_multiple(trajectories, ax=ax)
    
    # Add boundedness regions
    plot_boundedness_regions(m, r, r_SN, ystar, ax=ax)
    
    # # Show plot
    plt.show()


    '''
    Finding the smallest trapping region
    
    # TODO: 
    1. Write a SDP function to compute the size of trapping region for a given coordinate shift m
    
    2. Write a wrapper function 
    '''
    # 