import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import eigs
# from pathlib import Path

from Model.model_Lorenz import model_Lorenz
from BoundednessAnalysis import func_findNDShifting
from BoundednessAnalysis import func_TRSize
# from BoundednessAnalysis.func_ShiftSystem import func_ShiftSystem  # You'll need to convert this function too

def main():
    # Setups
    model = model_Lorenz()
    
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

    # Computing the ellipsoid Edot = 0
    lamb = np.diag(model_shifted.As)
    d = model_shifted.d.flatten()
    
    temp = np.sum(d**2/lamb)
    alpha = 1/2 * np.sqrt(temp/lamb)
    cE0 = 1/2 * np.divide(d, lamb)

    # Trajectories simulation
    if_sim = False
    
    if if_sim:
        num_traj = 100
        radius = r_SN * 2
        tspan = [0, 5]
        dt = 0.01  # Time resolution (seconds)
        num_points = int((tspan[1] - tspan[0]) / dt) + 1
        
        x0s = (np.random.rand(model.nx, num_traj) - 0.5) * radius
        Trajs = []
        
        for i in range(num_traj):
            # Define time points for solution
            t_eval = np.linspace(tspan[0], tspan[1], num_points)
            
            # Solve ODE using RK45 (equivalent to MATLAB's ode45)
            sol = solve_ivp(
                model_shifted.ode,
                t_span=tspan,
                y0=x0s[:,i].flatten(),
                method='RK45',
                t_eval=t_eval,
                rtol=1e-6,
                atol=1e-9
            )
            
            t = sol.t
            traj = sol.y.T  # Transpose to get shape (time_points, state_dim)
            Em = np.linalg.norm(traj - m.flatten(), axis=1)
            
            Trajs.append({
                't': t,
                'x': traj,
                'Em': Em
            })
            
        np.save("Data/Lorenze_traj.npy", {
            'Trajs': Trajs,
            'num_traj': num_traj,
            'tspan': tspan
        })
    else:
        data = np.load("Data/Lorenze_traj.npy", allow_pickle=True).item()
        Trajs = data['Trajs']
        num_traj = data['num_traj']
        tspan = data['tspan']
    
    # Visualization
    Tmax = 5
    
    fig = plt.figure(figsize=(10, 5.7))
    gs = fig.add_gridspec(1, 7)
    
    # State-space plot
    ax1 = fig.add_subplot(gs[0, :4])
    xlimits = [-120, 120]
    ylimits = [-80, 150]
    
    # Plot coordinate axes
    ax1.plot([xlimits[0]-10, xlimits[1]+10], [0, 0], 'k', linewidth=2)
    ax1.plot([0, 0], [ylimits[0]-10, ylimits[1]+10], 'k', linewidth=2)
    
    # Plot trajectories
    for i in range(10):
        traj = Trajs[i]
        t = traj['t']
        x = traj['x']
        x0 = x[0]
        
        idx = t <= Tmax
        ax1.plot(x[idx,1], x[idx,2], color='gray')
        ax1.scatter(x0[1], x0[2], s=20, c='gray')
    
    # Plot regions
    theta = np.linspace(0, 2*np.pi, 200)
    Y = np.cos(theta)
    Z = np.sin(theta)
    
    # Trapping regions
    ax1.plot(Y*r + m[1], Z*r + m[2], 'r--', linewidth=3, label='$B(m,R_m^*)$')
    ax1.plot(Y*r_SN + m[1], Z*r_SN + m[2], 'b--', linewidth=3, label='$B(m,R_m)$')
    
    # Ellipsoid
    ax1.fill(Y*alpha[1] + m[1,0] - cE0[1], 
             Z*alpha[2] + m[2,0] - cE0[2],
             color='g', alpha=0.2, label='$E$')
    
    # Critical points
    ystar = info_TR['ystar'] + m
    print(f"ystar = {ystar}")
    print(f"ystar shape = {ystar.shape}")
    ax1.scatter(ystar[:,1,0], ystar[:,2,0], s=100, c='purple', marker='D', label='$x^*$')
    
    # Formatting
    ax1.set_xlim(xlimits)
    ax1.set_ylim(ylimits)
    ax1.set_aspect('equal')
    ax1.grid(True)
    ax1.set_xlabel('$x_2$')
    ax1.set_ylabel('$x_3$')
    ax1.legend()
    
    # Energy plot
    ax2 = fig.add_subplot(gs[0, 4:])
    
    for i in range(10):
        traj = Trajs[i]
        t = traj['t']
        Em = traj['Em']
        
        idx = t <= Tmax
        ax2.plot(t[idx], Em[idx], color='gray')
    
    ax2.plot([0, Tmax], [r, r], 'r--', linewidth=2)
    ax2.plot([0, Tmax], [r_SN, r_SN], 'b--', linewidth=2)
    
    ax2.grid(True)
    ax2.set_xlabel('Time (sec)')
    ax2.set_ylabel('$K_m(x(t))$')
    
    plt.tight_layout()
    plt.show()

    # # Save figures
    # plt.savefig('Figure/Lorenze_TR_2D.pdf', bbox_inches='tight')
    # plt.savefig('Figure/Lorenze_TR_2D.png', bbox_inches='tight', dpi=300)
    
if __name__ == "__main__":
    main() 