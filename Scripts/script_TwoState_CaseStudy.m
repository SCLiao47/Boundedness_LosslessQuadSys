init;

%% Setup

CaseNum = 5;
mag_La = 0.2;

% SDP options
option.round_Ndigit = 3;
option.tol = 1e-6;
option.verbose = false;

% Simulation options
if_sim = true;
half_length = 3;
mag_x0 = 1;

%% initialize model
model = model_TwoState_CaseStudy(CaseNum, mag_La);

fprintf('=== [Running Case%i] === \n', CaseNum);

disp('Ls = ');
disp(model.Ls);
disp('[Q(:,:,1), Q(:,:,2)] = ');
disp([model.Q(:,:,1), model.Q(:,:,2)]); 

%% check TR existence of the model

[m_ND, info_ND] = func_findNDShifting(model, option);
[m_PD, info_PD] = func_findPDShifting(model, option);

% 
% [isBounded, info] = func_boundedness_SDP(model, option);

astar = info_ND.a;
bstar = info_PD.b;

assert(abs(astar - bstar) < 1e-8, ...
       "[ERROR] astar ~= bstar, Need to check the implementation"); 

disp('[a*, b*] = ');
disp([astar, bstar]);

%% Simulation
xbound = half_length*[-1, 1];
ybound = half_length*[-1, 1];

Param = struct('xbound', xbound, 'ybound', ybound);
odeOptions = odeset('Event', @(t,x) ode_EventFunc_sim(t,x,Param));

x0s = [1,  1, -1, -1, 1, -1, 0,  0;
       1, -1,  1, -1, 0,  0, 1, -1];
x0s = mag_x0*x0s;

if if_sim
    tspan = [0, 10];

    num_traj = size(x0s, 2);
    Trajs = cell(num_traj, 1);
    
    for i = 1:num_traj
        [t, traj] = ode45(model.ode, tspan, x0s(:,i), odeOptions);

        Em = vecnorm((traj'))';
        Trajs{i} = struct('t', t, 'x', traj, 'Em', Em);
    end
    
%     save("Data/", "Trajs", "num_traj", "tspan");
else
%     load("Data/");
end

%% Plotting
fig1 = figure(1); clf;
fig1.WindowStyle = 'normal';
fig1.DockControls = 'off';

width = 700;
height = 700;
set(fig1,'units','pixels','position',[1 1 width height]);
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile([1,1]);

% plot axis
plot(xbound, [0, 0], 'k', 'linewidth', 2); hold on;
plot([0, 0], ybound, 'k', 'linewidth', 2);

% plot quiver
plot_TwoState_quiver(model.ode, xbound, ybound, false, 1);
hold on;

% plot trajectory 
for i = 1:num_traj
    t_NAl = Trajs{i}.t;
    x = Trajs{i}.x;

    plot(x(:,1),x(:,2), 'color', [0.5,0.5,0.5], 'linewidth', 2);
    hold on;
    scatter(x(1,1), x(1,2), 30, 'filled', 'MarkerFaceColor', 0.5*[1,1,1]);
end

xlim(xbound);
ylim(ybound);