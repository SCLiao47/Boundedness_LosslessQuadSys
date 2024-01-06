init;

%% setups

% ifAlign = false;

caseNum = 2;
model = model_TwoState_Unbounded(caseNum);


% options
option.round_Ndigit = 3;
option.tol = 1e-6;
option.verbose = true;

%% check TR existence of the model
[isBounded, info] = func_boundedness_SDP(model, option);

astar = info.info_mND.a;
assert(astar >= 0, "ERROR: There is a TR and the system is bounded");

bstar = info.info_mPD.b;
assert(bstar >= 0, "ERROR: There is a TR and the system is bounded"); 

%% unbounded study
mPD = info.mPD;
model_PDshifted = func_ShiftSystem(model, mPD);

delta_SDP = info.delta_SDP;     % radius of escape region

% energy-decresing ellispoid
As = model_PDshifted.As;
d = model_PDshifted.d;

[U, E] = eig(As);
lamb = diag(E);

% center of ellipsoid at the standard coordinate
vE0 = -0.5*U'/As*d;

% semi-axis lengths
temp = d'/As*d;
alpha = 1/2 * sqrt(temp./lamb);

%% Trajectory
half_length = delta_SDP*4;
xbound = half_length*[-1, 1];
ybound = half_length*[-1, 1];

Param = struct('xbound', xbound-mPD(1), 'ybound', ybound-mPD(2));
odeOptions = odeset('Event', @(t,x) ode_EventFunc_sim(t,x,Param));

if_sim = ~true;

if if_sim
    num_traj = 100;
    radius = delta_SDP*1.5;
    tspan = [0, 3];

    x0s = 2*(rand(model.nx, num_traj)-0.5)*radius + mPD;
    Trajs = cell(num_traj, 1);
    
    for i = 1:num_traj
        [t, traj] = ode45(model.ode, tspan, x0s(:,i), odeOptions);

        Em = vecnorm((traj'-m))';
        Trajs{i} = struct('t', t, 'x', traj, 'Em', Em);
    end
    
    save("Data/TS_Unbounded_traj", "Trajs", "num_traj", "tspan");
else
    load("Data/TS_Unbounded_traj");
end

%%
Plot = figure(1); clf
Plot.WindowStyle='normal';
Plot.DockControls='off';

width = 700;
% width=700*2+10;
height=700;
set(gcf,'units','pixels','position',[1 1 width height]);
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

ax1=nexttile([1,1]);

% plot coordinate
limsize = half_length;
xlimit = [-limsize, limsize] - mPD(1);
ylimit = [-limsize, limsize] - mPD(2);
plot(xlimit, [0,0], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
hold on;
plot([0,0], ylimit, 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
 
% plot trajectory 
for i = 1:num_traj
    t_NAl = Trajs{i}.t;
    x = Trajs{i}.x;

    plot(x(:,1),x(:,2), 'color', [0.5,0.5,0.5]);
    hold on;
    scatter(x(1,1), x(1,2), 20, 'filled', 'MarkerFaceColor', [0.5,0.5,0.5]);
end

% [Escape Region]
th = linspace(0,2*pi, 200);
X = cos(th);
Y = sin(th);

% plot coordinate shift mPD
ph_mPD = scatter(mPD(1), mPD(2), 100, 'filled');
set(ph_mPD, 'CData', [1 0 1], 'SizeData', 200);

% plot Erengy growing region
E_zcord = [X*alpha(1) + vE0(1); Y*alpha(2) + vE0(2)];   % form the ellipsoid in standard coordinate
E_xcord = U*E_zcord;    % project into x-coordinate using eigenvectors U and shift m
ph_E0 = patch(E_xcord(1,:) + mPD(1), E_xcord(2,:) + mPD(2), 'g');
set(ph_E0, 'FaceAlpha', 0.2, 'LineStyle', '--', 'EdgeColor', 'g', 'LineWidth', 2);

% plot Escape Region
ph_ER = plot(X*delta_SDP + mPD(1), Y*delta_SDP + mPD(2));
set(ph_ER, 'Color', 'r', 'LineWidth', 2);

% plot critical point
ystar = info.info_ER.ystar + mPD;
ph_ystar = scatter(ystar(1,:), ystar(2,:), 100, 'filled', 'square');
set(ph_ystar, 'CData', [0.4940 0.1840 0.5560], 'SizeData', 200);

% formatting
axis equal;
grid on;
xlim(xlimit);
ylim(ylimit);
