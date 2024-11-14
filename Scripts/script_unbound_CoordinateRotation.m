

c = zeros(2,1);

%
Ls = diag([-1,2]);
La = [0 1; -1, 0];

L = Ls + 2*La;
%
Q = zeros(2,2,2);
Q(:,:,1) = [0, 1; 1, 2];
Q(:,:,2) = [-2, -1; -1, 0];

%
model = class_Model_LosslessQuad([], c, L, Q);

%%
% options
option.round_Ndigit = 3;
option.tol = 1e-6;
option.verbose = false;

[isBounded, info] = func_boundedness_SDP(model, option);
m = info.mND;
astar = info.info_mND.a;
W = info.info_mND.W;

%%
% Rotation

[V,D] = eig(W);

% swap the indices s.t. lagerest e-value is the x1 direction
[~, I] = sort(diag(D), 'descend');
D = D(I,I);
V = V(:,I);    

U = V';     % transpose it to match the notation in notes

%
model_rot = func_RotateSystem(model, U);

%% Simulation
x0_a = [5;0];
x0_b = [0;5];

func_plot(model, 1, U'*x0_a);
func_plot(model_rot, 2, x0_a);


%% unbounded direction
% ZBAR
M = 2*model_rot.Q(1,2,2);
zbar = - model_rot.L(2,1) / M;

linebar = [-20, 20; zbar, zbar];
line = U'*linebar;

figure(2);
plot(linebar(1,:), linebar(2,:), '--', 'linewidth', 1.5);

figure(1);
plot(line(1,:), line(2,:), '--', 'linewidth', 1.5);

%%
function func_plot(model, fig_num, x0_a, x0_b)
    figure(fig_num); clf

    half_length = 20;
    xbound = half_length*[-1, 1];
    ybound = half_length*[-1, 1];

    plot(xbound, [0, 0], 'k', 'linewidth', 2); hold on;
    plot([0, 0], ybound, 'k', 'linewidth', 2);

    plot_TwoState_quiver(model.ode, xbound, ybound)

    % simulation
    Param = struct('xbound', xbound, 'ybound', ybound);
    odeOptions = odeset('Event', @(t,x) ode_EventFunc_sim(t,x,Param));

    Tinf = 10;

    [tp_a, xp_a] = ode45(model.ode, [0,Tinf], 1*x0_a, odeOptions);
    [tn_a, xn_a] = ode45(model.ode, [0,Tinf], -x0_a, odeOptions);

    plot(xp_a(:,1), xp_a(:,2), 'r', 'linewidth', 2)
    plot(xn_a(:,1), xn_a(:,2), 'b', 'linewidth', 2);
    
    title(model.name);
    grid on;
    axis equal;
    xlim(xbound);
    ylim(ybound);
    
%     disp(xp_a(end,2));
%     disp(xn_a(end,2));
end