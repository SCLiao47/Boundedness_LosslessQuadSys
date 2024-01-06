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
m = info.mND;
astar = info.info_mND.a;

assert(astar >= 0, "There is a TR and the system is bounded");

model_shifted = func_ShiftSystem(model, m);

%% dual variable and directions
W = info.info_mND.W;

[V,D] = eig(W);

d = diag(D);
idx = find(d>0);

% energy growing direction with scaling by eigenvalues
V = V(:,idx)*diag(sqrt(d(idx)));

% solve for scaling
% s = norm(model_shifted.d) ./ astar;
ss = abs(model_shifted.d'*V) ./ astar;

s = ss;

% compute energy growing states
if model_shifted.d'*V >= 0    % constant dynamics aligned, energy grows more!
    xstar = m + s*V;
else
    xstar = m - s*V;
end

%% check rank 1
assert(size(V,2) == 1);

%% plot energy_dot growing direction for rank(W)==1 case
figure(1); clf;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

ax1 = nexttile([3, 3]);

%% Visualize energy on the circle
axes(ax1); cla;

half_length = 20;
xbound = half_length*[-1, 1];
ybound = half_length*[-1, 1];

plot(xbound, [0, 0], 'k', 'linewidth', 2); hold on;
plot([0, 0], ybound, 'k', 'linewidth', 2);

y_Vdir_atm = @(x) V(2)/V(1)*(x - m(1)) + m(2);
y_Vdir_at0 = @(x) V(2)/V(1)*x;
ph_V = plot(xbound, y_Vdir_atm(xbound), 'm--', 'linewidth', 2);
% ph_V = plot(xbound, y_Vdir_at0(xbound), 'g--', 'linewidth', 2);

% [Null line]
% model.La = model.L - model.Ls;
y_NLine = @(x) -1/model.Q(2,2,1)*(2*model.Q(1,2,1)*x + model.La(1,2));
ph_NL = plot(xbound, y_NLine(xbound), 'g--', 'linewidth', 2);

plot_TwoState_quiver(model.ode, xbound, ybound);
ph_circle = plot_TwoState_Kmdot_onCircle(model, V, s, m);

set_axes('Unbounded 2-State', xbound, ybound)

%% Trajectory along s*V
axes(ax1);

Param = struct('xbound', xbound-m(1), 'ybound', ybound-m(2));
odeOptions = odeset('Event', @(t,x) ode_EventFunc_sim(t,x,Param));

scale = 2;
Tinf = 10;

[tp, xp] = ode45(model_shifted.ode_shifted, [0,Tinf], 1*s*V, odeOptions);
[tsp, xsp] = ode45(model_shifted.ode_shifted, [0,Tinf], scale*s*V, odeOptions);

[tn, xn] = ode45(model_shifted.ode_shifted, [0,Tinf], -s*V, odeOptions);
[tsn, xsn] = ode45(model_shifted.ode_shifted, [0,Tinf], -scale*s*V, odeOptions);


ph_xp = plot(m(1)+xp(:,1), m(2)+xp(:,2), 'r', 'linewidth', 2);
ph_xn = plot(m(1)+xn(:,1), m(2)+xn(:,2), 'b', 'linewidth', 2);

% ph_xsp = plot(m(1)+xsp(:,1), m(2)+xsp(:,2), 'r--', 'linewidth', 2);
% ph_xsn = plot(m(1)+xsn(:,1), m(2)+xsn(:,2), 'b--', 'linewidth', 2);

legend([ph_circle.center, ph_V, ph_NL, ph_xp, ph_xn], ...
        '$m^*$', 'Dual', 'Rot-Null', '$x_p$', '$x_n$', ...
    'Interpreter','latex', 'fontsize', 20);

%% plot constant, linear, quadratic dynamics seperately
ode_const_shifted = @(t,x) model_shifted.d;
ode_linear_shifted = @(t,x) model_shifted.A*x;
ode_quad_shifted = @(t,x) [x'*model_shifted.Q(:,:,1)*x; x'*model_shifted.Q(:,:,2)*x];

figure(2); clf;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax2 = nexttile([1, 1]);
ax3 = nexttile([1, 1]);
ax4 = nexttile([1, 1]);

% constant dynamics
axes(ax2); cla
plot(xbound, y_Vdir_atm(xbound), 'm--'); hold on;
plot_TwoState_quiver(ode_const_shifted, xbound, ybound);

set_axes('Constant', xbound, ybound)


% linear dynamics
axes(ax3); cla
plot(xbound, y_Vdir_atm(xbound), 'm--'); hold on;
plot_TwoState_quiver(ode_linear_shifted, xbound, ybound);

set_axes('Linear', xbound, ybound)


% quadratic dynamics
axes(ax4); cla
plot(xbound, y_Vdir_atm(xbound), 'm--'); hold on;
plot_TwoState_quiver(ode_quad_shifted, xbound, ybound);

set_axes('Quadratic', xbound, ybound)


%%
figure(3); clf;
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax31 = nexttile([1, 1]);
ax32 = nexttile(4, [1, 1]);
ax33 = nexttile([2, 2]);

% model.La = model.L - model.Ls;

ode_linear_asym = @(t,x) model.La*x;
ode_quad = @(t,x) [x'*model.Q(:,:,1)*x; x'*model.Q(:,:,2)*x];

y_NLine_at0 = @(x) -1/model.Q(2,2,1)*2*model.Q(1,2,1)*x;
y_NLine_wLa = @(x) y_NLine_at0(x) - 1/model.Q(2,2,1)*model.La(1,2);

axes(ax31); cla; 
plot_TwoState_quiver(ode_linear_asym, xbound, ybound);
set_axes('Asymmetric linear dynamics', xbound, ybound)

axes(ax32); cla;
plot_TwoState_quiver(ode_quad, xbound, ybound); hold on;
plot(xbound, y_NLine_at0(xbound), '--', 'linewidth', 2);
set_axes('Lossless quadratic dynamics', xbound, ybound)

axes(ax33); cla;
plot_TwoState_quiver(@(t,x) ode_quad(t,x) + ode_linear_asym(t,x), ...
                     xbound, ybound); hold on;

ph_NL = plot(xbound, y_NLine_wLa(xbound), 'g--', 'linewidth', 2);

set_axes('Rotation field (AsymLinear + Quad)', xbound, ybound)


%%
model.La = model.L - model.Ls;

ode_c = @(t,x) model.c;
ode_Lsx = @(t,x) model.Ls*x;
ode_Lax = @(t,x) model.La*x;
ode_Lx = @(t,x) ode_Lsx(t,x) + ode_Lax(t,x);
ode_xQx = @(t,x) [x'*model.Q(:,:,1)*x; x'*model.Q(:,:,2)*x];

ode_full = @(t,x) ode_c(t,x) + ode_Lx(t,x) + ode_xQx(t,x);
ode_roatation = @(t,x) ode_Lax(t,x) + ode_xQx(t,x);
ode_SDP = @(t,x) ode_Lsx(t,x) + ode_xQx(t,x);

figure(4); clf;
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
ax4_UL = nexttile([3, 3]);
ax4_R1 = nexttile([1, 1]);
ax4_R2 = nexttile([1, 1]);
ax4_R3 = nexttile([1, 1]);
ax4_D1 = nexttile([1, 1]);
ax4_D2 = nexttile([1, 1]);
ax4_D3 = nexttile([1, 1]);
ax4_D4 = nexttile([1, 1]);

% full dyanmics
axes(ax4_UL); cla; 
plot_TwoState_quiver(ode_full, xbound, ybound); hold on;
set_axes('Full', xbound, ybound);

% constant dynamics
axes(ax4_R1); cla;
plot_TwoState_quiver(ode_c, xbound, ybound); hold on;
set_axes('c', xbound, ybound);

% linear dynamics
axes(ax4_R2); cla;
plot_TwoState_quiver(ode_Lx, xbound, ybound); hold on;
set_axes('(Ls + La)*x', xbound, ybound);

% quadratic dynamics
axes(ax4_R3); cla;
plot_TwoState_quiver(ode_xQx, xbound, ybound); hold on;
set_axes('f(x)', xbound, ybound);

% Linear symmetric
axes(ax4_D1); cla;
plot_TwoState_quiver(ode_Lsx, xbound, ybound); hold on;
set_axes('Ls*x', xbound, ybound);

% Linear assymmetric
axes(ax4_D2); cla;
plot_TwoState_quiver(ode_Lax, xbound, ybound); hold on;
set_axes('La*x', xbound, ybound);

% Rotational dynamics
axes(ax4_D3); cla;
plot_TwoState_quiver(ode_roatation, xbound, ybound); hold on;
set_axes('La*x + f(x)', xbound, ybound);

% SDP analysis dynamics
axes(ax4_D4); cla;
plot_TwoState_quiver(ode_SDP, xbound, ybound); hold on;
set_axes('Ls*x + f(x)', xbound, ybound);

%% utility
function set_axes(title_name, xbound, ybound)
    title(title_name, 'fontsize', 20);
    grid on;
    axis equal;
    xlim(xbound);
    ylim(ybound);
end