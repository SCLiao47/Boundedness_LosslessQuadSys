init;

%% setups
ifAlign = false;
model_NAl = model_TwoState_Alignment(ifAlign);

% options
option.round_Ndigit = 3;
option.tol = 1e-6;
option.verbose = true;

%% verify the model is bounded
eig_values_NAL = eig(model_NAl.Ls);
assert(all(eig_values_NAL<0), 'Some eigenvalue of Ls is not negative. Check model!');

%% compute the radius of TR
% Shift the model by 0 vector to append elements {d, A, As} to model.
% This allows computing the size using existing functions.
m = [0;0];
model_NAl = func_ShiftSystem(model_NAl, m);    

% SDP analysis (proposed)
[r_NAl, info_TR_NAl] = func_TRSize_SDP(model_NAl, option);

% spectrum analysis (Schlegel and Noack)
r_SN_NAl = func_TRSize_SN(model_NAl, option);

%% Computing the ellipsoid Edot = 0
% In this problem, As is already a diagonal matrix. Hence no need to
% transform coordinate
lamb_NAl = diag(model_NAl.As);
d_NAl = model_NAl.d;

temp = sum(d_NAl.^2./lamb_NAl);
alpha_NAl = 1/2 * sqrt(temp./lamb_NAl);
cE0_NAl = 1/2 * d_NAl./lamb_NAl;

%% Trajectory
% trajectories
if_sim = ~true;

if if_sim
    num_traj = 100;
    radius = max(r_NAl,r_SN_NAl)*2.5;
    tspan = [0, 3];

    x0s = (rand(model_NAl.nx, num_traj)-0.5)*radius + m;
    Trajs_NAl = cell(num_traj, 1);
    
    for i = 1:num_traj
        [t_NAl, traj_NAl] = ode45(model_NAl.ode, tspan, x0s(:,i));

        Em_NAl = vecnorm((traj_NAl'-m))';
        Trajs_NAl{i} = struct('t', t_NAl, 'x', traj_NAl, 'Em', Em_NAl);
    end
    
    save("Data/TwoStateFig_traj", "Trajs_NAl", "num_traj", "tspan");
else
    load("Data/TwoStateFig_traj");
end

%% Preliminary Section plot
Plot_Prelim = figure(1); clf
Plot_Prelim.WindowStyle='normal';
Plot_Prelim.DockControls='off';

width = 700;
height=700;
set(gcf,'units','pixels','position',[1 1 width height]);
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

ax1=nexttile([1,1]);

% plot coordinate
limsize = 1.3;
plot([-limsize, limsize], [0,0], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
hold on;
plot([0,0], [-limsize, limsize], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
 
% plot trajectory 
for i = 1:num_traj
    t_NAl = Trajs_NAl{i}.t;
    x = Trajs_NAl{i}.x;

    plot(x(:,1),x(:,2), 'color', [0.5,0.5,0.5]);
    hold on;
    scatter(x(1,1), x(1,2), 20, 'filled', 'MarkerFaceColor', [0.5,0.5,0.5]);
end

% [Trapping region]
th = linspace(0,2*pi, 200);
X = cos(th);
Y = sin(th);

% [Energy growing region]
ph_E0 = patch(X*alpha_NAl(1) - cE0_NAl(1) + m(1), Y*alpha_NAl(2) - cE0_NAl(2) + m(2), 'g');
set(ph_E0, 'FaceAlpha', 0.2, 'LineStyle', '-', 'EdgeColor', 'g', 'LineWidth', 2);

% SDP estimate
ph_S1 = plot(X*r_NAl + m(1), Y*r_NAl + m(2));
set(ph_S1, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% SN estimate
ph_S2 = plot(X*r_SN_NAl + m(1), Y*r_SN_NAl + m(2));
set(ph_S2, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

% equilibrium point
ph_yeq = scatter(0, 0.25, 100, 'filled', 'square');
set(ph_yeq, 'CData', [1 0 1], 'SizeData', 200);

% formatting
formatSetting = {'interpreter', 'latex', 'fontsize', 20};

text(limsize-0.15, 0.1, '$x_1$', formatSetting{:})
text(0.05, limsize-0.1, '$x_2$', formatSetting{:})
legend([ph_E0, ph_yeq, ph_S2, ph_S1], ...
       {'$E$', 'Equilibrium', '$B(0,R_0)$', 'Proposed'}, ...
       formatSetting{:}, 'fontsize', 16);

formatting_PhasePlot(ax1, limsize);

print('Figure/TwoState_Prelim','-depsc')

%% Plot for example section
% [visualization 1, State-space plot]
Plot_Apdx = figure(2); clf
Plot_Apdx.WindowStyle='normal';
Plot_Apdx.DockControls='off';

width = 600;
height=600+200;
set(gcf,'units','pixels','position',[1 1 width height]);
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

ax1=nexttile([3,1]);
% plot coordinate
plot([-limsize, limsize], [0,0], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
hold on;
plot([0,0], [-limsize, limsize], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');

% plot trajectory 
for i = 1:num_traj
    t_NAl = Trajs_NAl{i}.t;
    x = Trajs_NAl{i}.x;

    plot(x(:,1),x(:,2), 'color', [0.5,0.5,0.5]);
    hold on;
    scatter(x(1,1), x(1,2), 20, 'filled', 'MarkerFaceColor', [0.5,0.5,0.5]);
end

% [Trapping region]
th = linspace(0,2*pi, 200);
X = cos(th);
Y = sin(th);

% SDP estimate
ph_S1 = plot(X*r_NAl + m(1), Y*r_NAl + m(2));
set(ph_S1, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% SN estimate
ph_S2 = plot(X*r_SN_NAl + m(1), Y*r_SN_NAl + m(2));
set(ph_S2, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

% [Energy growing region]
ph_E0 = patch(X*alpha_NAl(1) - cE0_NAl(1) + m(1), Y*alpha_NAl(2) - cE0_NAl(2) + m(2), 'g');
set(ph_E0, 'FaceAlpha', 0.2, 'LineStyle', '-', 'EdgeColor', 'g', 'LineWidth', 2);

% equilibrium point
ph_yeq = scatter(0, 0.25, 100, 'filled', 'square');
set(ph_yeq, 'CData', [1 0 1], 'SizeData', 200);

% [ciritcal points]
ystar = info_TR_NAl.ystar + m;
ph_ystar = scatter(ystar(1,:), ystar(2,:), 100, 'filled', 'Diamond');
set(ph_ystar, 'CData', [0.4940 0.1840 0.5560], 'SizeData', 200);

% formatting
formatSetting = {'interpreter', 'latex', 'fontsize', 20};

text(limsize-0.15, 0.1, '$x_1$', formatSetting{:})
text(0.05, limsize-0.1, '$x_2$', formatSetting{:})
legend([ph_E0, ph_yeq, ph_S2, ph_S1, ph_ystar], ...
       {'$E$', 'Equilibrium', '$B(0, R_0)$', '$B(0, R_0^*)$', '$x^*$'}, ...
       formatSetting{:}, 'fontsize', 16);

formatting_PhasePlot(ax1, limsize)


% [visualization 2, Energy plot]
ax2=nexttile([1,1]);

% trajectory
for i = 1:num_traj
    t_NAl = Trajs_NAl{i}.t;
    Em_NAl = Trajs_NAl{i}.Em;
    
    plot(t_NAl, Em_NAl, 'color', [0.5,0.5,0.5]);
    hold on;
end

plot(tspan, [r_NAl; r_NAl], 'r--', 'linewidth', 3);
plot(tspan, [r_SN_NAl; r_SN_NAl], 'b--', 'linewidth', 3);
xlabel('Time (sec)', formatSetting{:});
ylabel('$K_0(x(t))$', formatSetting{:});

grid on;

print('Figure/TwoState_Example','-depsc')

%% function Phase Plot setup
function formatting_PhasePlot(ax1, limsize)
    axis equal;
    grid on;
    % box off;

    ax1.XTick = [-1:0.25:1];
    ax1.YTick = [-1:0.25:1];
    xlim([-limsize, limsize]);
    ylim([-limsize, limsize]);
end