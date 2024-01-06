% This script demonstrate verify the boundedness of Lorenze system

init;

%% setups

model = model_Lorenz();

% options
option.round_Ndigit = 3;
option.verbose = true;
option.tol = 1e-6;

%% verify the existence of boundedness region
[m, info_m] = func_findNDShifting(model, option);
model_shifted = func_ShiftSystem(model, m);

%% solve for the size of boundedness region
[r, info_TR] = func_TRSize_SDP(model_shifted, option);

r_SN = norm(model_shifted.d) / abs(eigs(model_shifted.As, 1, 'sm'));

%% Computing the ellipsoid Edot = 0
% In this problem, As is already a diagonal matrix. Hence no need to
% transform coordinate
lamb = diag(model_shifted.As);
d = model_shifted.d;

temp = sum(d.^2./lamb);
alpha = 1/2 * sqrt(temp./lamb);
cE0 = 1/2 * d./lamb;

%% trajectories
if_sim = ~true;

if if_sim
    num_traj = 100;
    radius = r_SN*2;
    tspan = [0, 5];

    x0s = (rand(model.nx, num_traj)-0.5)*radius;
    Trajs = cell(num_traj, 1);

    for i = 1:num_traj
        [t, traj] = ode45(model.ode, tspan, x0s(:,i));

        Em = vecnorm((traj'-m))';
        Trajs{i} = struct('t', t, 'x', traj, 'Em', Em);
    end
    
    save("Data/Lorenze_traj", "Trajs", "num_traj", "tspan");
else
    load("Data/Lorenze_traj");
end

%% visualization 1, State-space plot

Tmax = 2;

%
Plot_LorTraj = figure(1); clf
Plot_LorTraj.WindowStyle='normal';
Plot_LorTraj.DockControls='off';

% width = 600; height = 800;
% width = 400; height = width*4/3;
width = 700; height = width/7*4;
% width = 700; height = 385.6;
set(gcf,'units','pixels','position',[1 1 width height])

% tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
% ax1=nexttile([3,1]);
tiledlayout(1,7, 'TileSpacing','compact','Padding','compact');
ax1 = nexttile([1,4]);

% plot coordinate
xlimits = [-120, 120];
ylimits = [-80, 150];
plot(xlimits+[-10, 10], [0,0], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
hold on;
plot([0,0], ylimits+[-10, 10], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');

% plot trajectory 
for i = 1:10
    t = Trajs{i}.t;
    x = Trajs{i}.x;
    x0 = x(1,:);
    
%     % 3D
%     plot3(x(:,1),x(:,2),x(:,3), 'color', [0.5, 0.5, 0.5]);
%     hold on;
%     scatter3(x0(1), x0(2), x0(3), 20, 'filled','MarkerFaceColor', [0.5,0.5,0.5]);
    
    idx = t <= Tmax;

    % yz-plan
    plot(x(idx,2),x(idx,3), 'color', [0.5, 0.5, 0.5]);
    hold on;
    scatter(x0(2), x0(3), 20, 'filled','MarkerFaceColor', [0.5,0.5,0.5]);
end

% % ===[ 3D ]===
% % Trapping region
% [X,Y,Z] = sphere;
% 
% ph_S1 = surf(X*r + m(1), Y*r + m(2), Z*r + m(3));
% set(ph_S1, 'FaceColor', 'None', 'EdgeColor', 'r');
% 
% ph_S2 = surf(X*r_SN + m(1), Y*r_SN + m(2), Z*r_SN + m(3));
% set(ph_S2, 'FaceColor', 'None', 'EdgeColor', 'b');
% 
% ph_E0 = surf(X*alpha(1) - cE0(1) + m(1), ...
%              Y*alpha(2) - cE0(2) + m(2), ...
%              Z*alpha(3) - cE0(3) + m(3));
% set(ph_E0, 'FaceColor', 'None', 'EdgeColor', 'g');
% 
% % ciritcal points
% ystar = info_TR.ystar + m;
% ph_y = scatter3(ystar(1,:), ystar(2,:), ystar(3,:), 100, 'filled', '^');

% ===[ 2D ]===
% Trapping region
th = linspace(0,2*pi, 200);
Y = cos(th);
Z = sin(th);

ph_S1 = plot(Y*r + m(2), Z*r + m(3));
set(ph_S1, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3);

ph_S2 = plot(Y*r_SN + m(2), Z*r_SN + m(3));
set(ph_S2, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 3);

% ph_E0 = plot(Y*alpha(2) - cE0(2) + m(2), Z*alpha(3) - cE0(3) + m(3));
% set(ph_E0, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 3);

ph_E0 = patch(Y*alpha(2) - cE0(2) + m(2), Z*alpha(3) - cE0(3) + m(3), 'g');
set(ph_E0, 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'LineWidth', 2);
% set(ph_E0, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 3);

% ciritcal points
ystar = info_TR.ystar + m;
ph_y = scatter(ystar(2,:), ystar(3,:), 100, 'filled', 'Diamond');
set(ph_y, 'MarkerFaceColor', [0.4940 0.1840 0.5560]);


% ===[ formatting figure ]===
% % [ 3D ]
% xlabel('x1');
% ylabel('x2');
% zlabel('x3');

% [ 2D ]
conf = config();
size_label = conf.plotting.size_label;
size_legend = conf.plotting.size_legend;
size_tick = conf.plotting.size_tick;
formatSetting = {'interpreter', 'latex', 'fontsize', size_label};

set(gca,'TickLabelInterpreter','latex','FontSize', size_tick);
ax1.XTick = [-100:25:100];
ax1.YTick = [-75:25:150];

% text(xlimits(2)-16, -8, '$x_2$', formatSetting{:})
% text(-17, ylimits(2)-1, '$x_3$', formatSetting{:})
xlabel('$x_2$',formatSetting{:});
ylabel('$x_3$',formatSetting{:});


xlim(xlimits);
ylim(ylimits);
axis equal
grid on;
box on;

lgd=legend([ph_E0, ph_S2, ph_S1, ph_y], ...
       {'$E$', '$B(m,R_m)$', '$B(m,R_m^*)$', '$x^*$'}, ...
       'fontsize', size_legend, 'interpreter', 'latex');
lgd.Location = 'NorthEast';

%% visualization 2, Energy plot
% ax2=nexttile([1,1]);
ax2 = nexttile([1,3]);

% plot
for i = 1:10
    t = Trajs{i}.t;
    Em = Trajs{i}.Em;
    
    idx = t<=Tmax;
    
    plot(t(idx), Em(idx), 'color', [0.5,0.5,0.5]);
    hold on;
end

plot([0, Tmax], [r; r], 'r--', 'linewidth', 2);
plot([0, Tmax], [r_SN; r_SN], 'b--', 'linewidth', 2);

set(gca,'TickLabelInterpreter','latex','FontSize', size_tick);
xlabel('Time (sec)', formatSetting{:});
ylabel('$K_m(x(t))$', formatSetting{:});

grid on;
box on;

width = 700; height = 379;
set(gcf,'units','pixels','position',[1 1 width height])

%%

print('Figure/Lorenze_TR_2D','-depsc')
exportgraphics(Plot_LorTraj, 'Figure/Lorenze_TR_2D.pdf')

%% for github repo
% width = 700; height = 500;
set(gcf,'units','pixels','position',[1 1 width height])

lgd=legend([ph_E0, ph_S2, ph_S1, ph_y], ...
       {'States with Non-Decreasing Energy', 'Trapping Region: Prior State-of-the-Art', 'Trapping Region: Proposed Method', 'Critical States'}, ...
       'fontsize', size_legend);

exportgraphics(Plot_LorTraj, 'Figure/Lorenze_TR_2D.png')