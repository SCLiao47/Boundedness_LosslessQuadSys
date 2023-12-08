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
Plot_LorTraj = figure(1); clf
Plot_LorTraj.WindowStyle='normal';
Plot_LorTraj.DockControls='off';
width=900;
height=900;
set(gcf,'units','pixels','position',[1 1 width height])

% figure(1); clf
% subplot(4,1,[1,2,3]);
tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
ax1=nexttile([3,3]);

% plot coordinate
xlimits = [-150, 150];
ylimits = [-85, 150];
plot(xlimits+[-10, 10], [0,0], 'k', ...
     'linewidth', 2, 'HandleVisibility','off');
hold on;
plot([0,0], ylimits, 'k', ...
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
    
    % yz-plan
    plot(x(:,2),x(:,3), 'color', [0.5, 0.5, 0.5]);
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
ph_y = scatter(ystar(2,:), ystar(3,:), 150, 'filled', 'Diamond');
set(ph_y, 'MarkerFaceColor', [0.4940 0.1840 0.5560]);


% ===[ formatting figure ]===
% % [ 3D ]
% xlabel('x1');
% ylabel('x2');
% zlabel('x3');

% [ 2D ]
fontsize = 20;
fontsize2 = 15;
formatSetting = {'interpreter', 'latex', 'fontsize', 20};

text(xlimits(2)-10, 6, '$x_2$', formatSetting{:})
text(10, ylimits(2)-5, '$x_3$', formatSetting{:})
% xlabel('$x_2$','Interpreter','latex','fontsize', fontsize);
% ylabel('$x_3$','Interpreter','latex','fontsize', fontsize);
set(gca,'TickLabelInterpreter','latex','FontSize', fontsize2);

xlim(xlimits);
ylim(ylimits);
axis equal
grid on;

%% visualization 2, Energy plot
% subplot(4,1,4);
ax2=nexttile([1,2]);

% plot
for i = 1:10
    t = Trajs{i}.t;
    Em = Trajs{i}.Em;
    
    plot(t, Em, 'color', [0.5,0.5,0.5]);
    hold on;
end

plot(tspan, [r; r], 'r--', 'linewidth', 3);
plot(tspan, [r_SN; r_SN], 'b--', 'linewidth', 3);
xlabel('Time (sec)','fontsize', 16);
ylabel('$K_m(x(t))$','Interpreter','latex','fontsize', 16);

% set(gcf, 'Position', [462.6000 13 831.2000 749]);

lgd=legend([ph_E0, ph_S2, ph_S1, ph_y], ...
       {'$E$', '$B(m,R_m)$', '$B(m,R_m^*)$', '$y^*$'}, ...
       'fontsize', fontsize2, 'interpreter', 'latex');
lgd.Layout.Tile=8;

%%
print('Figure/Lorenze_TR_2D','-depsc')