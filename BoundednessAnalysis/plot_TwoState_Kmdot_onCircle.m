function [ph] = plot_TwoState_Kmdot_onCircle(model, V, scale, center)
% Compute Kmdot on a circel at "center" with radius "scale*V", where Kmdot
% is the time derivative of kinetic energy of model define at "center".
%
% The energy Km is compute by shifting the system to "center"
% Plot is done on the original coordinate

%% input
if nargin < 4
    center = [0;0];
end

%% options
num_th = 50;

%% computing Kmdot on circle
th_onC = linspace(-pi, pi, num_th);
unitDisk = [cos(th_onC); sin(th_onC)];

ring_shifted = unitDisk*scale;

Kdot_onC = 0*th_onC;
xdot_onC = 0*unitDisk;

% shift model to center 
model_shifted = func_ShiftSystem(model, center);

% use new shifted for energy and ode
for i = 1:num_th
    Kdot_onC(i) = model_shifted.dKm(ring_shifted(:,i));
    
    xdot_onC(:,i) = model_shifted.ode_shifted(0, ring_shifted(:,i));
end

%% plotting circle on the original coordinate
qX = center(1) + ring_shifted(1,:);
qY = center(2) + ring_shifted(2,:);
qU = Kdot_onC.*unitDisk(1,:);
qV = Kdot_onC.*unitDisk(2,:);

% plot center
ph.center = scatter(center(1), center(2), 100, 'o', 'filled');

% patch([ring_shifted(1,:), nan], [ring_shifted(2,:), nan], [Kdot_onC, nan],...
%       'EdgeColor','interp', 'LineWidth', 2)
ph.ring = patch([qX, nan], [qY, nan], [Kdot_onC, nan],...
      'EdgeColor','interp', 'LineWidth', 2);

colormap jet;
colorbar;

% cmap = colormap;
% cK = @(K) floor(interp1([min(Kdot_onC), max(Kdot_onC)], [1, size(cmap,1)], K));
% 
% % normal direction quiver on the circle
% for i = 1:num_th
%     quiver(qX(i), qY(i), qU(i), qV(i), .02, ...
%            'Color', cmap(cK(Kdot_onC(i)),:), 'linewidth', 1.5);
% end

% dynamics on the circle
ph.dynamics = quiver(qX, qY, xdot_onC(1,:), xdot_onC(2,:));
end