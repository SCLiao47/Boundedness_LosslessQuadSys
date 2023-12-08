% This script is to verify the implementation of the rotation of quadratic
% systems. The script use Lorenz system as a test. The coordinate "x" is in
% the original cooridnate of Lorenz system, while coordinate "z" is the
% rotated coordintate with z=R*x.

init;

%% Lorenz in original coordinate
L = model_Lorenz();

x0 = [0, 100, 100]';
tSpan = [0, 20];

[tL, xL] = ode45(L.ode, tSpan, x0);

%% Lorenz in rotated coordinate
rng(0);

[M, R] = model_KStackedLorenz(1,true);

z0 = R*x0;
[tM, zM] = ode45(M.ode, tSpan, z0);

% translate zM back to xM
xM = 0*zM;
for i = 1:length(tM)
    xM(i,:) = (R'*zM(i,:)')';
end

disp('Rotation:');
disp(R)

%% visualization
figure(1); clf

plot3(xL(:,1), xL(:,2), xL(:,3), 'b');

hold on;
plot3(xM(:,1), xM(:,2), xM(:,3), 'r--');

set(gca, 'View', [-14.7210 19.4215]);
set(findobj(gcf,'type','line'), 'linewidth', 2)
grid on

fontsetting = {'fontsize', 16};
legend(["x(t)", "R'.*z(t)"], 'Location', 'northwest', fontsetting{:});
xlabel('x1', fontsetting{:});
ylabel('x2', fontsetting{:});
zlabel('x3', fontsetting{:});

