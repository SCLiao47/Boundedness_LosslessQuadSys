% This script is to verify the implementation of the rotation of quadratic
% systems. The script use Lorenz system as a test. The coordinate "x" is in
% the original cooridnate of Lorenz system, while coordinate "z" is the
% rotated coordintate with z=R*x.

init;

%% L:   Lorenz in original coordinate
L = model_KStackedLorenz(1, false);

x0 = [0, 100, 100]';
tSpan = [0, 20];

[tL, xL] = ode45(L.ode, tSpan, x0);

%% M:   Lorenz in rotated coordinate
% rng(0);

[M, R] = model_KStackedLorenz(1, true);

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
set(gcf, 'Position', [744 400.2000 1.1714e+03 649.8000])

fontsetting = {'fontsize', 16};

subplot(1,2,1);
plot3(zM(:,1), zM(:,2), zM(:,3), 'r');

title('Rotated coordinate', fontsetting{:});
xlabel('z1', fontsetting{:});
ylabel('z2', fontsetting{:});
zlabel('z3', fontsetting{:});
plot_setting()

%
subplot(1,2,2);
plot3(xL(:,1), xL(:,2), xL(:,3), 'b');

hold on;
plot3(xM(:,1), xM(:,2), xM(:,3), 'r--');

legend(["x(t)", "R'.*z(t)"], 'Location', 'northwest', fontsetting{:});
title('Original coordinate', fontsetting{:});
xlabel('x1', fontsetting{:});
ylabel('x2', fontsetting{:});
zlabel('x3', fontsetting{:});
plot_setting()

%% Verifying on the boundedness analysis
% 1. The rotated coordinate shift can be different from the one in the
%    original coordinate. As the shift along x1=0 does not change the
%    symmetric part of the shifted dynamics As, only change the asymmetric
%    part Aa. This implies sum(R'*mM - mL) might not equals to 0.
% 2. As R'*mM and mL are different, the resulting TR regions would be
%    different as well.  
% 3. The things would be the same is the eigenvalues of As for both models.
%    Also, the largest (most positive) eignvalues a*.
%
% [Note]
% 1. SDP analysis of findingNDShifting containts rounding to the system
%    data. This needs to be turn off for testing rotation, as the rotation
%    makes the data sensitive to the small digits. 

% options
option.round_Ndigit = nan;
option.verbose = false;
option.tol = 1e-6;

% TR 
[mM, info_M] = func_findNDShifting(M, option);
[mL, info_L] = func_findNDShifting(L, option);

assert(abs(info_M.a - info_L.a) < option.tol, ...
       "a^* is not the same for models with/without rotation");

% test the rotated shifting 
mM_x = R'*mM;

M_shift = func_ShiftSystem(M, info_M.m);
L_shift = func_ShiftSystem(L, info_L.m);

eig_AsM = sort(eig(M_shift.As));
eig_AsL = sort(eig(L_shift.As));

assert(all(abs(eig_AsM - eig_AsL) < option.tol), ...
       "Evalues does not match for models with/without rotation"); 

%%
[L2, ~] =  model_KStackedLorenz(2, false);
[M2, R2] = model_KStackedLorenz(2, true);

x02 = [x0; R*x0];
z02 = R2*x02;

[tL2, xL2] = ode45(L2.ode, tSpan, x02);
[tM2, zM2] = ode45(M2.ode, tSpan, z02);

% translate zM back to xM
xM2 = 0*zM2;
for i = 1:length(tM2)
    xM2(i,:) = (R2'*zM2(i,:)')';
end

figure(2); clf;
set(gcf, 'Position', [700 100 1120 350])

subplot(1,2,1); 
plot3(xL2(:,1), xL2(:,2), xL2(:,3), 'b');  hold on;
plot3(xM2(:,1), xM2(:,2), xM2(:,3), 'r--');
xlabel('x1', fontsetting{:});
ylabel('x2', fontsetting{:});
zlabel('x3', fontsetting{:});
plot_setting()

subplot(1,2,2);
plot3(xL2(:,4), xL2(:,5), xL2(:,6), 'b');  hold on;
plot3(xM2(:,4), xM2(:,5), xM2(:,6), 'r--');
xlabel('x4', fontsetting{:});
ylabel('x5', fontsetting{:});
zlabel('x6', fontsetting{:});
plot_setting()

sgtitle('2-Stacked Lorenz System')
legend(["x(t)", "R'.*z(t)"], 'Location', 'northeast', fontsetting{:})

%% function 
function plot_setting()
    set(gca, 'View', [-14.7210 19.4215]);
    set(findobj(gcf,'type','line'), 'linewidth', 2)
    grid on
end