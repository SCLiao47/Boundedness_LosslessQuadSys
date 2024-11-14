%% Trapping Region Example: System with 2 states

%% System
% dx/dt = c + L x + [x' Q1 x; x' Q2 x]

ifPete = ~true;

if ifPete
    % Quadratic Terms:
    Q1_22 = 1;
    Q1_23 = -4;
    Q1_33 = -3;

    Q2_13 = 2;
    Q2_33 = 0; %Q2_33 = -2;

    Q3_22 = 0; %Q3_22 = 4;

    Q1 = [0 0 0; 0 Q1_22 Q1_23; 0 Q1_23 Q1_33];
    Q2 = [0 -Q1_22/2 Q2_13; -Q1_22/2 0 -Q3_22/2; Q2_13 -Q3_22/2 Q2_33];

    Q3_12 = -(Q1_23+Q2_13);
    Q3 = [0 Q3_12 -Q1_33/2; Q3_12 Q3_22 -Q2_33/2; -Q1_33/2 -Q2_33/2 0];

    % Verify result
    if false
        syms x1 x2 x3 real;
        x = [x1; x2; x3];
        f = [x'*Q1*x; x'*Q2*x; x'*Q3*x];
        simplify(x'*f)
    end

    M = 2*[-Q1_22/2 Q2_13; Q3_12 -Q1_33/2];
    eig(M)

    % Constant term
    c = 0;

    % Linear term
    %L = diag([0.02,-0.2,-0.3]);
    %L = [0.02 0 0; 0.1 -0.2 0; 0 0 -0.3];
    L = [0.2 0 0; 1 -2 0; 0 0 -3];
else
    M = [1, 0; 0, -1];
    
    Q = func_MtoQs_N3(M);
    Q1 = Q(:,:,1);
    Q2 = Q(:,:,2);
    Q3 = Q(:,:,3);
    
    L = [1 0 0; -1 -5 -1; 1 0 -2];

    c = zeros(3,1);
end

Ls = (L+L')/2;

isEN = func_checkEffectiveNon(L, Q);
if isEN
    fprintf('Is the nonlinearity effective? YES \n');
else
    fprintf('Is the nonlinearity effective? NO \n');
end
    
%% Phase plane
% Integrate backwards in time from initial conditions on a small circle
% near the origin.

% System time derivative
fh = @(t,x) c + L*x + [x'*Q1*x; x'*Q2*x; x'*Q3*x];

% Simulation parameters
Nth = 3;
xmax = 100;
x10 = xmax*linspace(-1,1,Nth);
x20 = xmax*linspace(-1,1,Nth);
x30 = xmax*linspace(-1,1,Nth);
Tf = 50;
ax = 2*xmax*[-1 1 -1 1 -1 1];

Param = struct('xbound', ax(1:2), 'ybound', ax(3:4), 'zbound', ax(5:6));
odeOptions = odeset('Event', @(t,x) ode_EventFunc_sim(t,x,Param));

% Generate plots
figure(1); clf
figure(2); clf
% figure(3); clf
figure(4); clf
for i1=1:Nth
    for i2 = 1:Nth
        for i3=1:Nth
            x0 = [x10(i1); x20(i2); x30(i3)];
            [t,x] = ode45(fh,[0 Tf],  x0, odeOptions);

            figure(1); subplot(1, 5, [1:3]);
%             plot3(x(:,1), x(:,2),x(:,3),'b'); hold on;
            plot3(x(:,1), x(:,2),x(:,3)); hold on;

            figure(2)
            subplot(3,1,1)
            plot(t,x(:,1),'b'); hold on;
            subplot(3,1,2)
            plot(t,x(:,2),'b'); hold on;
            subplot(3,1,3)
            plot(t,x(:,3),'b'); hold on;
        end
    end
end

% One more IC
x0 = [20; 0; 0];
[tp,xp] = ode45(fh,[0 Tf],  x0, odeOptions);
x0 = [-20; 0; 0];
[tn,xn] = ode45(fh,[0 Tf],  x0, odeOptions);

figure(1); subplot(1, 5, [1:3]);
plot3(xp(:,1), xp(:,2),xp(:,3),'r'); hold on;
plot3(xn(:,1), xn(:,2),xn(:,3),'m'); hold on;

figure(1); subplot(1, 5, [4:5]);
plot3(xp(:,1), xp(:,2),xp(:,3),'r'); hold on;
plot3(xn(:,1), xn(:,2),xn(:,3),'m'); hold on;

figure(2)
subplot(3,1,1)
plot(tp,xp(:,1),'r'); hold on;
plot(tn,xn(:,1),'m'); hold on;
subplot(3,1,2)
plot(tp,xp(:,2),'r'); hold on;
plot(tn,xn(:,2),'m'); hold on;
subplot(3,1,3)
plot(tp,xp(:,3),'r'); hold on;
plot(tn,xn(:,3),'m'); hold on;


% Add labels
figure(1); subplot(1, 5, [1:3]);
grid on;
xlabel('x1');
ylabel('x2');
zlabel('x3');
% axis(ax)
axis equal
hold off

figure(2)
subplot(3,1,1)
grid on;
ylabel('x1');
hold off

subplot(3,1,2)
ylabel('x2');
grid on;
hold off;

subplot(3,1,3)
ylabel('x3');
xlabel('Time');
grid on;
hold off;


figure(1); subplot(1, 5, [4:5]);
grid on;
xlabel('x1');
ylabel('x2');
zlabel('x3');
axis equal;

figure(4)
plot(tp, vecnorm(xp'), 'r'); hold on;
plot(tn, vecnorm(xn'), 'm');
title('Energy');
xlabel('time');
ylabel('Energy');


%% Solve Primal SDP

cvx_begin sdp quiet
    % Create variables
    variables m1 m2 m3 a;
    dual variable W;
    
    % Create LMI
    W: Ls + m1*Q1 + m2*Q2 + m3*Q3 - a*eye(3) <= 0;
    
    % Objective function:  min trace(P)
    minimize( a );
cvx_end

disp('TRSDP:');
[a, m1, m2, m3]
% astar = a

cvx_begin sdp quiet
    % Create variables
    variables m1 m2 m3 b;
    dual variable Z;
    
    % Create LMI
    Z: Ls + m1*Q1 + m2*Q2 + m3*Q3 - b*eye(3) >= 0;
    
    % Objective function:  min trace(P)
    maximize( b );
cvx_end

disp('ERSDP:');
[b, m1, m2, m3]
% bstar = b

%%

%% Symbolic verifications

% Quadratic term
syms x1 x2 x3 real;
x = [x1; x2; x3];
f = simplify( [x'*Q1*x; x'*Q2*x; x'*Q3*x] );

% Solutions to f(x)=0:
%  1) 0 = x2^2 - 8*x2*x3 - 3*x3^2
%  2) 0 = -x1*(x2 - 4*x3)
%  3) 0 = = x1*(4*x2 + 3*x3)
% Solving these equations give:
%  1) x2 = (4+/- sqrt(19))*x3
%  2) x1=0 and/or x2=4*x3
%  3) x1=0 and/or x2 = -(3/4)*x3
% Thus the solutions of f(x)=0 have the form
%  x = [1; 0; 0]*c for any c
%  x = [0; 4+/-sqrt(19); 1]*c for any c
% Neither of these directions is in an eigenspace of L so f is an
% "effective nonlinearity" as defined by S&N.

% Dynamics
xdot = c + L*x + f;

% Eq. Points: It looks like there are five different eq. points
% and including the origin.
s = solve(xdot);
eqpts = double([s.x1 s.x2 s.x3])'

figure(1); 
subplot(1, 5, [1:3]); hold on;
scatter3(s.x1, s.x2, s.x3, 50, 'fill');
subplot(1, 5, [4:5]); hold on;
scatter3(s.x1, s.x2, s.x3, 50, 'fill');
    