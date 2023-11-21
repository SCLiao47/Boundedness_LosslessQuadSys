

function model = model_TwoState_Unbounded(caseNum)
% example 2D models to illustrate unbounded behaviors

if nargin < 1
    caseNum = 1;
end


%% model info
model.name = "Case"+num2str(caseNum); 
model.nx = 2;
model.m = zeros(model.nx,1);

%% dynamics
% dx/dt = c + L*x + [x'*Q1*x, ..., x'*Qn*x]'
[c, L, Q] = get_dynamicsCoef(caseNum);

model.c = c;    

% L = [-10 -5; -5 -1];

model.L = L;
model.Ls = 1/2*(L+L');

model.Q = Q;

% ode
model.ode = @(t,x) ode_quadraticDyn(model.c, model.L, model.Q, t, x);

%% Energy
model.dK0 = @(x) model.c'*x + x'*model.Ls*x;

end

function [c, L, Q] = get_dynamicsCoef(caseNum)
    c = [0;0];

    switch caseNum
        case 1
            % unbounded: 1 unstable EQ
            L = diag([1,1]);

            Q121 = 1;
            Q221 = 1;
            
        case 2
            % unbounded: 1 stabel, 2 unstable
                
%             c = [10;0];
            
            s = -10;
%             a = 5;
            a = 0;
            
            Ls = [-10 s; s -1];            
            La = -a*[0 1; -1 0];
            
%             Ls = [11 5; 5 3.5]; % to have the same eigenvector of null space
            
            L = Ls + La;
            
            Q121 = 1;
            Q221 = 1;
            
        otherwise % case 1
            disp('Case not specified, set to Case1');
            L = diag([1,1]);

            Q121 = 1;
            Q221 = 1;
    end
    
    Q111 = 0;
    Q112 = -2*Q121;
    Q122 = -Q221/2;
    Q222 = 0;
    
    Q1 = [Q111 Q121; Q121 Q221];
    Q2 = [Q112 Q122; Q122 Q222];
    Q = cat(3,Q1,Q2); 
end
