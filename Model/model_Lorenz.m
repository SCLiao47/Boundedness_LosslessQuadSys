

function model = model_Lorenz()
    % dx/dt = c + L*x + [x'*Q1*x, ..., x'*Qn*x]'

    name = "Lorenz_Chaotic";
    
    nx = 3;
    c = zeros(nx,1);
    
    sig = 10;
    rho = 28;
    bet = 8/3;
    L = [-sig sig  0; rho -1 0; 0 0 -bet;];
    
    Q = zeros(nx,nx,nx);
    Q(:,:,1) = 0;
    Q(:,:,2) = [0 0 -1/2; 0 0 0; -1/2 0 0];
    Q(:,:,3) = [0 1/2 0; 1/2 0 0; 0 0 0];
    
    model = class_Model_LosslessQuad(name, c, L, Q);
    
%     nx = 3;
%     model.nx = nx;
%     model.m = zeros(model.nx,1);
%     
%     c = [0; 0; 0];
%     model.c = c;
%     
%     sig = 10;
%     rho = 28;
%     bet = 8/3;
% 
%     L = [-sig sig  0; rho -1 0; 0 0 -bet;];
%     model.L = L;
%     model.Ls = 1/2*(L+L');
% 
%     Q = zeros(nx,nx,nx);
%     Q(:,:,1) = 0;
%     Q(:,:,2) = [0 0 -1/2; 0 0 0; -1/2 0 0];
%     Q(:,:,3) = [0 1/2 0; 1/2 0 0; 0 0 0];
%     model.Q = Q;
%     
%     model.ode = @(t,x) func_xdotNL(c,L,Q,t,x);
end