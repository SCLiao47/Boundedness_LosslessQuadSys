function model = model_TwoState_CaseStudy(CaseNum, mag_La)   
% a function to create 2D lossless model to study influence of the 
% boundedness by Ls, and the results from {TRSDP and ERSDP}. The 2D case
% has the same two feasible points for both TRSDP and ERSDP. The two points 
% are different only by a negative sign, resulting the same objective value
% in both SDP. 
%   
% [Model] 
%   xdot = 0 + (Ls+La)x + f(x)
%       - zero constant dynamics
%       - Ls is symmetric part of linear dynamics
%       - La is asymmetric part of linear dynamics
%       - f(x) is energy-preserving nonlinearties
%
% [Factor]
%   - symmetric part of linear matrix
%   - quadratic function, influencing the feasible solution of SDPs, a*
%
% [Cases]
% 
%   - Case 1: Ls is PD
%   - Case 2: Ls is ND
%   - Case 3: Ls is indefinite, a*>0
%   - Case 4: Ls is indefinite, a*<0
%   - Case 5: Ls is indefinite, a*=0
%   - Case 6: Ls is PSD, a*>0
%   - Case 7: Ls is PSD, a*=0
%   - Case 8: Ls is NSD, a*<0
%   - Case 9: Ls is NSD, a*=0
%
% To discuss more:
%   - how asymmetric part affect the system?

    switch nargin
        case 0
            CaseNum = 1;
            mag_La = 0;
            
        case 1
            mag_La = 0;
    end
    
    % bound the CaseNum argument
    CaseNum = max(min(CaseNum, 9), 1);
    
    % dynamics
    c = [0; 0];
    [L, Q] = get_dynamicsCoef(CaseNum, mag_La);
    

    name = "TwoState_CaseStudy_Case"+num2str(CaseNum);
    model = class_Model_LosslessQuad(name, c, L, Q);
end

function [L, Q] = get_dynamicsCoef(CaseNum, mag_La)
    % set asymmetric part
    La = [0, 1; -1, 0];
    
    % set Ls and Q
    switch CaseNum
        case 1 % Ls is PD
            Ls = diag([1,1]);
            qa = 1;
            qb = 1;
            
        case 2 % Ls is ND
            Ls = diag([-1,-1]);
            qa = 1;
            qb = 1;
            
        case 3 % Ls is indefinite, a*>0
            Ls = diag([1, -1]);
            qa = 0;
            qb = 1;
            
        case 4 % Ls is indefinite, a*<0
            Ls = diag([1, -1]);
            qa = 1;
            qb = 0;
            
        case 5 % Ls is indefinite, a*=0
            Ls = diag([1, -1]);
            qa = 1;
            qb = 1;
            
        case 6 % Ls is PSD, a*>0
            Ls = diag([1, 0]);
            qa = 0;
            qb = 1;

        case 7 % Ls is PSD, a*=0
            Ls = diag([1, 0]);
            qa = 1;
            qb = 0;
            
        case 8 % Ls is NSD, a*<0
            Ls = diag([-1, 0]);
            qa = 0;
            qb = 1;
            
        case 9 % Ls is NSD, a*=0
            Ls = diag([-1, 0]);
            qa = 1;
            qb = 0;
            
        otherwise
            error(['CaseNum #', num2str(CaseNum), ' not implemented']);
    end
    
    L = mag_La*La + Ls;
    
    Q1 = [0 qa; qa 2*qb];
    Q2 = [-2*qa -qb; -qb 0];
    Q = cat(3,Q1,Q2); 
end