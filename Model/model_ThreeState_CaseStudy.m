function model = model_ThreeState_CaseStudy(CaseNum, mag_La)
% example 3D models to illustrate unbounded behaviors

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
    c = [0; 0; 0];
    [L, Q] = get_dynamicsCoef(CaseNum, mag_La);
    
    name = "ThreeState_CaseStudy_Case"+num2str(CaseNum);
    model = class_Model_LosslessQuad(name, c, L, Q);
end

function [L, Q] = get_dynamicsCoef(CaseNum, mag_La)
    % set asymmetric part
    % La = [0, 1; -1, 0];
    La = zeros(3);
    
    % set Ls and Q
    switch CaseNum
        case 1 % Pete's example on Jan. 21
            Ls = [0.2,  0,   0;
                  1,   -2,   0;
                  0,    0,  -3];
             
            Q1_22 = 1;
            Q1_23 = -4;
            Q1_33 = -3;
            
            Q2_13 = 2;
            Q2_33 = 0;
            
            Q3_22 = 0;
            
            Q1 = [0 0 0; 0 Q1_22 Q1_23; 0 Q1_23 Q1_33];
            Q2 = [0 -Q1_22/2 Q2_13; -Q1_22/2 0 -Q3_22/2; Q2_13 -Q3_22/2 Q2_33];

            Q3_12 = -(Q1_23+Q2_13);
            Q3 = [0 Q3_12 -Q1_33/2; Q3_12 Q3_22 -Q2_33/2; -Q1_33/2 -Q2_33/2 0];
            
        otherwise
            error(['CaseNum #', num2str(CaseNum), ' not implemented']);
    end
    
    L = mag_La*La + Ls;
    
    Q = zeros(3,3,3);
    Q(:,:,1) = Q1;
    Q(:,:,2) = Q2;
    Q(:,:,3) = Q3;
end