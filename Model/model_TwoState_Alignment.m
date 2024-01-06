function model = model_TwoState_Alignment(ifAlign)
    % an example model to illustrate conservative bound by Schlegel 
    % and Noack. The model is a two-state linear-quadratic system 
    % with negative definte linear part. The system satisfies the TR
    % condition without any coordinate shift. 
    
    % xdot = c + Lx + f(x)
    %   - f(x) is an energy preserving nonlinearity.
    %   - L is choosen as a negative definite symmetric matrix with two
    %     different eigenvalues. L = diag([lam1, lam2]), where lam1 > lam2.
    %   - c is choosen to be ||c|| be a constant with [ifAlign] to the
    %     eigenvector correspond to the largest (least negative)
    %     eigenvalues.

    if nargin < 1
        ifAlign = false;
    end
    
    if ifAlign
        name = "TwoStateBD_Aligned";
        c = [1; 0];
    else
        name = "TwoStateBD_Aligned";
        c = [0; 1];
    end
    
    L = diag([-1, -4]);
    
    Q1 = [0 -1/2; -1/2 0];
    Q2 = [1 0; 0 0];
    Q = cat(3, Q1, Q2);
    
    model = class_Model_LosslessQuad(name, c, L, Q);
end