function [m, info] = func_findNDShifting(model,option)
% solve trapping region condition using SDP. 
%
% Corresponding to Section 3.1 of Liao et. al, 2024

    if nargin < 2
        option.round_Ndigit = 3;
        option.verbose = false;
        option.tol = 1e-6;
    end
    
    % initialize output
    m = nan;
    info.feasibility = false;
    
    %% system param
    nx = model.nx;
    Ls = model.Ls;
    Q = model.Q;
    
    Inx = eye(nx);
    
    %% filtering the data matrix
    % if don't filter, set option.round_Ndigit = nan.
    if ~isnan(option.round_Ndigit)
        Ls = round(Ls, option.round_Ndigit);
        Q = round(Q, option.round_Ndigit);
    end
    
    %% solve by SDP
    % This SDP push As most further to negative definite (the most positive
    % eigenvalue is minimized). 
    cvx_begin sdp quiet
    cvx_solver mosek
%     cvx_solver SeDuMi
        variable a(1,1);
        variable m(nx,1);
        dual variable W;

        minimize(a);

        subject to 
            As = Ls;
            for i = 1:nx
                As = As - m(i)*Q(:,:,i);
            end
            W : As <= Inx*a;
    cvx_end
    
    % The SDP can be unbounded, i.e., a*=-inf. 
    % Redo the SDP with regularization of the coordinate shifting 
    if strcmp(cvx_status, 'Unbounded')
        warning('Unregularized TRSDP is unbounded below(a*=-inf). Rerun TRSDP with regularization');
        cvx_begin sdp quiet
        cvx_solver mosek
    %     cvx_solver SeDuMi
            variable m_bounded(nx,1);
            dual variable W;

            minimize(norm(m_bounded));

            subject to 
                As = Ls;
                for i = 1:nx
                    As = As - m_bounded(i)*Q(:,:,i);
                end
                W : As <= -1e-3*Inx;
        cvx_end
        
        m = m_bounded;
        a = max(eig(As));
    end
    
    % setting output
    if strcmp(cvx_status, 'Solved')        
        info.feasibility = true;
        info.a = a;
        info.As = As;
        info.W = W;
        
        if a < 0
            info.existTR = true;
            if option.verbose
                disp('Coordinated shift s.t. As<0 is found:');
                disp(m')
                fprintf('Most positive eigenvalue a = %.4f\n', a)
            end
        else
            info.existTR = false;
        end
    else
        if option.verbose
            disp(['cvx_status: ', cvx_status]);
            disp('could not find m such that As<0');
        end
    end
    
    info.m = m;
end