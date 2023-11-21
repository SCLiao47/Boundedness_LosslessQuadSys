

function [m, info] = func_findNDShifting(model,option)
    if nargin < 2
        option.verbose = false;
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
    if ~isnan(option.round_Ndigit)
        Ls = round(Ls, option.round_Ndigit);
        Q = round(Q, option.round_Ndigit);
    end
    
    %% solve by SDP
    % This SDP push As most further to negative definite (the most positive
    % eigenvalue is minimized). 
    % [!!!] We conjecture that this eigenvalue is lower bounded.
    
    cvx_begin sdp quiet
    cvx_solver mosek
%     cvx_solver SeDuMi
        variable a(1,1);
        variable m(nx,1);
        dual variable W;

        minimize(a);

        subject to 
%             a <= 0;

            As = Ls;
            for i = 1:nx
                As = As - m(i)*Q(:,:,i);
            end
    %         As <= -1e-6*Inx;
            W : As <= Inx*a;
    cvx_end
    
    % Regularize the coordinate shifting 
    if strcmp(cvx_status, 'Unbounded')
        warning('Unregularized SDP is unbounded (a*=-inf). Rerun SDP with regularization');
        cvx_begin sdp quiet
        cvx_solver mosek
    %     cvx_solver SeDuMi
%             variable a;
            variable m_bounded(nx,1);
            dual varabile W;

            minimize(norm(m_bounded));

            subject to 
                As = Ls;
                for i = 1:nx
                    As = As - m_bounded(i)*Q(:,:,i);
                end
                W : As <= -1e-6*Inx;
%                 As <= a*Inx;
        cvx_end
        
        m = m_bounded;
    end
    
    % setting output
    if strcmp(cvx_status, 'Solved')
%         if option.verbose
%             disp('Coordinated shift s.t. As<0 is found:');
%             disp(m')
%             fprintf('Most positive eigenvalue a = %.4f\n', a)
%         end
        
        info.feasibility = true;
        info.a = a;
        info.As = As;
        info.W = W;
        
        if a < 0
            info.existTR = true;
        else
            info.existTR = false;
        end
    else
        if option.verbose
            disp(['cvx_status: ', cvx_status]);
%             disp('could not find m such that As<0');
        end
    end
    
    info.m = m;
end