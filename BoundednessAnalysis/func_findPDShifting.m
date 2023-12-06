
function [m, info] = func_findPDShifting(model, option)

    disp('==== Inverse-time study: START ====');
    [m, info] = func_findPDShifting_invTime(model, option);
    
    % [for debugging]
    [m_direct, info_direct] = func_findPDShifting_direct(model, option);
    
    assert( sum(abs(m - m_direct)) == 0, ...
        "INVTIME and DIRECT methods are not consistent!")
    
    disp('==== Inverse-time study: END ====');
end

%implement finding PD shift m with findNDShifting with inverstime model
function [m, info] = func_findPDShifting_invTime(model, option)
    model_invT = get_inverseTimeModel(model);

    [mND, info_ND] = func_findNDShifting(model_invT, option);
    
    % convert to forward-time solution
    m = mND;
    
    info.feasibility = info_ND.feasibility;
    
    if info.feasibility
        % primal solution
        info.b = - info_ND.a; 
        info.As = - info_ND.As;

        % dual solution
        info.Z = info_ND.W;
        
        % ER
        info.existER = info_ND.existTR;
    end
end

% implement finding PD shift m with directly using SDP
function [m, info] = func_findPDShifting_direct(model, option)
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
    % This SDP push As most further to positive definite (the most negative
    % eigenvalue is maximize). 
    
    cvx_begin sdp quiet
    cvx_solver mosek
%     cvx_solver SeDuMi
        variable b(1,1);
        variable m(nx,1);
        dual variable Z;

        maximize(b);

        subject to 

            As = Ls;
            for i = 1:nx
                As = As - m(i)*Q(:,:,i);
            end
            Z : As >= Inx*b;
    cvx_end
    
    % Regularize the coordinate shifting 
    if strcmp(cvx_status, 'Unbounded')
        warning('Unregularized SDP is unbounded (b*=inf). Rerun SDP with regularization');
        cvx_begin sdp quiet
        cvx_solver mosek
    %     cvx_solver SeDuMi
            variable m_bounded(nx,1);
            dual varabile Z;

            minimize(norm(m_bounded));

            subject to 
                As = Ls;
                for i = 1:nx
                    As = As - m_bounded(i)*Q(:,:,i);
                end
                Z : As >= 1e-6*Inx;
        cvx_end
        
        m = m_bounded;
    end

    if strcmp(cvx_status, 'Solved')       
        info.feasibility = true;
        info.b = b;
        info.As = As;
        info.Z = Z;
        
        if b > 0
            info.existER = true;
            
            if option.verbose
                disp('Coordinated shift s.t. As>0 is found:');
                disp(m')
                fprintf('Most negative eigenvalue b = %.4f\n', b)
            end
        else
            info.existER = false;
        end
    else
        if option.verbose
            disp(['cvx_status: ', cvx_status]);
%             disp('could not find m such that As<0');
        end
    end
end