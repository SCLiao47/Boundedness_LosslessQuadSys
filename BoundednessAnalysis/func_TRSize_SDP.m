function [rTrap, info]= func_TRSize_SDP(model, option)
% solving the size of trapping region by QCQP through dual SDP
%
% % Corresponding to Section 3.2 and 3.3 of Liao et. al, 2024

    if nargin < 2
        option.verbose = false;
        option.tol = 1e-6;
    end
    
    % extract parameter
    As = model.As;
    d = model.d;

    % check if As is negative definite
    [~, flag] = chol(-As);
    if ~(flag == 0)
        disp('As is not negative definite!');
        
        rTrap = inf;
        info = [];
        
        return 
    end
    
    % utility matrix
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    Znx1 = zeros(nx,1);

    %% SDP to find max y'y s.t. y'Asy+d'y>=0
    flag_Lag = false;
    % flag_Rank = false;
    
    %% [Lagrangian dual SDP formulation]
    cvx_begin sdp quiet
    cvx_solver mosek       
        variable gam(1,1);
        variable lam(1,1)       nonnegative;
        
        minimize(gam);
        
        subject to 
            [gam Znx1'; Znx1 -Inx] + lam*[0 d'/2; d/2 -As] >= 0;
    cvx_end
    cvx_status_Lag = cvx_status;
    
    if strcmp(cvx_status_Lag, 'Solved')  
        flag_Lag = true;
        
        % check Astar <= 0
        Astar = Inx + lam*As;
        % check d \in range(Astar)
        assert(norm(Astar*d) >= eps, "d is not in range of Astar!"); 
        
        %% [compute ystar using KKT condition] 
        % ystar = y0 + v*c, where y0 = -lam/2*(Astar\d) and v is the
        %   nullspace of Astar
        %
        % The solution is characterized by a (n-r)-dimensional sphere,
        % where r=rank(Astar). Note that r>0 as d is in range(Astar).
        %
        % By complementary slackness (CS), y'*As*y+d'*y == 0 since lam>0. 
        %  Solve the coefficient c by CS.
        y0 = -lam/2 * (Astar\d);
        
        [~,S,V] = svd(Astar);
        r = sum(abs(diag(S)) > option.tol);
        v = V(:, (r+1):nx);
        
        if r == nx
            ystar = y0;
        elseif r == nx-1           
            % solve CS, which is a quadratic equation in c
            coef1 = v'*As*v;
            coef2 = 2*v'*As*y0 + d'*v;
            coef3 = y0'*As*y0 + d'*y0;
            
            s = sqrt(coef2^2 - 4*coef1*coef3);
            c = (-coef2 + [s, -s]) / (2*coef1);
            
            ystar = y0 + v*c;
        else
            % rank(Astar) <= nx - 2, the solution is characterized by a
            % high-dimensional shpere. Do not solve. 
            ystar = ['A ', num2str(nx-r), '-dimensional shpere.'];
        end
        
        % [check ystar solutions]
        if r >= nx-1
            Lag = @(y, lam) y'*y + lam*(y'*As*y + d'*y); 

            for i = 1:size(ystar,2)
                ys = ystar(:,i);

                % 3. gam* = L(ystar, lam*)
                assert(check_RelTol(gam, Lag(y0, lam), option), ... 
                        "Error: gam* == L(ystar, lam*)");
                assert(check_RelTol(gam, Lag(ys, lam), option), ...
                        "Error: gam* == L(ystar, lam*)");

                % 4. ystar'*ystar == L(ystar, lam*)
                [ifPass, RelErr] = check_RelTol(gam, ys'*ys, option);
                if ~ifPass      % if ~check_RelTol(gam, ys'*ys, option)
                    warning(strcat("RelTol not satisified: gam* == ystar'*ystar with relative error ", num2str(RelErr)));
                end
            end
        end
    end
    
    %% Set output
    % initialize output
    rTrap = nan;
    info = struct('feasibility',false,'cvx',nan);
    
%     if all([flag_Lag, flag_Rank])
    if flag_Lag
        rTrap = sqrt(gam);
        
        info.cvx_Lag = cvx_status_Lag;
        info.feasibility = true;
        
        % dual solutions
        info.gam = gam;
        info.lam = lam;
        % primal solutions
        info.ystar = ystar;
        info.rank = r;
        
        if option.verbose
            disp('Trapping region found!');
            fprintf('TR size = %.3f \n', rTrap);
            fprintf('y* =\n')
            disp(ystar');
        end
    end
end

function [ifPass, RelErr] = check_RelTol(a, b, option)
    RelErr = abs(a-b) / max(a,b);
    ifPass = RelErr < option.tol;
end


%% second formulation
%     cvx_begin sdp quiet
%     cvx_solver mosek       
%         variable gam_ti(1,1);
%         variable lam_ti(1,1)    nonnegative;
%         
%         minimize(-gam_ti);
%         
%         subject to
%             [-gam_ti Znx1'; Znx1 -lam_ti*Inx] ...
%                 + [0 -d'/2; -d/2 -As] >= 0
%     cvx_end
%     gam2 = gam_ti/lam_ti;
    
%     disp(gam);
%     disp(gam2);

% function [flag_Rank, cvx_status_RR, y, Y] = RankRelaxation(As, d, gam)
%     %% [Rank-relaxation formulation]
%     cvx_begin sdp quiet
%     cvx_solver mosek
%         variable y(nx);
%         variable Y(nx,nx) symmetric;
%         
%         minimize(-trace(Y));
%         
%         subject to 
%             trace(As*Y) + d'*y >= 0;
%             [Y y; y' 1] >=0;
%     cvx_end
%     cvx_status_RR = cvx_status;
%     
%     if strcmp(cvx_status_RR, 'Solved')   
%         flag_Rank = true;
%         
%         % [check solution] Lagrangian dual == Rank relaxation
%         assert(check_RelTol(gam, trace(Y), option), ...
%                 "Optimums of the two SDP are not equal");
%     end
% end