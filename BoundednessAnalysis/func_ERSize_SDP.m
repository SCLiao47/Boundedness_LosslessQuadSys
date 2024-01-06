function [rER, info]= func_ERSize_SDP(model, option)

    % get inverse-time model
    model_invT = get_inverseTimeModel(model);
    
    % do trapping region analysis on inverse-time model
    disp('==== Inverse-time study: START ====');
    [rTR, info_TR] = func_TRSize_SDP(model_invT, option);
    disp('==== Inverse-time study: END ====');
    
    % convert back to forward-time model
    info.feasibility = info_TR.feasibility;
    
    if info.feasibility
        rER = rTR;
        
        % dual solution
        info.gam = info_TR.gam;
        info.lam = info_TR.lam;
        
        % primal solution
        info.ystar = info_TR.ystar;
        info.rank = info_TR.rank;
    end
end