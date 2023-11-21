function [rTrap, info]= func_TRSize_SN(model, option)
% estimate the size of trapping region by Schlegel and Noack's method:
% worst-case spectral analysis

    lam1 = eigs(model.As, 1, 'lr');
    
    if lam1 >= 0 
        % TR not exists, set rTrap = inf
        rTrap = inf;
    else
        % exists Trapping region. 
        rTrap = norm(model.d) / abs(lam1);
    end
    
    info.lam1 = lam1;
end