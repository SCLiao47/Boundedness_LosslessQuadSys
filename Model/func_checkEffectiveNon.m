function [isEN, flag] = func_checkEffectiveNon(L, Q)

    nx = size(L, 1);
    [V, ~] = eig(L);
    
    flag = zeros(nx,nx);    % each row is for a given Nonlinearity
    for idx_v = 1:nx
        v = V(:, idx_v);
        for idx_Q = 1:nx
            flag(idx_Q, idx_v) = v'*Q(:,:,idx_Q)*v;
        end
    end

    if all(any(abs(flag) > 0))
        isEN = true;
    else
        isEN = false;
    end

end