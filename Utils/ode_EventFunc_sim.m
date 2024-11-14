function [value, isterminal, direction] = ode_EventFunc_sim(t, x, Param)

switch size(x,1)
    case 2
        % boundary event
        xbound = Param.xbound;
        ybound = Param.ybound;

        dist = [x(1) - xbound(1), xbound(2) - x(1), ...
                x(2) - ybound(1), ybound(2) - x(2)];
    case 3
        xbound = Param.xbound;
        ybound = Param.ybound;
        zbound = Param.zbound;
        
        dist = [x(1) - xbound(1), xbound(2) - x(1), ...
                x(2) - ybound(1), ybound(2) - x(2), ...
                x(3) - zbound(1), zbound(2) - x(3)];
end
        
value = min(dist);
isterminal = true;
direction = 0;
end