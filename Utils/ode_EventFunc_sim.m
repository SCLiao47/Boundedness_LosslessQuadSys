function [value, isterminal, direction] = ode_EventFunc_sim(t, x, Param)



% boundary event
xbound = Param.xbound;
ybound = Param.ybound;

dist = [x(1) - xbound(1), xbound(2) - x(1), ...
        x(2) - ybound(1), ybound(2) - x(2)];

value = min(dist);
isterminal = true;
direction = 0;
end