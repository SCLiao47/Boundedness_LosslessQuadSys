function plot_TwoState_quiver(ode, xlim, ylim, ifMag, scale)
% plot in the original coordinate

if nargin < 4
    ifMag = true;
    scale = 2;
end

%% options
num_grid = 40;

%% grid
[X, Y] = meshgrid(linspace(xlim(1), xlim(2), num_grid), ...
                  linspace(ylim(1), ylim(2), num_grid));
              
Xdot = X*0;
Ydot = X*0;

%% extract xdot
for i = 1:num_grid
    for j = 1:num_grid
        xdot = ode(0, [X(i,j); Y(i,j)]);
        
        % plot without magnitude if scale == -1
        if ifMag == false
            xdot = xdot./norm(xdot);
        end
        
        Xdot(i,j) = xdot(1);
        Ydot(i,j) = xdot(2);
    end
end

%% plot 
quiver(X, Y, Xdot, Ydot, scale);

end