function plot_TwoState_Kdot_alongLine(model_shifted, V, scale)

num_th = 100;

%% computing
th_alongV = linspace(-1.2*scale, 1.2*scale, num_th);

Kdot_alongV = 0*th_alongV;
for i = 1:num_th
    Kdot_alongV(i) = model_shifted.dKm(V * th_alongV(i));
end

%% plotting
plot(th_alongV, Kdot_alongV, 'linewidth', 2)

grid on;
hold on;

rect_x = [-scale, scale, scale, -scale];
rect_y = [min(Kdot_alongV)*[1,1], max(Kdot_alongV)*[1,1]];
ph_r = patch(rect_x, rect_y, 'red','FaceAlpha',.3);

end