function [output] = config()
%CONFIG Summary of this function goes here
%   Detailed explanation goes here


size_label = 10;
size_legend = 10;
size_tick = 10;
output.plotting = struct('size_label', size_label, ...
                         'size_legend', size_legend, ...
                         'size_tick', size_tick ...
                        );


end

