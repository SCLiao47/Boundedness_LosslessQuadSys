% This script demonstrate verify the boundedness of Lorenze system

init;

%% setups

% SDP analysis options
option.round_Ndigit = nan;
option.verbose = false;
option.tol = 1e-4;
 
% if rerun scalability analysis. False for display-only
if_sim = ~true;

%% test KStackedLorenz implementation
disp("[KStackedLorenz test with K=2]");

K = 2;
model = model_KStackedLorenz(K, true);

% verify the existence of boundedness region
[m, info_m] = func_findNDShifting(model, option);
model_shifted = func_ShiftSystem(model, m);

assert(abs(info_m.a - -1) < option.tol); 

%% Scalability studies
disp("[Scalability Test]");

if if_sim
    % random number generator for reproducibility
    seed = 0;
    rng(seed);
    
    % setup
    numK = 20;
    maxK = 300;
    numTest = 10;
    
    gridK = ceil(logspace(0, log10(maxK), numK));
    gridK = unique(gridK);
    numK = length(gridK);

    % memory
    T_shift = zeros(numK, numTest);
    T_size = zeros(numK, numTest);

    % Scalability test
    for idx_K = 1:numK              % for each state-dimension
        K = gridK(idx_K);
        fprintf("K = %04i: ", K);

        for i = 1:numTest       % repeate numTest of analysis for each dimension.
            fprintf("T%02i ", i);

            model = model_KStackedLorenz(K, true);

            % run and time the boundedness analysis
            t0 = tic;
            [m, info_m] = func_findNDShifting(model, option);
            texc_shfit = toc(t0);
            
            model_shift = func_ShiftSystem(model, m);
            
            t1 = tic;
            [r, info_TR] = func_TRSize_SDP(model_shifted, option);
            texc_size = toc(t1);
            
            % validation 
            val_shift = and(info_m.existTR, abs(info_m.a - -1) < option.tol);
            val_size = info_TR.feasibility;
            
            if and(val_shift, val_size)
                T_shift(idx_K, i) = texc_shfit;
                T_size(idx_K, i) = texc_size;
            else
                warning("Validation step failed, check");
                T_shift(idx_K, i) = nan;
                T_size(idx_K, i) = nan;
            end
        end

        fprintf("\n")
    end

    % save the data
    save("Data/KStackedLorenz_RunTime", "gridK", "T_shift", "T_size", "seed");
else
%     load("Data/KStackedLorenz_RunTime");

    load("Data/KStackedLorenz_RunTime_1_223.mat");
    idx_remove = 19;                % This set of K was out of memories
    T_shift(idx_remove,:) = [];
    T_size(idx_remove,:) = [];
    gridK(idx_remove) = [];
    
    [numK, numTest] = size(T_shift);
end

%% Data analysis
T_total = T_shift + T_size;

% statistics
mean_T = mean(T_total');
std_T = std(T_total');

% Complexity analysis via asymptote for large K (3K>100):
%   T = c* (3K)^u => log(T) = log(c) + u*log(3K)
%
% Linear regression on x = [log(c); u]:
%   Ax = b 
% with 
%   A = [1 log(3K_i); ...; 1 log(3K_end)];
%   b = [log(T_i); ... log(T_end)];

idx_asy = find(3*gridK > 150);
num_asy = size(idx_asy,2);

A_asy = [ones(num_asy,1) log(3*gridK(idx_asy)')];
b_asy = log(mean_T(idx_asy)');

x_asy = A_asy\b_asy;
c_asy = exp(x_asy(1));
u_asy = x_asy(2);

%% plotting
Plot_complexity = figure(2); clf
Plot_complexity.WindowStyle = 'normal';
Plot_complexity.DockControls = 'off';
width = 700;            % for IJRNC format, almost fit the 1 column
height = 350;
set(gcf, 'units', 'pixels', 'position', [1, 1, width, height]);
set(gca, 'XScale', 'log', 'YScale','log');
hold on;

% standard deviation
ph_std = fill([3*gridK, flip(3*gridK)], [mean_T+std_T, flip(mean_T-std_T)], 0.7*[1 1 1]);

% mean
ph_mean = loglog(3*gridK, mean_T, '*-', 'linewidth', 2, 'color', 'b');

% asymptote
grid_asy = [60, 900];
ph_asy = loglog(grid_asy, c_asy*(grid_asy).^u_asy, '--', 'linewidth' , 2, ...
                                        'color', 'r');

formatSetting = {'interpreter', 'latex', 'fontsize', 16};

text_asy =sprintf("log(T)=%.2flog(n)%.2f", u_asy, log(c_asy));

% title(sprintf("Computation Time w.r.t. System Dimension. log(T) = %.2f + %.2f * log(n)", log(c_asy), u_asy), ...
%       formatSetting{:});
xlabel("State Dimension, $n=3K$", formatSetting{:});
ylabel("Execution Time, $T$ (sec)", formatSetting{:});
legend([ph_mean, ph_std, ph_asy], {'Mean'; 'Standard Deviation'; text_asy}, ...
       'location', 'northwest', formatSetting{:}, 'fontsize', 12);

box on;
grid on;

print('Figure/StackedLorenz_TimeComplexity','-depsc')