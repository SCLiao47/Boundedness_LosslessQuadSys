% This script demonstrate verify the boundedness of Lorenze system

init;

%% setups

% SDP analysis options
option.round_Ndigit = 3;
option.verbose = true;
option.tol = 1e-6;

% if rerun scalability analysis. False for display-only
if_sim = true;

%% test KStackedLorenz implementation
disp("[KStackedLorenz test with K=2]");

K = 2;
model = model_KStackedLorenz(K);

% verify the existence of boundedness region
[m, info_m] = func_findNDShifting(model, option);
model_shifted = func_ShiftSystem(model, m);

%% Scalability studies
disp("[Scalability Test]");

if if_sim
    % random number generator for reproducibility
    seed = 0;
    rng(seed);
    
    % setup
    numK = 10;
    maxK = 350;
    
    gridK = ceil(logspace(0, log10(maxK), numK));
    gridK = unique(gridK);
    numK = length(gridK);
    
    numTest = 3;
    option.verbose = false;

    % memory
    T = zeros(numK, numTest);

    % Scalability test
    for idx_K = 1:numK              % for each state-dimension
        K = gridK(idx_K);
        fprintf("K = %04i: ", K);

        for i = 1:numTest       % repeate numTest of analysis for each dimension.
            fprintf("T%02i ", i);

            model = model_KStackedLorenz(K);

            % run and time the boundedness analysis
            % ONLY solve for coordinate?
            t0 = tic;
            [m, info_m] = func_findNDShifting(model, option);
            texc = toc(t0);
            
            if and(info_m.existTR, norm(m) == sqrt(38^2*K))
                T(idx_K,i) = texc;
            else
                warning("info_m.existTR = %i", info_m.existTR);
                warning("norm(m) = %f, sqrt(38^2*K) = ", norm(m), sqrt(38^2*K));
                T(idx_K,i) = nan;
            end
        end

        fprintf("\n")
    end

    % save the data
    save("Data/KStackedLorenz_RunTime", "gridK", "T", "seed");
else
    load("Data/KStackedLorenz_RunTime");
    
    [numK, numTest] = size(T);
end

%% plotting
figure(2); clf

average = sum(T,2)/numTest;

loglog(3*gridK, average, '*-', 'linewidth', 2);

title("Computation Time w.r.t. System Dimension");
xlabel("# of states");
ylabel("Execution Time (sec)");

grid on