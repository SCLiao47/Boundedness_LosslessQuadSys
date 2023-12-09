% This script demonstrate verify the boundedness of Lorenze system

init;

%% setups

% SDP analysis options
option.round_Ndigit = nan;
option.verbose = false;
option.tol = 1e-6;

% if rerun scalability analysis. False for display-only
if_sim = true;

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
    numK = 3;
    maxK = 100;
    numTest = 1;
    
    gridK = ceil(logspace(0, log10(maxK), numK));
    gridK = unique(gridK);
    numK = length(gridK);

    % memory
    T = zeros(numK, numTest);

    % Scalability test
    for idx_K = 1:numK              % for each state-dimension
        K = gridK(idx_K);
        fprintf("K = %04i: ", K);

        for i = 1:numTest       % repeate numTest of analysis for each dimension.
            fprintf("T%02i ", i);

            model = model_KStackedLorenz(K, true);

            % run and time the boundedness analysis
            % ONLY solve for coordinate?
            t0 = tic;
            [m, info_m] = func_findNDShifting(model, option);
            texc = toc(t0);
            
            % validation 
            if and(info_m.existTR, abs(info_m.a - -1) < option.tol)
                T(idx_K,i) = texc;
            else
                warning("Validation step failed, check");
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