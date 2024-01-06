function [isBounded, info] = func_boundedness_SDP(model, options)
% [isBounded, info] = func_boundedness_SDP(model, options)
%
% DESCRIPTION
%   This functino analyze the boundedness of given model using the SDP
%   approach based on S&N Trapping Theorem.
%

%% check option
if nargin < 2
    %         options.round_Ndigit = 6;
    options.round_Ndigit = nan;
    
    options.verbose = true;
    options.tol = 1e-6;
end

%% initialize output
isBounded = [];
info = struct('existTR', [], 'mND', [], 'info_mND', [], ...
    'r_SDP', [], 'info_TR', [], 'r_SN', []);

%% verify the existence of boundedness region
[mND, info_mND] = func_findNDShifting(model, options);

info.mND = mND;
info.existTR = info_mND.existTR;
info.info_mND = info_mND;

if info.existTR
    % TR exists.
    isBounded = true;
    
    % Solve for the size of boundedness region
    model_shifted = func_ShiftSystem(model, mND);
    
    % proposed SDP approach
    [r_SDP, info_TR] = func_TRSize_SDP(model_shifted, options);
    
    % Schlegel and Noack's approach
    r_SN = func_TRSize_SN(model_shifted, options);
    
    if options.verbose
        disp('Trapping region exists. The model is bounded.');
    end
    
    % setting output
    info.r_SDP = r_SDP;
    info.info_TR = info_TR;
    info.r_SN = r_SN;
else
    % there is no trapping region exists
    
    if options.verbose
        disp('No trapping region exists. Futher analysis is needed to conclude boundedness.');
    end
end
end
