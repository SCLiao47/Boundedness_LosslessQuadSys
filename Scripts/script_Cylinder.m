% This script analyze the boundedness of various Cynlinder model

init;

%% setups

% options
option.round_Ndigit = 3;
option.verbose = true;
option.tol = 1e-6;

%% Mean-field model
fprintf("\n===[ Mean-field Model ]===\n");

model_MF = model_Cylinder_MeanField();

[isBounded_MF, info_MF] = func_boundedness_SDP(model_MF, option);

report_Boundedness(isBounded_MF, info_MF)

%% 3-mode model from data
fprintf("\n===[ 3-mode Model ]===\n");

model_3mode = model_Cylinder_Noack(3);

[isBounded_3mode, info_3mode] = func_boundedness_SDP(model_3mode, option);

report_Boundedness(isBounded_3mode, info_3mode)

%% 8-mode
fprintf("\n===[ 8-mode Model ]===\n");

model_8mode = model_Cylinder_Noack(8);

[isBounded_8mode, info_8mode] =     func_boundedness_SDP(model_8mode, option);

report_Boundedness(isBounded_8mode, info_8mode)

% CHECK WHY AS>0 found but no ER

%% 9-mode
fprintf("\n===[ 9-mode Model ]===\n");

model_9mode = model_Cylinder_Noack(9);

[isBounded_9mode, info_9mode] = func_boundedness_SDP(model_9mode, option);

report_Boundedness(isBounded_9mode, info_9mode)

% FIX OUTPUT of the analysis with REGULAIRZATION
% FIX output SDP size


%% utilities
function report_Boundedness(isBounded, info)
    if isBounded
        disp("The system is BOUNDED");
        
        if info.existTR
            disp("A Trapping Region B(m,R*) exists with:");
            disp("   m = [" + regexprep(num2str(info.mND'),' +',' ;') + "]");
            disp("   R* = " + num2str(info.r_SDP));
        end
    else
        disp("The system is UNBOUNDED");
    end
end