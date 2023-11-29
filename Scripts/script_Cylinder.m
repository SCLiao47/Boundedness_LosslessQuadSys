% This script analyze the boundedness of various Cynlinder model

init;

%% setups

% options
option.round_Ndigit = 3;
option.verbose = true;
option.tol = 1e-6;

%% Mean-field model
model_MF = model_Cylinder_MeanField();

[isBounded_MF, info_MF] = func_boundedness_SDP(model_MF, option);

%% 3-mode model from data
model_3mode = model_Cylinder_Noack(3);

[isBounded_3mode, info_3mode] = func_boundedness_SDP(model_3mode, option);


%% 8-mode
model_8mode = model_Cylinder_Noack(8);

[isBounded_8mode, info_8mode] = func_boundedness_SDP(model_8mode, option);

% CHECK WHY AS>0 found but no ER

%% 9-mode
model_9mode = model_Cylinder_Noack(9);

[isBounded_9mode, info_9mode] = func_boundedness_SDP(model_9mode, option);

% FIX OUTPUT of the analysis with REGULAIRZATION