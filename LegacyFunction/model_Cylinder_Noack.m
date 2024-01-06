function model = model_Cylinder_Noack(numState)
% return the n-state Cylinder model constructe from the data shared by
% Prof. Noack.

% load data
data = load('Model\Noack2003Cylinder_Refined.mat');

name = ["Cylinder_", num2str(numState), "mode"];

switch numState
    case 3
        c = data.c3;
        L = data.L3;
        Q = data.Q3;
        
    case 8
        c = data.c8;
        L = data.L8;
        Q = data.Q8;
        
    case 9
        c = data.c9;
        L = data.L9;
        Q = data.Q9;
        
    otherwise
        error([num2str(numState), "-state model not implemented"]);
end

model = class_Model_LosslessQuad(name, c, L, Q);
end