function [model, R] = model_KStackedLorenz(K, toRotate)
% This model stack K Lorenz systems to get a large system with 3K states.
% The flag toRotate indicate if do a random coordinate rotation to couple
% the stacked systems. The boundedness property should remainds for the 
% stacked system.
%
% Note that the data in this example is sparse both constant, linear, and
% quadratic terms. Customized SDP solver could utilize the structure and
% accelerate the analysis accordingly.

switch nargin
    case 0
        K = 1;
        toRotate = false;
    case 1
        toRotate = false;
    case 2
        % nothing
    otherwise
        error(["Number of input is incorrect. Specify [K] for number of", ... 
        "stacked system and [toRotate] for if rotate the coordinate"]);
end

name = "K-Stacked_Lorenz";
nx = 3*K;

% get data of Lorenz system
Lorenz = model_Lorenz();

% create matrices for the stacked system
c = repmat(Lorenz.c, K, 1);

LCell = repmat({Lorenz.L}, K, 1);
L = blkdiag(LCell{:});

Q = zeros(nx, nx, nx);
for i = 1:K
    idx = (i-1)*3 + (1:3);
    Q(idx, idx, idx) = Lorenz.Q;
end

% Rotate the system
if toRotate
%     warning("Rotation TBI");

    % get randomized orthogonal matrix, R
    flag = false;
    while ~flag
        R = orth(rand(nx));
        
        flag = (rank(R) == nx);
    end

    c = R*c;
    L = R*L*R';

    % [Rotation of quadratic terms]
    % rotation of "xdot=...+phi(x)" to "xdot=...+phi(z)"
    for i = 1:nx
        Qi = Q(:, :, i);
        Q(:, :, i) = R*Qi*R';
    end
    
    % rotation of "xdot=...+phi(z)" to "zdot=...+phi(z)"
    Qtemp = zeros(nx, nx, nx);
    for i = 1:nx
        temp = zeros(nx, nx);
        for j = 1:nx
            temp = temp + R(i,j)*Q(:,:,j);
        end
        Qtemp(:, :, i) = temp;
    end
    Q = Qtemp;
end

model = class_Model_LosslessQuad(name, c, L, Q);
end