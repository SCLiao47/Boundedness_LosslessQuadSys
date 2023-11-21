function [xdot] = ode_quadraticDyn(c, L, Q, t, x)
%% ODE for quadratic nonlinear system given a model
%
% Given the system:
%   xdot = c + Lx + [x'Qix]_i^n
% where c is the vector of affine dynamics, L is linear dynamics, and Q is
% a tensor with Qi being the quardic dynamics for i-th state. 
%
% INPUTS:
%   c:  	affine dynamics
%   L:      linear dynamics
%   Q:      quadratic dynamics
%   t:      current time
%   x:      current state
%
% OUTPUT:
%   xdot:   time-derivative of states given by system dynamics
%
% ---
% Shih-Chi Liao, 2023

%% compute dynamics
Qterms = 0*x;

for i = 1:length(x)
    Qterms(i) = x'*Q(:,:,i)*x;
end

xdot = c + L*x + Qterms;
end
