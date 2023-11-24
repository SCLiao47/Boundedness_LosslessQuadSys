% get inverse time model for unbounded SDP analysis
% 
% TODO
% - should be implemented in class for maintanability

function model_invT = get_inverseTimeModel(model)
    model_invT = model;
    
    % dynamics
    model_invT.c = - model.c;
    model_invT.L = - model.L;
    model_invT.Ls = - model.Ls;
    model_invT.Q = - model.Q;

    % utility functions
    model_invT.ode = @(t,x) -model.ode(t,x);
    model_invT.dK0 = @(x) -model.dK0(x);
    
    % if shifted dynamics
    if isfield(model, 'd')
        % dynamics
        model_invT.d = - model.d;
        model_invT.A = - model.A;
        model_invT.As = - model.As;
        
        % utility functions
        model_invT.ode = @(t,x) -model.ode_shifted(t,x);
        model_invT.dKm = @(x) -model.dKm(x);
    end
end