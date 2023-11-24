classdef class_Model_LosslessQuad < matlab.mixin.Copyable
    properties
        name    % name of the model
        nx      % number of states
        
        % dynamics at the origin: xdot = c + Lx + phi(x)
        c       % constant dynamics
        L       % linear dynamics
        Ls      % symmetric part of linear dynamics
        La      % asymmetric part of linear dynamics
        Q       % quadratic dynamics
        
        % dynamics at the coordinate shift ydot = d + Ay + phi(y)
        m       % shifted coordinate, default to be 0. 
        d       % constant dynamics
        A       % linear dynamics
        As      % symmetric part of linear dynamics
        Aa      % asymmetric part of linear dynamics
        
        % 
        ode
        dK0 
        ode_shifted
        dKm
    end
    
    methods
        function obj = class_Model_LosslessQuad(name, c, L, Q, m)
            obj.name = name;
            obj.nx = size(c,1);
            
            obj.c = c;
            obj.L = L;
            obj.Ls = 1/2*(L+L');
            obj.La = L - obj.Ls;
            obj.Q = Q;
            
            obj = setup_dyn(obj);
            
            if nargin ~= 5
                m = zeros(obj.nx,1);
            end
            obj = func_ShiftSystem(obj, m);                
        end
        
        function obj = setup_dyn(obj)
            obj.ode = @(t,x) ode_quadraticDyn(obj.c, obj.L, obj.Q, t, x);
            obj.dK0 = @(x) obj.c'*x + x'*obj.Ls*x;
        end
        
        function obj = func_ShiftSystem(obj, m)
            obj.m = m;
            
            % update constant part
            obj.d = zeros(obj.nx,1);
            %     d = sym('d',[nx,1],'real');
            for i = 1:obj.nx
                obj.d(i) = obj.c(i) + obj.L(i,:)*m + m'*obj.Q(:,:,i)*m;
            end
           
            % update linear part
            obj.A = zeros(obj.nx,obj.nx);
            %     A = sym('A',[nx,nx],'real');
            for i = 1:obj.nx
                obj.A(i,:) = obj.L(i,:) + 2*m'*obj.Q(:,:,i);
            end
            
            obj.As = 1/2*(obj.A+obj.A');
            %             % alternative:
            %             obj.As = 1/2*(obj.L+obj.L');
            %             for i = 1:obj.nx
            %                 obj.As = obj.As - m(i)*obj.Q(:,:,i);
            %             end
            obj.Aa = obj.A - obj.As;
            
            % update ode and power function
            obj.ode_shifted = @(t,x) ode_quadraticDyn(obj.d, obj.A, obj.Q, t, x);
            obj.dKm = @(x) obj.d'*x + x'*obj.As*x;
        end
        
        function model_invTime = get_inverseTimeModel(obj)
            cInv = -obj.c;
            LInv = -obj.L;
            QInv = -obj.Q;
            
            model_invTime = class_Model_LosslessQuad(obj.name, ...
                cInv, LInv, QInv, obj.m);
        end
    end
end