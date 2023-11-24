function model_shifted = func_ShiftSystem(model, m)

    % initialize model
    model_shifted = model;
    
    nx = model.nx;
    c = model.c;
    L = model.L;
    Q = model.Q;
        
    % update constant part
    d = zeros(nx,1);
%     d = sym('d',[nx,1],'real');
    for i = 1:nx
        d(i) = c(i) + L(i,:)*m + m'*Q(:,:,i)*m;
    end
    
    % update linear part
    A = zeros(nx,nx);
%     A = sym('A',[nx,nx],'real');
    for i = 1:nx
%         for j = 1:nx
%             A(i,j) = L(i,j) +  2*m'*Q(:,j,i);
%         end
        
        A(i,:) = L(i,:) + 2*m'*Q(:,:,i);
    end
    As = 1/2*(A+A');
    
    As1 = 1/2*(L+L');
    for i = 1:nx
        As1 = As1 - m(i)*Q(:,:,i);
    end
    
    % quadratic part stays the same
    
    % update model class object
    model_shifted.m = m;
    
    model_shifted.d = d;
    model_shifted.A = A;
    model_shifted.As = As;
    
    model_shifted.ode_shifted = @(t,x) ode_quadraticDyn(d,A,Q,t,x);
    model_shifted.dKm = @(x) d'*x + x'*As*x;
end