function [x, residual] = lsp(A, b)
    % Solve overdetermined system Ax = b when A has full rank.
    
    [n, m] = size(A);
    [~, R, gamma] = reflectqr_modified(A);
    c = b;
    
    for i = 1 : m
        u = ones(n - i + 1, 1);
        u(2 : end) = R(i + 1 : n, i);
        c(i : n) = appreflect(u, gamma(i), c(i : n));
    end
    
    R = triu(R);
    R1 = R(1 : m, :);
    x = colbackward(R1, c(1 : m));
    residual = norm(c(m + 1 : end));
end


function [Q, R, gamma_val] = reflectqr_modified(A)
    [n, m] = size(A);
    gamma_val = zeros(n, 1);
    Q = eye(n);
    
    iter = m;
    if m == n
        iter = m - 1;
    end
    
    for k = 1 : iter
       [u, gamma, tau] = reflect(A(k : n, k));
        A(k : n, k+1 : m) = appreflect(u, gamma, A(k : n, k+1 : m));
        A(k : n, k) = u;
        A(k, k) = -tau;
        gamma_val(k) = gamma;
        
        for i = k : -1 : 1
            u = ones(n - i + 1, 1);
            u(2 : end) = A(i + 1 : n, i);
            Q(i : n, k) = appreflect(u, gamma_val(i), Q(i : n, k));
        end
    end
    
    for k = iter + 1 : n
        for i = m : -1 : 1
            u = ones(n - i + 1, 1);
            u(2 : end) = A(i + 1 : n, i);
            Q(i : n, k) = appreflect(u, gamma_val(i), Q(i : n, k));
        end
    end
    R = A;
end