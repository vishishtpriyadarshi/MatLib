function [x, residual] = lsp_qr(A, b)
    % Solve overdetermined system Ax = b.
    
    [n, m] = size(A);
    [~, R, p, gamma, r] = reflectqr_rank_revealing_modified(A);
    c = b;
    
    for i = 1 : m
        u = ones(n - i + 1, 1);
        u(2 : end) = R(i + 1 : n, i);
        c(i : n) = appreflect(u, gamma(i), c(i : n));
    end
    
    R = triu(R);
    R1 = R(1 : r, 1 : r);
    [c1, c2] = deal(c(1 : r), c(r + 1 : end));
    
    residual = norm(c2);
    y = zeros(m, 1);
    y(1 : r) = rowforward(R1, c1);
    x = p * y;
end

function [Q, R, p, gamma_val, r] = reflectqr_rank_revealing_modified(A)
    [n, m] = size(A);
    gamma_val = zeros(n, 1);
    Q = eye(n);
    
    iter = m;
    if m == n
        iter = m - 1;
    end
    
    col_norm = zeros(m, 1);
    for i = 1 : m
        col_norm(i) = A(:, i)' * A(:, i);
    end
    
    p = (1 : m);
    for k = 1 : iter      
        [~, idx] = max(col_norm(k : end));
        idx = idx + k - 1;
        permutation = [idx, k];
        
        A(:, permutation) = A(:, permutation([2, 1]));
        p(permutation) = p(permutation([2, 1]));       
        
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
        
        for i = k + 1 : m
            if p(i) == i
                col_norm(i) = col_norm(i) - A(k, i) ^ 2;
            else
                col_norm(i) = col_norm(k) - A(k, k) ^ 2;
            end
        end
    end
    
    for k = iter + 1 : n
        for i = m : -1 : 1
            u = ones(n - i + 1, 1);
            u(2 : end) = A(i + 1 : n, i);
            Q(i : n, k) = appreflect(u, gamma_val(i), Q(i : n, k));
        end
    end
    
    % R = triu(A)
    R = A;
    I = eye(m);
    p = I(:, p);
    
    [r, tolerance] = deal(m, n * abs(R(1, 1)) * eps(1));
    for i = 2 : m
        if abs(R(i, i)) < tolerance
            r = i - 1;
            break;
        end
    end
    % R(r + 1 : end, :) = zeros(n - r, m);
end