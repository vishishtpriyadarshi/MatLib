function [Q, R, p, gamma_val, m] = reflectqr_rank_revealing(A)
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
    
    R = triu(A);
    I = eye(m);
    p = I(:, p);
end