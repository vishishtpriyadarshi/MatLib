function [H, Q] = hessenberg(A)
    % Find upper hessenberg matrix H such that Q*AQ = H
    
    [n, ~] = size(A);
    gamma_val = zeros(n, 1);
    Q = eye(n);
     
    for k = 2 : n
       [u, gamma, tau] = reflect(A(k : n, k-1));
        A(k : n, k : n) = appreflect(u, gamma, A(k : n, k : n));
        A(:, k : n) = appreflect_modified(u, gamma, A(:, k : n));
        A(k : n, k-1) = u;
        A(k, k-1) = -tau;
       gamma_val(k) = gamma;
        
       for i = k : -1 : 2
            u = ones(n - i + 1, 1);
            u(2 : end) = A(i + 1 : n, i - 1);
            Q(i : n, k) = appreflect(u, gamma_val(i), Q(i : n, k));
       end
    end
    H = triu(A, -1);
end

% BUG - Last col gets multiplied by (-1)