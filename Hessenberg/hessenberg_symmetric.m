function [H, Q] = hessenberg_symmetric(A)
    % Find upper hessenberg matrix H such that Q*AQ = H of a symmetric
    % matrix efficiently
    
    [n, ~] = size(A);
    gamma_val = zeros(n, 1);
    Q = eye(n);
     
    for k = 2 : n
       [u, gamma, tau] = reflect(A(k : n, k-1));
        A(k : n, k : n) = appreflect_symmetric(u, gamma, A(k : n, k : n));
      
        A(k : n, k-1) = u;
        A(k, k-1) = -tau;
        
        A(k-1, k : n) = u;
        A(k-1, k) = -tau;
        gamma_val(k) = gamma;
        
        for i = k : -1 : 2
            u = ones(n - i + 1, 1);
            u(2 : end) = A(i + 1 : n, i - 1);
            Q(i : n, k) = appreflect(u, gamma_val(i), Q(i : n, k));
        end
    end
    H = triu(A, -1);
    H = triu(H', -1)';
end


function B = appreflect_symmetric(u, gamma, A)
    v = -gamma * A * u;
    alpha = -0.5 * gamma * (u' * v);
    w = v + alpha * u;
    B = A + w * u' + u * w';
end