function [Q, R] = hessqr(A)
    % Perform QR decomposition of upper hessenberg matrix H in O(n2) flops
    
    [n, ~] = size(A);
    gamma_val = zeros(n, 1);
    Q = eye(n);
    
    for k = 1 : n - 1
       [u, gamma, tau] = reflect(A(k : k+1, k));
        A(k : k+1, k+1 : n) = appreflect(u, gamma, A(k : k+1, k+1 : n));
        A(k : k+1, k) = u;
        A(k, k) = -tau;
        gamma_val(k) = gamma;
        
        for i = k : -1 : 1
            u = ones(2, 1);
            u(2) = A(i + 1, i);
            Q(i : i+1, k) = appreflect(u, gamma_val(i), Q(i : i+1, k));
        end
    end
    
    for i = n-1 : -1 : 1
        u = ones(2, 1);
        u(2) = A(i + 1, i);
        Q(i : i+1, n) = appreflect(u, gamma_val(i), Q(i : i+1, n));
    end
    
    R = triu(A);
end