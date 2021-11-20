function [iter, lambda] = Rayleigh(A, x, k)
    [n, ~] = size(A);
    [~, H] = hess(A);
    
    [~, idx] = max(abs(x));
    iter = [x / x(idx)];
    
    q0 = iter(:, 1);
    rho = q0' * H * q0 / (q0' * q0);
    
    for j = 1 : k
        A_prime = H - rho * eye(n);
        [L, U, P] = lu(A_prime); 
        b = P * iter(:, end);
        y = rowforward(L, b);
        q = colbackward(U, y);
        
        [~, idx] = max(abs(q));
        scaling_factor = q(idx);
        q = q / scaling_factor;
       
        iter(:, end + 1) = q;
        rho = q' * H * q / (q' * q);
    end
    
    lambda = rho + 1/scaling_factor;
end