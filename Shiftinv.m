function [iter, lambda] = Shiftinv(A, x, s, k)
    [n, ~] = size(A);
    A_prime = A - s * eye(n);
    [L, U, P] = lu(A_prime);
    
    iter = [x / max(abs(x))];
    for j = 1 : k
       b = P * iter(:, end);
       y = rowforward(L, b);
       q = colbackward(U, y);
        
       [~, idx] = max(abs(q));
       scaling_factor = q(idx);
       q = q / scaling_factor;
       
       iter(:, end + 1) = q;
    end
    
    lambda = s + 1/scaling_factor;
end