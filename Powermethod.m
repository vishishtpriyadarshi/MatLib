function [iter, lambda] = Powermethod(A, x, k)
    iter = [x / max(abs(x))];
    
    for j = 1 : k
       q = A * iter(:, end);
       
       [~, idx] = max(abs(q));
       scaling_factor = q(idx);
       q = q / scaling_factor;
       
       iter(:, end + 1) = q;
    end
    
    lambda = scaling_factor;
end