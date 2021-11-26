function x = rowforward(L, b)
    % This function is to be called as x = rowforward(L, b).
    % It solves a lower triangular system Lx = b, by row oriented forward substitution.
    
    [n, ~] = size(L);
    
    for k = 1: n
        for j = 1: k - 1
            b(k) = b(k) - L(k, j) * b(j);
        end
        
        if L(k, k) ~= 0
            b(k) = b(k) / L(k, k);
        else
            error("Matrix is Singular");
        end
    end
    
    x = b;
end