function [L, U, p] = gepp(A)
    % This function is to be called as [L, U, p] = gepp(A).
    % It finds an LU factorization A = LU of an n-by-n matrix A and a column vector p satisfying
    % A(p, :) = LU via Gaussian Elimination with partial pivoting (GEPP).
    
    [n, ~] = size(A);
    p = (1: n);
    for k = 1: n - 1
        [~, idx] = max(abs(A(k: end, k)));
        idx = idx + k - 1;
        
        if idx ~= k
            permutation = [idx, k];
            A(permutation, :) = A(permutation([2, 1]), :);
            p(permutation) = p(permutation([2, 1]));
        end
        
        if A(k, k) ~= 0
            A(k + 1: n, k) = A(k + 1: n, k) / A(k, k);
        end
        
        A(k + 1: n, k + 1: n) = A(k + 1: n, k + 1: n) - A(k + 1: n, k) * A(k, k + 1: n);
    end
    
    L = eye(n) + tril(A, -1);
    U = triu(A);
end