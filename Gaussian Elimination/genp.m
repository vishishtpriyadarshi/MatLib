function [L, U] = genp(A)
    % This function is to be called as [L, U] = genp(A).
    % It finds an LU factorization A = LU of an n-by-n matrix A by
    % performing Gaussian Elimination with no pivoting (GENP).
    
    [n, ~] = size(A);
    for k = 1: n - 1
        if A(k, k) ~= 0
            for i = k + 1: n
                A(i, k) = A(i, k) / A(k, k);
            end
        else
            error("Zero pivot encountered");
        end
        
        for i = k + 1: n
            for j = k + 1: n
                A(i, j) = A(i, j) - A(i, k) * A(k, j);
            end
        end
    end
    
    L = eye(n) + tril(A, -1);
    U = triu(A);
end