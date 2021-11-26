function G = mychol(A)
    % This function is to be called as G = mychol(A).
    % It executes the 'inner product form' of the Cholesky Decomposition for
    % finding the Cholesky factor of an n x n positive definite matrix A
    
    [n, ~] = size(A);
    G = zeros(n);
    
    for j = 1 : n
        G(j, j) = sqrt(A(j, j) - sum(G(1 : j - 1, j) .^ 2));
        for k = j + 1 : n
            G(j, k) = (A(j, k) - sum(G(1 : j - 1, j) .* G(1 : j - 1, k))) / G(j, j);
        end
    end
end