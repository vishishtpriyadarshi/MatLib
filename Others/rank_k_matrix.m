function [A, rk] = rank_k_matrix(m, n, k)
    % create a mxn matrix with rank = k
    
    P = orth(randn(m, k));
    Q = orth(randn(n, k))';
    A = P * Q;
    rk = rank(A);
end