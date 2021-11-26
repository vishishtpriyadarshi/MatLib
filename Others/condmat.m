function A = condmat(n, kappa)
    % Generate a random nxn positive definite matrix A with condition
    % number as kappa
    
    X = randn(n, n);
    [U, ~] = qr(X);
    val = kappa .^ ((1:n) ./ (n-1));
    D = diag(val);
    A = U * D * U';
end