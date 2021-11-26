function H = Hamiltonian(n)
    % This function is to be called as H = Hamiltonian(n). 
    % It creates a random real 2n X 2n Hamiltonian matrix.
    
    H_11 = randn(n);
    H_12 = randn(n);
    H_12 = tril(H_12) + tril(H_12)' - eye(n) .* diag(H_12);
    
    H_21 = randn(n);
    H_21 = tril(H_21) + tril(H_21)' - eye(n) .* diag(H_21);
    
    H = [H_11, H_12; H_21, -1 .* H_11'];
end