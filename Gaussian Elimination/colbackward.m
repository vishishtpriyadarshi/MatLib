function x = colbackward(U, b)
    % This function is to be called as x = colbackward(U, b).
    % It solves a upper triangular system Ux = b, by column oriented back substitution.
    
    [n, ~] = size(U);
    
    for i = n: -1: 1
        if U(i, i) ~= 0
            b(i) = b(i) / U(i, i);
        else
            error("Matrix is Singular");
        end
     
        for j = i - 1: -1: 1
            b(j) = b(j) - U(j, i) * b(i);
        end
    end
   
    x = b;
end