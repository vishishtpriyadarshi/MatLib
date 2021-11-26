function x = geppsolve(A, b)
    % This function is to be called as x = geppsolve(A, b).
    % It solves the system Ax = b via GEPP.
    
    [L, U, p] = gepp(A);
   
    y = rowforward(L, b(p, :));
    x = colbackward(U, y);
end