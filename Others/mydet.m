function d = mydet(A)
    % This function is to be called as d = mydet(A).
    % It computes the determinant of the matrix in O(n^3) flops using an 
    % efficient version of LU factorisation 
    
    [~, U, ~, sign] = gepp_modified(A);
    [n, ~] = size(A);
    d = sign;
    
    for i = 1: n
        d = d * U(i, i);
    end
end