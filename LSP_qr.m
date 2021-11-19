function x = LSP_qr(A, b)
    [Q, R, p, gamma, r] = reflectqr_rank_revealing(A);
    c = b;
    
    for i = 1 : r
        c = appreflect(c, gamma(i), Q(:, i));
    end
    disp(norm(A * p - Q*R));
    
    [c1, c2] = deal(c(1 : r), c(r + 1 : end));
    residual = norm(c2);
    y = rowforward(R(1 : r, 1 : r), c1);
    x = p * y;
    
    disp(residual);
end

