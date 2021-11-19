function [c, s, v]= rotator(xi, xj)
    % Q = [c -s; s c]    x = [xi xj]
    
    beta = max(abs(xi), abs(xj));
    if beta == 0
        [c, s, v] = deal(1, 0, 0);
    else
        [xi, xj] = deal(xi/beta, xj/beta);
        v = sqrt(xi^2 + xj^2);
        [c, s] = deal(xi/v, xj/v);
        v = v * beta;
    end
end

