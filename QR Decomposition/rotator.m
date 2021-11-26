function [c, s, v]= rotator(xi, xj)
    % Q = [c -s; s c]    x = [xi xj]
    
    if abs(xj) <= abs(xi)
        t = xj / xi;
        c = sign(xi) / sqrt(1 + t*t);
        s = -c * t;
        v = xi * sqrt(1 + t*t);
    else
        k = xi / xj;
        s = -sign(xj) / sqrt(1 + k*k);
        c = -s * k;
        v = xj * sqrt(1 + k*k);
    end
end

