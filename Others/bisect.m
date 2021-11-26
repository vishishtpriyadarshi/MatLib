function x = bisect(p, x0, x1, tol)
    x = (x0 + x1) / 2;
    while x1 - x0 > 2*tol
        [a, b, c] = deal(Horner(p, x0), Horner(p, x), Horner(p, x1));
                
        if a * b < 0
            x1 = x;
        else
            x0 = x;
        end
        
        x = (x0 + x1) / 2;
    end
end

