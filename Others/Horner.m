function y = Horner(p, x)
    n = length(p);
    y = zeros(size(x));
    
    for i = 1 : n
        y = p(i) + x .* y;
    end
end

