function W = Wilkinson(n)
    % This function is to be called as W = Wilkinson(n). 
    % It creates a square Wilkinson matrix without using for loops.
    
    W = (tril(ones(n)) .* -1) + (eye(n) .* 2);
    W(:, end) = 1;
end