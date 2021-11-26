function [Q, R] = mgs(V)
    % Orthonormalize the columns of an mxn matrix by Modified Gram Schmidt
    % so that Q is an isometry
    
    [~, m] = size(V);
    R = zeros(m);

    for k = 1 : m
        for i = 1 : k-1
            R(i, k) = V(:, i)' * V(:, k);
            V(:, k) = V(:, k) - V(:, i) * R(i, k);
        end
    
        R(k, k) = norm(V(:, k));
        if R(k, k) == 0
            error("Columns of V are dependent");
        end
        V(:, k) = V(:, k) / R(k, k);
    end
    Q = V;
end

