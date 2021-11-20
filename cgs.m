function [Q, R] = cgs(V)
    % Orthonormalize the columns of an mxn matrix by Classical Gram Schmidt
    % so that Q is an isometry
    
    [~, m] = size(V);
    R = zeros(m);

    for k = 1 : m
        R(1 : k-1, k) = V(:, 1 : k-1)' * V(:, k);
        V(:, k) = V(:, k) - V(:, 1 : k-1) * R(1 : k-1, k);
    
        R(k, k) = norm(V(:, k));
        if R(k, k) == 0
            error("Columns of V are dependent");
        end
        V(:, k) = V(:, k) / R(k, k);
    end
    Q = V;
end