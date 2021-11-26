function [Q, R] = cgsrep(V)
    % Orthonormalize the columns of an mxn matrix by Classical Gram Schmidt
    % with reorthognalization so that Q is an isometry
    
    [~, m] = size(V);
    R = zeros(m);

    for k = 1 : m
%         for i = 1 : k-1
%             R(i, k) = V(:, i)' * V(:, k);
%             V(:, k) = V(:, k) - V(:, i) * R(i, k);
%         end
%         for i = 1 : k-1
%             alpha = V(:, i)' * V(:, k);
%             V(:, k) = V(:, k) - V(:, i) * alpha;
%             R(i, k) = R(i, k) + alpha;
%         end
        
        R(1 : k-1, k) = V(:, 1 : k-1)' * V(:, k);
        V(:, k) = V(:, k) - V(:, 1 : k-1) * R(1 : k-1, k);  % Orthogonalization
        s = V(:, 1 : k-1)' * V(:, k);
        
        V(:, k) = V(:, k) - V(:, 1 : k-1) * s;  % Re-orthogonalization
        R(1 : k-1, k) = R(1 : k-1, k) + s;

        R(k, k) = norm(V(:, k));
        if R(k, k) == 0
            error("Columns of V are dependent");
        end
        V(:, k) = V(:, k) / R(k, k);
    end
    Q = V;
end

