function [Q, R] = qr_update(Q, R, u, v)
    % Find R of QR decomposition of M = QR + uv' in O(n2) flops
    
    n = size(Q, 1);
    w = Q' * u;
     
    for k = n - 1 : -1 : 1
        [c, s, r] = rotator(w(k), w(k+1));
        [w(k), w(k+1)] = deal(r, 0);
        J = [c -s; s c];
        
        R(k : k+1, :) = J * R(k : k+1, :);
        Q(:, k : k+1) = Q(:, k : k+1) * J';
    end
    
    R(1, :) = R(1, :) + w(1)*v';
    [Q1, R1] = hessqr(R);
    [Q, R] = deal(Q * Q1, R1);
end