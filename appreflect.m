function B = appreflect(u, gamma, A)
    % Performs multiplication QA where, Q = I - gamma * u * u'
    
    v = gamma * u;
    w = u' * A;
    C = v * w;
    B = A - C;
end