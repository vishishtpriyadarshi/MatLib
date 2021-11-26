function [u, gamma, tau] = reflect(x)
    beta = max(abs(x));
    if beta == 0
        [gamma, tau] = deal(0);
    else
       x = x ./ beta;
       tau = norm(x);
       if x(1) < 0
           tau = -tau;
       end
       x(1) = tau + x(1);
       gamma = x(1) / tau;
       x(2 : end) = x(2 : end) / x(1);
       x(1) = 1;
       tau = tau * beta;
    end
    u = x;
end

