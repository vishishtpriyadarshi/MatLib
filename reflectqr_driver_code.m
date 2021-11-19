iter = input("Enter no of iterations = ");

for i = 1 : iter
    fprintf("\n\n===========  Iteration - %d  ===========\n", i);
    n = input("Enter n = ");
    m = input("Enter m = ");
    
    if n < m
        error("n < m not allowed");
    end
    
    A = rand(n, m);

    [Q, R] = reflectqr(A);
    [Qhat, Rhat] = qr(A, 0);

    fprintf("|| Q*R - A ||\t\t= ");
    disp(norm(Q*R - A));

    fprintf("|| Q'*Q - I || \t\t= ");
    disp(norm(Q'*Q - eye(m)));

    fprintf("tril(R, -1) \t\t= ");
    disp(norm(tril(R, -1)));

    fprintf("|| R - R_hat || \t= ");
    disp(norm(R - Rhat));

    fprintf("|| Q - Q_hat || \t= ");
    disp(norm(Q - Qhat));
end