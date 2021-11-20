format long e;
x = [1 1 1]';

array = zeros(3, 3, 3);
array(:, :, 1) = [1 1 1; -1 9 2; 0 -1 2];
array(:, :, 2) = [1 1 1; -1 9 2; -4 -1 2];
array(:, :, 3) = [1 1 1; -1 3 2; -4 -1 2];

for i = 1 : 3
    fprintf("\n*************  Sub-part %d  *************\n\n", i);
    A = array(:, :, i);

    [iter, lambda] = Rayleigh(A, x, 15);
    [V, D] = eig(A, 'vector');
    [D, I] = sort(D, 'descend');
    V = V(:, I);
   
    [lambda_1, lambda_2] = deal(D(1), D(2));
    theoretical_rate = abs(lambda_2) / abs(lambda_1);
    
    [val, idx] = max(abs(V(:, 1)));
    v = V(:, 1);
    v = v / v(idx);
    experimental_rate = norm(iter(:, end) - v) / norm(iter(:, end - 1) - v);
   
    fprintf("Theoretical rate of convergence = ");
    disp(theoretical_rate);
   
    fprintf("Experimental rate of convergence = ");
    disp(experimental_rate);
    
%     fprintf("Theoretical eigen vector = \n");
%     disp(v);
%    
%     fprintf("Experimental eigen vector = \n");
%     disp(iter(:, end));
end