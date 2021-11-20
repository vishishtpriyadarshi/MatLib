kappa = [10^5, 10^7, 10^11];
idx = 1;
data = [];

for k = kappa
    % mat = condmat(50, k);
    % save(sprintf('gram_schmidt-%d.mat', idx), 'mat');
    load(sprintf('gram_schmidt-%d.mat', idx));
    
    [Q1, R1] = cgs(mat);
    [Q2, R2] = mgs(mat);
    [Q3, R3] = cgsrep(mat);
    
    % Condensed QR using matlab command
    [Q4, R4] = qr(mat);
    Q4 = Q4(:, 1 : size(mat, 2));
    R4 = R4(1 : size(mat, 2), :);

    I = eye(50);
    
    res = [idx, k, norm(Q1' * Q1 - I), norm(Q2' * Q2 - I), norm(Q3' * Q3 - I), norm(Q4' * Q4 - I)];
    data = [data; res];
    idx = idx + 1;
end

fprintf("***********  Table representing departure from orthogonality (i.e., ||Q'Q - I||)  ***********\n\n")
T = array2table(data);
T.Properties.VariableNames = {'SI No.', 'Condition Number', 'cgs', 'mgs', 'cgsrep', 'qr'};
disp(T);