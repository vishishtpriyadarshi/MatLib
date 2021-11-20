function B = FrancisQRS(A)
    [n, ~]= size(A);
    
    for i = 1 : n - 1
        rho = A(end, end);
        idx = 1;
        if i > 1
           idx = i - 1;
        end
        
        p = A(i : i + 1, idx);
        if(i == 1)
            p(1) = p(1) - rho;
        end    
        
        [u, gamma, ~] = reflect(p);
        A(i : i + 1, idx : n) = appreflect(u, gamma, A(i : i + 1, idx : n));
        A(idx : n, i : i + 1) = appreflect_modified(u, gamma, A(idx : n, i : i + 1));       
    end
    
    B = A;
end