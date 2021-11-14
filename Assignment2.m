function Assignment2
    
    f = @(x1,x2) 0.25*(x1^4) - 0.5*(x1^2) + 0.1*x1 + 0.5*x2;
    N = 30;
    m = 6;
    Imax = 45;
    n = 2;
    x = zeros(1,2);
    y = zeros(1,2);
    A = zeros(30,3);
    C = zeros(6,3);
    L = [-5,-5];
    U = [6,6];
    
    
    for i = 1:N
        for j = 1:n
            A(i,j) = L(j) + [(U(j)-L(j))*rand];
            x(j) = A(i,j);
        end
        A(i,n+1) = f(x(1),x(2));
    end
    A = sortrows(A,-(n+1));
    
    for I = 1:Imax
        for k = 1:2:m-1
            
            for j = 1:2
               n1 = 1;
               n2 = 1;
                    while n1==n2
                    n1 = 1 + floor(N*rand);
                    n2 = 1 + floor(N*rand);
                    end
                    
                    if A(n1,n+1) < A(n2,n+1)
                       L(j) = n1;
                    else
                       L(j) = n2;
                    end
            end
            
            for j = 1:n
            
                alpha1 = -0.5 + 2*rand;
                x(j) = alpha1*A(L(1),j) + (1 - alpha1)*A(L(2),j);
                
                if x(j) <= L(j) || x(j) >= U(j)
                    x(j) = L(j) + rand*(U(j) - L(j));
                end
                
                alpha2 = -0.5 + 2*rand;
                y(j) = alpha2*A(L(2),j) + (1 - alpha2)*A(L(1), j);
                
                if y(j) <= L(j) || y(j) >= U(j)
                
                    y(j) = L(j) + rand*(U(j) - L(j));
                    
                end
                
            C(k,j) = x(j);
            C(k+1,j) = y(j);
            
            end
        C(k,n+1) = f(x(1),x(2));
        C(k+1,n+1) = f(y(1),y(2));
        
        end
        
        for i = 1:m
         for j = 1:n+1
             A(i,j) = C(i,j);
         end
        end
    A = sortrows(A,-(n+1));
    fmin = A(N,n+1)
    xstar = A(N,1:n)
end