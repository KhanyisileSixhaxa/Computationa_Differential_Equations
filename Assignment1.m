 function Assignment1
     
    %p = [10,7,11,15,9,7,8,9,4,7,14,7,9,12,10];
    %w = [10,5,7,8,14,9,9,9,12,15,8,8,10,8,10];
    N = 20;
    n = 15;
    m = 6;
    %R = 100;
    Imax = 25;
    A = zeros(20,16);
    x = zeros(1,15);
    y = zeros(1,15);
    l = zeros(1,3);
    C = zeros(16,16);
    
    function [fcall] = Fcall(x)
        R = 100;
        p = [10,7,11,15,9,7,8,9,4,7,14,7,9,12,10];
        w = [10,5,7,8,14,9,9,9,12,15,8,8,10,8,10];
        sum = 0;
        sum1 = 0;
        n = 15;
        
        for s = 1:n
            sum = sum + (w(1,s)*x(1,s));
            sum1 = sum1 + (p(1,s)*x(1,s));
        end
        
        sum = sum - 100;
        
        if sum <= 0
            s2=0;
        else
            s2 = sum;
        end
        fcall = sum1 - (R*s2);
end
    
    for i = 1:N
        for j = 1:n
            if rand >=0.5
                x(j)=0;
            else
                x(j)=1;
            end
         A(i,j) = x(j);
        end
        A(i,n+1) = Fcall(x);
    end
    A = sortrows(A,n+1);
    
    for I = 1:Imax
        for k = 1:2:m-1
          for j = 1:2
            n1 = 1;
            n2 = 1;
            
            while n1==n2
                n1 = 1 + floor(rand*N);
                n2 = 1 + floor(rand*N);
            end
            
            if A(n1,n+1) > A(n2,n+1)
                l(j) = n1;
            else
                l(j) = n2;
            end
          end
          
          for i=1:n
              C(k,i) = A(l(1),i);
              C(k+1,i) = A(l(2),i);
          end
          
          n1 = 1;
          n2 = 2;
          
          while n1==n2
             n1 = 1 + floor(rand*n);
             n2 = 1 + floor(rand*n);
          end
          
          if n1>n2
              for i=n2:n1
                  C(k,i) = A(l(2),i);
                  C(k+1,i) = A(l(1),i);
              end
          else
              for i=n1:n2
                  C(k,i) = A(l(2),i);
                  C(k+1,i) = A(l(1),i);
              end
          end
          
          for j=1:n
              x(j) = C(k,j);
              y(j) = C(k+1,j);
          end
          
          C(k,n+1) = Fcall(x);
          C(k+1,n+1) = Fcall(y);
          
        end
        
        for i=1:m
            for j=1:n+1
                A(i,j) = C(i,j);
            end
        end
        A = sortrows(A,n+1);
    end
    
    xstar = A(N,1:n)
    fmin = A(N,n+1)
    
end