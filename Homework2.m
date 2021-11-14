function h2
    
    N=20;
    M=25;
    dx=0.05;
    dt=0.013;
    U = zeros(N+1,M+1);
    A = zeros(N+1,N+1);
    c = 1;
    lambda = (c*dt)/(dx^2);
    xvalues=0.05:dx:0.95;
%     disp(xvalues)
    
    
    syms f(x)
        f(x) = piecewise(0<=x<=0.5,2*x,0<=x<=1,2-2*x);
        U(2:N,1) = f(xvalues);
    
    for m=2:M
        for n=2:N
            if(n==m)
                A(m,n) = 1-2*lambda;
                A(m,n+1) = lambda;
                A(m,n-1) = lambda;
            end
        end
    end
    
%     disp(U(1:N+1,1))
    
    for m=2:M
        U(1:N+1,m) = A*U(1:N+1,m-1);
    end
    
   disp(U)
end