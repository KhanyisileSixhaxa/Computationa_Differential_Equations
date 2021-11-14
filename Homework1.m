%1590202
%Khanyisile Sixhaxa
%CDEs
function CDEHomework2
    
    N = 20;
    M = 50;
    dx = 0.05;
    c = 1;
    dt = 0.0012;
    lambda = (c*dt)/dx^2;
    U = zeros(M+1,N+1);
    xval = 0:dx:1;
    U(:,1)=0;
    U(:,N+1)=0;
    Ut1 = zeros(1,N+1);
    Ut2 = zeros(1,N+1);
    
    syms f(x)
        f(x) = piecewise(0<=x<=0.5,2*x,0.5<=x<=1,2-2*x);
        
    U(1,:) = f(xval);
    
    for m=1:M
        for n=2:N
            U(m+1,n) = lambda*U(m,n-1) + (1-(2*lambda))*U(m,n) + lambda*U(m,n+1);
        end
    end
    Ut0 = U(1,:);
    Ut1 = U(26,:);
    Ut2 = U(51,:);
    
    ExactUt0 = sin(pi*xval)*exp(-c*(pi^2)*0);
    ExactUt1 = sin(pi*xval)*exp(-c*(pi^2)*0.03);
    ExactUt2 = sin(pi*xval)*exp(-c*(pi^2)*0.06);
    
    plot(xval,Ut0);
    hold on
    plot(xval,ExactUt0);
    legend({'Numerical Approximation','Exact Approximation'});
    xlabel('x');
    ylabel('U(x,t)');
    title(['Comparison of the numerical and exact solution at t=0 and dt=0.0012']);
    hold off
    
    
 end