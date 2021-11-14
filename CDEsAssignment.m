function CDEsAssignment
    
   L = 1;
   D = 0.01;
   u = 0.1;
   dx = 0.1;
   dt = 0.01;
   xval = 0:dx:L;
   N = length(xval)
   M = 1000;
   beta = (D*dt)/(dx^2);
   alpha = (u*dt)/(2*dx);
   C = zeros(M,N);
   
   C(:,1)=0;
   C(:,N)=100;
   
%    C0=zeros(1,N+1);
%    C1=zeros(1,N+1);
%    C2=zeros(1,N+1);
%    C3=zeros(1,N+1);
%    
%    Abs0=zeros(1,N+1);
%    Abs1=zeros(1,N+1);
%    Abs2=zeros(1,N+1);
%    Abs3=zeros(1,N+1);
%    
  length(C(1,:))
   C(1,:) = (100/L)*xval;
   P = (u*L)/D;
   A = zeros(N);

   
   for m=1:M-1
        for n=2:N-1
            C(m+1,n) = (1-2*beta)*C(m,n)+(beta-alpha)*C(m,n+1)+(alpha+beta)*C(m,n-1);
        end
   end
   
   C0 = C(1,:);
   C1 = C((2/dt)+1,:);
   C2 = C((5/dt)+1,:);
   C3 = C(M,:);
   

syms k x;
        betak = ((P/2)^2)+(k*pi)^2;
        lambdak = (D*betak)/(L^2);
        Ak =@(x,t) symsum((((-1)^k)*(k/betak)*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);
        Bk =@(x,t) symsum(((((-1)^(k+1))*(k/betak)*(1+(P/betak))*exp((-P)/2)+((k*P)/((betak)^2)))*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);

for t = 1:N    
    for x = 1:N
        A(x,t) = 100*(((exp(P*xval(x))-1)/(exp(P)-1))+ ...
            ((4*pi*exp((P*xval(x))/2)*sinh(P/2))/(exp(P)-1))*Ak(xval(x),t-1)+(2*pi*exp((P*xval(x))/(2*L)))*Bk(xval(x),t-1));
    end
end
    
 B = A.';
  Abs0 = abs(B(1,:)-C0);
  Abs1 = abs(B(3,:)-C1);
  Abs2 = abs(B(6,:)-C2);
  Abs3 = abs(B(11,:)-C3);
  
  figure
  plot(xval,C0,'-o')
  hold on
  plot(xval,B(1,:),'-*');
  hold on
  plot(xval,C1,'-o')
  hold on
  plot(xval,B(3,:),'-*');
  hold on
  plot(xval,C2,'-o')
  hold on
  plot(xval,B(6,:),'-*');
  hold on
  plot(xval,C3,'-o')
  hold on
  plot(xval,B(11,:),'-*');
  hold on
  legend('Numerical Approximation at t=0','Analytical Approximation at t=0','Numerical Approximation at t=2','Analytical Approximation at t=2','Numerical Approximation at t=5','Analytical Approximation at t=5', ...
    'Numerical Approximation at t=10','Analytical Approximation at t=10'  );
  xlabel('x');
  ylabel('C(x,t)');
  title('Compatison of the numerical and analytical approximation at t=0, t=2, t=5, t=10');
  hold off

  
  figure
  plot(xval,Abs0,'-s');
  hold on
  plot(xval,Abs1,'-o');
  hold on
  plot(xval,Abs2,'-x');
  hold on
  plot(xval,Abs3,'-*');
  legend('t=0','t=2','t=5','t=7');
  xlabel('x');
  ylabel('Absolut Error');
  title(['Absolute Error of the numerical and analytical solution at t=0,t=2,t=5,t=10']);
  hold off
    
end