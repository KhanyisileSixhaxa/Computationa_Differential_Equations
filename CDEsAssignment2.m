function CDEsAssignment2
   %initializing all constants% 
   L = 1;
   D = 0.01;
   u = 0.1;
   dx = 0.1;
   dt = 0.01;
   N = 10;
   M = 500;
   alpha = (u*dt)/(dx);
   beta = (D*dt)/(dx^2);
   %an array of all the x values 
   xval = 0:dx:1; 
   %declaring my C matrix which is our Numerical Solution Matrix and imposing the boundary conditions%
   C = zeros(M+1,N+1);
   C(:,1)=0;
   C(:,N+1)=100;
   C(1,:) = (100/L)*xval;
   %declaring the U(m+1) which is the RHS and U(m) matrix which is the LHS
   RHS = zeros(N+1);
   LHS = zeros(N+1);
   P = (u*L)/D;%Peclect Number
   A = zeros(N+1); %Exact Solution Matrix
   
   %Initiliazing the coefficients for the exact solution
   syms k x;
        betak = ((P/2)^2)+(k*pi)^2;
        lambdak = (D*betak)/(L^2);
        Ak =@(x,t) symsum((((-1)^k)*(k/betak)*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);
        Bk =@(x,t) symsum(((((-1)^(k+1))*(k/betak)*(1+(P/betak))*exp((-P)/2)+((k*P)/((betak)^2)))*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);
%Generating the exact solution matrix
for t = 1:N+1    
    for x = 1:N+1
        A(x,t) = 100*(((exp(P*xval(x))-1)/(exp(P)-1))+ ...
            ((4*pi*exp((P*xval(x))/2)*sinh(P/2))/(exp(P)-1))*Ak(xval(x),t-1)+(2*pi*exp((P*xval(x))/(2*L)))*Bk(xval(x),t-1));
    end
end
    
 TA = A.';%Transposing the exact solution
 
 
%Filling in the RHS and LHS matrix with constants
for i=2:N
    RHS(i,i-1)=-(alpha+2*beta);
    RHS(i,i+1)=alpha-2*beta;
    RHS(i,i)=4+4*beta;
    LHS(i,i-1)=alpha+2*beta;
    LHS(i,i+1)=-(alpha-2*beta);    
    LHS(i,i)=4-4*beta;
end

RHS(1,1) = 1;
RHS(N+1, N+1) = 1;
LHS(1,1) = 1;
LHS(N+1, N+1) = 1; 
   %Generating the Numerical Solution matrix
   for i=2:M+1
       b = C(i-1,:);
       C(i,:) = inv(RHS)*LHS*(b.');  
       C(i,1) = 0;
       C(i,N+1) = 100;
   end
   
   %Plotting the Numerical and Analytical solution
   figure
   plot(xval,C(1,:),'r')
   hold on
   plot(xval,TA(1,:),'b*')
   plot(xval,C(101,:),'r')
   plot(xval,TA(2,:),'b*')
   plot(xval,C(201,:),'r')
   plot(xval,TA(3,:),'b*')
   plot(xval,C(501,:),'r')
   plot(xval,TA(6,:),'b*')
   legend('Numerical ','Analytical');
   xlabel('x');
   ylabel('C(x,t)')
   hold off
   


end