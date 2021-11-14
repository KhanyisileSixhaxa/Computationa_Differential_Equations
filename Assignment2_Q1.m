%Modified Code for plotting the Absolute Error at diffirent changes of x%
function [xval,AbsE] = Assignment2_Q1(dx)
   %initializing all constants%  
   L = 1;
   D = 0.01;
   u = 0.1;
   M = 500;
   xval = 0:dx:1;%an array of all the x values 
   N = length(xval);
   alpha = (u*dt)/(dx); 
   beta = (D*dt)/(dx^2);
   %declaring my C matrix which is our Numerical Solution Matrix and imposing the boundary conditions%
   C = zeros(M+1,N); 
   C(:,1)=0;
   C(:,N)=100;
   C(1,:) = (100/L)*xval;
   %declaring the U(m+1) which is the RHS and U(m) matrix which is the LHS
   RHS = zeros(N);
   LHS = zeros(N);
   P = (u*L)/D;%Peclect Number
   A = zeros(N);%Exact Solution Matrix  
   %Initiliazing the coefficients for the exact solution
   syms k x;
        betak = ((P/2)^2)+(k*pi)^2;
        lambdak = (D*betak)/(L^2);
        Ak =@(x,t) symsum((((-1)^k)*(k/betak)*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);
        Bk =@(x,t) symsum(((((-1)^(k+1))*(k/betak)*(1+(P/betak))*exp((-P)/2)+((k*P)/((betak)^2)))*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);
%Generating the exact solution matrix
for t = 1:N    
    for x = 1:N
        A(x,t) = 100*(((exp(P*xval(x))-1)/(exp(P)-1))+ ...
            ((4*pi*exp((P*xval(x))/2)*sinh(P/2))/(exp(P)-1))*Ak(xval(x),t-1)+(2*pi*exp((P*xval(x))/(2*L)))*Bk(xval(x),t-1));
    end
end
    
 TA = A.';%Transposing the exact solution
%Filling in the RHS and LHS matrix with constants 
for i=2:N-1
    RHS(i,i-1)=-(alpha+2*beta);
    RHS(i,i+1)=alpha-2*beta;
    RHS(i,i)=4+4*beta;
    LHS(i,i-1)=alpha+2*beta;
    LHS(i,i+1)=-(alpha-2*beta);    
    LHS(i,i)=4-4*beta;
end

RHS(1,1) = 1;
RHS(N, N) = 1;
LHS(1,1) = 1;
LHS(N, N) = 1; 
    %Generating the Numerical Solution matrix
   for i=2:M+1
       b = C(i-1,:);
       C(i,:) = inv(RHS)*LHS*(b.');  
       C(i,1) = 0;
       C(i,N) = 100;
   end
    
   %Calculating the absolute error 
   AbsE = abs(C(501,:)-TA(6,:));

end

