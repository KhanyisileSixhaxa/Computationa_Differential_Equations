function u = AnalyticalSol(y,K,L,t,u,alpha)
M = length(y);
B = zeros(M,1);
B(1,1) = 1;
B(M,1) = 0;
for i = 2:M-1
    a = y(i);
    sum = 0;
    for k = 1:K
%         sum = sum + (sin(((2*j-1)*pi*(a - u*t))/L)*exp(-alpha*((2*j-1)^2)*(pi^2)*t/(L^2)))/(2*j-1);
          sum = sum + sin((2*k-1)*((pi*(a-u*t))/L))*exp((-alpha*(2*k-1)^(2)*pi^(2)*(t/(L^2))))/(2*k-1);
    end
    B(i,1) = 0.5 - (2/pi)*sum;
end
u = B;
end