function Analytical
   format long
   L = 1;
   D = 0.01;
   u = 0.1;
   dx = 0.1;
   xval = 0:dx:L;
   P = (u*L)/D;
   N = 10;
   A = zeros(N+1);
syms k x;
        betak = ((P/2)^2)+(k*pi)^2;
        lambdak = (D*betak)/(L^2);
        Ak =@(x,t) symsum((((-1)^k)*(k/betak)*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);
        Bk =@(x,t) symsum(((((-1)^(k+1))*(k/betak)*(1+(P/betak))*exp((-P)/2)+((k*P)/((betak)^2)))*sin((k*pi*x)/L)*exp(-(lambdak*t))),k,1,50);

for t = 1:N+1    
    for x = 1:N+1
        A(x,t) = 100*(((exp(P*xval(x))-1)/(exp(P)-1))+((4*pi*exp((P*xval(x))/2)*sinh(P/2))/(exp(P)-1))*Ak(xval(x),t-1)+(2*pi*exp((P*xval(x))/(2*L)))*Bk(xval(x),t-1));
    end
end
   
B = A.';

plot(xval,B(1,:))
hold on
plot(xval,B(3,:),'-o')
hold on
plot(xval,B(6,:),'-x')
hold on
plot(xval,B(11,:),"-.")
hold off    
end