lambda = 0.2;
t=-1:0.01:1;
c=(1/(lambda*sqrt(2*pi)))*exp(-1*t.^2/(2*(lambda)^2));

plot(t,c);