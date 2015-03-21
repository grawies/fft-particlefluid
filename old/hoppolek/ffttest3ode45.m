N = 100;
L = 2*pi;

x = 0:L/(N-1):L;

y1 = 0.0;
y2 = 0.5;

u0 = y2*0.6*exp(-4*(x-3).^2);

plot(x,u0);
hold all;

f = dydt(t,y) -0.5*ddx(y);

y0 = [10;5;4;3;1];
[t,y] = ode45(dydt,[0,5],y0);

surf(y);
