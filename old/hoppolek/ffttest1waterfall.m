N = 100;
L = 2*pi;

x = 1/N:L/N:L;

y1 = 0.0;
y2 = 0.5;

%u0 = y2*0.6*exp(-4*(x-3).^2);
%u0 = y2*0.6*(abs(sin(x)-0.5)<0.2);
u0 = y2*0.3*exp(-100*(x-1).^2);
c = ffts(u0);
k = -N/2:1:N/2-1;
dt = 1.9/N;

k2 = k.*k;
nu = 0.06;
% u_t = u_xx
%dcdt = -k2.*c;
% u_t = -u*u_x + nu*u_xx    Burger's equation
%dcdt = @(t, c) -ffts(iffts(c).*iffts(1j*k.*c)) - nu*k2.*c;
% variable coefficient wave eq u_t + F(x)*u_x = 0
velprop = 0.2 + sin(x-1).^2;
dcdt = @(t, c) -ffts(velprop.*iffts(1j*k.*c));
Nplots = 50;
data = [u0; zeros(Nplots,N)];
Nt = 500;
plotInterval = Nt/Nplots;
t = dt*(1:Nt);
ploti=1;
for i = 1:Nt
    tc = t(i);
    k1 = dcdt(tc,c);
    k2 = dcdt(tc+.5*dt, c+.5*k1*dt);
    k3 = dcdt(tc+.5*dt, c+.5*k2*dt);
    k4 = dcdt(tc+dt, c+k2*dt);
    c = c + dt * (k1 + 2*k2 + 2*k3 + k4) /6;
    if mod(i,plotInterval) == 0
        ploti = ploti + 1;
        data(ploti,:) = real(iffts(c));
    end
end
plotxi = 1:N;
waterfall(x(plotxi),dt*plotInterval*(1:Nplots+1),data(:,plotxi));
axis([0 L 0 Nt*dt 0 y2]);