N = 100;
L = 2*pi;

x = 1/N:L/N:L;

y1 = 0.0;
y2 = 0.5;

ffts = @(u) fftshift(fft(u));
iffts = @(u) ifft(ifftshift(u));

u0 = y2*0.6*exp(-4*(x-3).^2);
%u0 = y2*0.6*(abs(sin(x)-0.5)<0.2);

plot(x,u0);
hold all;

c = ffts(u0);
u = u0;
k = -N/2:1:N/2-1;
dt = 0.02;

k2 = k.*k;
nu = 0.1;
% u_t = u_xx
%dcdt = -k2.*c;
% u_t = -u*u_x + nu*u_xx    Burger's equation
%dcdt = @(t, c) -ffts(iffts(c).*iffts(1j*k.*c)) - nu*k2.*c;
dudt = @(t, u) -u.*iffts(1j*k.*ffts(u)) - nu*iffts(k2.*ffts(u));
% variable coefficient wave eq u_t + F(x)*u_x = 0
%dcdt = @(t, c) -ffts(iffts(c).*iffts(1j*k.*c)) - nu*k2.*c;

for t = 0:dt:500*dt
    k1 = dudt(t,u);
    k2 = dudt(t+.5*dt, u+.5*k1*dt);
    k3 = dudt(t+.5*dt, u+.5*k2*dt);
    k4 = dudt(t+dt, u+k2*dt);
    u = u + dt * real(k1 + 2*k2 + 2*k3 + k4) /6;
    
    %c = c + dt * dcdt;
    clf;
    plot(x,u,'b');
    axis([0 L y1 y2]);
    drawnow();
end
