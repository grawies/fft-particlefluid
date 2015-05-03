
L=2*pi; g = @(x,y,z) sin(L*cos(x*2*pi/L)).*cos(y*2*pi/L) + sin(z*2*pi/L);
[ns,re] = ShowConvergence(8,24,2*pi, g);
L=2.0; g = @(x,y,z) sin(x).*cos(y).*(2/L^2 * z.^3 - 3/L * z.^2 + z);
[nsn,ren] = ShowConvergence(8,24,2*pi, g);
figure(2)
loglog(ns,re,'r-')
hold on
loglog(nsn,ren,'b--')
legend('smooth','non-smooth');

loglog(ns,re,'ko','MarkerSize',12);
loglog(nsn,ren,'kx','MarkerSize',12);

grid on;
xlabel('Grid size N','FontSize',12);
ylabel('relative 2-norm error in velocity','FontSize',12);
title('convergence of solution for two forcing functions','FontSize',12);