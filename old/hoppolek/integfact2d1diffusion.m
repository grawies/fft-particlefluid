
N = 64;
t = 0; dt = 0.08;
x = (0:N-1)*2*pi/(N-1);
y = x;
[x,y] = meshgrid(x,y);

nu = 0.02;
k = -N/2:N/2-1;
l = k;
[k,l] = meshgrid(k,l);
k2 = k.*k;
l2 = k2;
k2l2 = k2+l2;

u = exp(-((x-pi).^2+(y-pi).^2));
U = fft2s(u);

% g = forcing func / rho
g = @(t,x,y) 0;
ghat = @(t,k,l) 0;

tmax = 8.0;
nmax = round(tmax/dt);
plotnum = 100;
nplt = nmax/plotnum;

for n=1:nmax
    t = t+dt;
    E = exp(nu * (k2+l2) * t);
    U = U;% + dt * (ghat(t,k) - .5i*k.*ffts(u.*u)).*E;
    u = ifft2s(U ./ E);
    if mod(n, nplt) == 0
        surf(x,y,real(u));
        axis([0 2*pi 0 2*pi 0 1]);
        drawnow;
        pause(0.1);
        %tdata = [tdata; t];
    end
end
%udata = real(udata);
%waterfall(x,tdata,udata), colormap(1e-6*[1 1 1]);
%view(-163,10)
%xlabel x, ylabel t, axis([0 2*pi 0 tmax 0 1]), grid off

