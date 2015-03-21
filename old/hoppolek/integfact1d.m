
N = 64;
t = 0; dt = 0.02;
x = (0:N-1)*2*pi/(N-1);

nu = 0.02;
k = -N/2:N/2-1;
k2 = k.*k;

u = exp(-4*(x-2).^2);
v = ffts(u);
udata = u;
tdata = 0;

% g = forcing func / rho
g = @(t,x) 0;
ghat = @(t,k) 0;

tmax = 4.0;
nmax = round(tmax/dt);
plotnum = 100;
nplt = nmax/plotnum;

% image properties
figure(1)
filename = 'burger.gif';
gifdelay = 0.05;
plot(x,u);
hold on;
drawnow;
frame=getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',gifdelay);

for n=1:nmax
    t = t+dt;
    E = exp(nu * k2 * t);
    v = v + dt * (ghat(t,k) - .5i*k.*ffts(u.*u)).*E;
    u = iffts(v ./ E);
    if mod(n, nplt) == 0
        udata = [udata; u];
        tdata = [tdata; t];
        clf;
        plot(x,real(u));
        axis([0 2*pi 0 1]);
        frame=getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',gifdelay);
    end
end
udata = real(udata);
waterfall(x,tdata,udata), colormap(1e-6*[1 1 1]);
view(-163,10)
xlabel x, ylabel t, axis([0 2*pi 0 tmax 0 1]), grid off

