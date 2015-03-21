
N = 128;
t = 0; dt = 0.005;
x = (0:N-1)*2*pi/(N-1);
y = x;
[x,y] = meshgrid(x,y);

nu = 0.04;
k = -N/2:N/2-1;
l = k;
[k,l] = meshgrid(k,l);
k2 = k.*k;
l2 = k2;
k2l2 = k2+l2;

u = exp(-((x-pi).^2+(y-pi).^2));
v = exp(-4*((x-0.5*pi).^2+(y-1.5*pi).^2));
U = fft2s(u);
V = fft2s(v);

% g = forcing func / rho
gx = @(t,x,y) 0;%-exp(-4*(x-pi).^2);
gxhat = @(t,k,l) 0;%fft2s(gx(t,x,y));
gy = @(t,x,y) 0;
gyhat = @(t,k,l) 0;

tmax = 2.0;
nmax = round(tmax/dt);
plotnum = 50;
nplt = nmax/plotnum;

% GIF INITIALIZATION
figure(1)
filename = 'burger2d.gif';
gifdelay = 0.01;
surf(x,y,abs(u)+abs(v), 'EdgeColor', 'none');
axis([0 2*pi 0 2*pi 0 1]);
view(0,90)
drawnow;
frame=getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',gifdelay);


for n=1:nmax
    t = t+dt;
    E = exp(nu * (k2+l2) * t);
    
    U = U + dt * (...
        ...%forcing fcn
        + gxhat(t,k,l)...
        ...% -.5 * (u^2)_x
        - .5i*k.*fft2s(u.*u)...
        ...% - v * u_y
        - fft2s(v.*ifft2s(1i*l.*fft2s(u)))...
        ...% integrating factor
        ).*E;
    
    V = V + dt * (...
        + gyhat(t,k,l)...
        - .5i*l.*fft2s(v.*v)...
        - fft2s(u.*ifft2s(1i*k.*fft2s(v)))...
        ).*E;
    
    u = ifft2s(U ./ E);
    v = ifft2s(V ./ E);
    if mod(n, nplt) == 0
        surf(x,y,abs(u)+abs(v), 'EdgeColor', 'none');
        axis([0 2*pi 0 2*pi 0 1]);
        view(0,90)
        drawnow;
        
        % GIF
        frame=getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',gifdelay);

        pause(0.1);
        %tdata = [tdata; t];
    end
end
%udata = real(udata);
%waterfall(x,tdata,udata), colormap(1e-6*[1 1 1]);
%view(-163,10)
%xlabel x, ylabel t, axis([0 2*pi 0 tmax 0 1]), grid off

