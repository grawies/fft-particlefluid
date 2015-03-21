
N = 64;
t = 0; dt = 0.0008;
x = (0:N-1)*2*pi/(N-1);
y = x;
[x,y] = meshgrid(x,y);

nu = 0.2;
k = -N/2:N/2-1;
l = k;
[k,l] = meshgrid(k,l);
k2 = k.*k;
l2 = k2;
k2l2 = k2+l2;

u = 0*x;%exp(-((x-pi).^2+(y-pi).^2));
v = 0*x;%exp(-4*((x-0.5*pi).^2+(y-1.5*pi).^2));
U = fft2s(u);
V = fft2s(v);

w = 5;
% g = forcing func / rho
gx = @(t,x,y) -0.1*sin(w*t)*exp(-4*((x-cos(w*t)-pi).^2 + (y-sin(w*t)-pi).^2));
gxhat = @(t,k,l) fft2s(gx(t,x,y));
gy = @(t,x,y) 0.1*cos(w*t)*exp(-4*((x-cos(w*t)-pi).^2 + (y-sin(w*t)-pi).^2));
gyhat = @(t,k,l) fft2s(gy(t,x,y));

tmax = 3.0;
nmax = round(tmax/dt);
plotnum = 100;
nplt = nmax/plotnum;

% GIF INITIALIZATION
figure(1)
filename = 'burger2d_swirl.gif';
gifdelay = 0.01;
surf(x,y,abs(u)+abs(v), 'EdgeColor', 'none');
axis([0 2*pi 0 2*pi 0 1]);
view(0,90)
drawnow;
frame=getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',gifdelay);

I = 1:N/16:N;

dudt = @(t,u,v) (...
        ...%forcing fcn
        + gxhat(t,k,l)...
        ...% -.5 * (u^2)_x
        - .5i*k.*fft2s(u.*u)...
        ...% - v * u_y
        - fft2s(v.*ifft2s(1i*l.*fft2s(u)))...
        ...% integrating factor
        );
dvdt = @(t,u,v) (...
        + gyhat(t,k,l)...
        - .5i*l.*fft2s(v.*v)...
        - fft2s(u.*ifft2s(1i*k.*fft2s(v)))...
        );

for n=1:nmax
    t = t+dt;
    E = exp(nu * (k2+l2) * t);
    
    U = U + dt * dudt(t,u,v).*E;
    V = V + dt * dvdt(t,u,v).*E;
    u = ifft2s(U ./ E);
    v = ifft2s(V ./ E);
    
    if mod(n, nplt) == 0
        clf;
        surf(x,y,sqrt(abs(u.^2)+abs(v.^2)), 'EdgeColor', 'none');
        axis([0 2*pi 0 2*pi 0 1]);
        view(0,90)
        %view(-20,45)
        %clf;
        %quiver(x(I,I),y(I,I),u(I,I),v(I,I));
        %axis([0 2*pi 0 2*pi]);
        hold on;
        plot([pi+cos(w*t)],[pi+sin(w*t)],'r*');
        drawnow;
        
        % GIF
        frame=getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',gifdelay);

        %pause(0.1);
        %tdata = [tdata; t];
    end
end
%udata = real(udata);
%waterfall(x,tdata,udata), colormap(1e-6*[1 1 1]);
%view(-163,10)
%xlabel x, ylabel t, axis([0 2*pi 0 tmax 0 1]), grid off

