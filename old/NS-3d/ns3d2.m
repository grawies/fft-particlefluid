%% PRELIMINARIES
% constants
N = 32; dt = .0002;
tmax = dt*1000;
L = 2*pi;
nu = .4;

% space vectors
x = (0:N-1)*L/N;
[x1,x2,x3] = meshgrid(x,x,x);

% k-space vectors
k = -N /2:N/2-1;
[k1,k2,k3] = meshgrid(k,k,k);
ksq = k1.^2+k2.^2+k3.^2;
ksqinv = 1./ksq;
ksqinv(N/2+1,N/2+1,N/2+1) = 0;

%% INITIAL VALUES AND EQUATION PROPERTIES

% initial values u0, v0, w0
e = exp(-2*((x1-L/2).^2 + (x2-L/2).^2));
u = -(x2-L/2) .* e;
v =  (x1-L/2) .* e;
w = 0*x1;
% ...and U0, V0, W0
U = fftns(u);
V = fftns(v);
W = fftns(w);

% forcing function
grav = 9.82; % actually g / rho  (rho = 1)
fox = @(t) 0 * x1;
foy = @(t) -grav * x2;
foz = @(t) 0 * x3;

% nonlinear virtual force (burger's term)
gx = 0*x1; % = - u.*grad u = 0 for now
gy = 0*x2;
gz = 0*x3;

% f + g
fpgx = @(t) fox(t) + gx;
fpgy = @(t) foy(t) + gy;
fpgz = @(t) foz(t) + gz;

% pressure term (-grad p / rho)
pterm = @(kdir, fx, fy, fz) - kdir .* (fx.*k1 + fy.*k2 + fz.*k3) .* ksqinv;

% derivatives for Euler steppin'
dUdt = @(t,U,V,W,fx,fy,fz) fx + pterm(k1,fx,fy,fz);
dVdt = @(t,U,V,W,fx,fy,fz) fy + pterm(k2,fx,fy,fz);
dWdt = @(t,U,V,W,fx,fy,fz) fz + pterm(k3,fx,fy,fz);

%% SIMULATION TIME STEPPING
% timestep constants
nmax = round(tmax/dt); t = 0;

% plot specs
plotnum = 500; nplt = nmax/plotnum;
I = 1:N/16:N;
figure(1);
% streamlines
[startx,startz] = meshgrid(x,x);
starty = L*ones(N,N);

% timestep
for n =1:nmax
    t = t + dt;
    E = exp(nu * ksq * t);
    
    fx = fftns(fpgx(t));
    fy = fftns(fpgy(t));
    fz = fftns(fpgz(t));
    
    U = U + dt * dUdt(t,U,V,W,fx,fy,fz) .* E;
    V = V + dt * dVdt(t,U,V,W,fx,fy,fz) .* E;
    W = W + dt * dWdt(t,U,V,W,fx,fy,fz) .* E;
    
    % plot
    if mod(n,nplt) == 0
        u = ifftns(U./E);
        v = ifftns(V./E);
        w = ifftns(W./E);
        
        clf;
        % QUIVER
        %{
        quiver(x1(I,I,N/2),x2(I,I,N/2),u(I,I,N/2),v(I,I,N/2),0);
        axis([0, L, 0, L]);
        %}
        % STREAMLINES
        streamline(stream3(x1,x2,x3,u,v,w,startx(I,I),starty(I,I),startz(I,I),[0.1]));
        axis([0 L 0 L 0 L]);
        %view(0,90);
        view(-25,50);
        %pause(0.1);
        % might have to pause if:
        %UIJ_AreThereWindowShowsPending - timeout waiting for window to show up
        drawnow;
    end
end

