
% constants
N = 64; dt = .0001;
tmax = dt*1000;
L = 2*pi; nu = .04;

% space vectors
x = (0:N-1)*L/N;
[x1,x2,x3] = meshgrid(x,x,x);

% k-space vectors
k = -N/2:N/2-1;
[k1,k2,k3] = meshgrid(k,k,k);
ksq = k1.^2+k2.^2+k3.^2;
ksqinv = 1./ksq;

% initial values u0, v0, w0
u = exp(-4*((x1-L/2).^2 + (x2-L/2).^2));
v = -u;
w = 0*x1;
% ...and U0, V0, W0
U = fftns(u);
V = fftns(v);
W = fftns(w);

% timestep data
nmax = round(tmax/dt); t = 0;
% plot specs
plotnum = 100; nplt = nmax/plotnum;
I = 1:N/16:N;

% save data
%Udata = zeros(N,N,N,nmax+1);
%Udata(:,:,:,1) = U;
%Vdata = zeros(N,N,N,nmax+1);
%Vdata(:,:,:,1) = V;
%Wdata = zeros(N,N,N,nmax+1);
%Wdata(:,:,:,1) = W;

% derivatives for Euler steppin'
dUdt = @(U,V,W) 0*U;
dVdt = @(U,V,W) 0*V;
dWdt = @(U,V,W) 0*W;

% plotting preparations
figure(1)
plotN = 20; % plot every plotN:th frame

% timestep
for n =1:nmax
    t = t + dt;
    E = exp(nu * ksq * t);
    
    U = U + dt * dUdt(U,V,W) .* E;
    V = V + dt * dVdt(U,V,W) .* E;
    W = W + dt * dWdt(U,V,W) .* E;
    
    %Udata(:,:,:,n+1) = U;
    %Vdata(:,:,:,n+1) = V;
    %Wdata(:,:,:,n+1) = W;
    
    % plot
    if mod(n,plotN) == 0
        u = ifftns(U./E);
        v = ifftns(V./E);
        clf;
        quiver(x1(I,I,N/2),x2(I,I,N/2),u(I,I,N/2),v(I,I,N/2),0);
        axis([0, L, 0, L]);
        drawnow;
    end
end

