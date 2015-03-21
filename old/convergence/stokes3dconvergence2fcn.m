function [u,v,w] = stokes3dconvergence1fcn(N)

%% PRELIMINARIES
% constants
L = 2*pi;
nu = 0.4;

% space vectors
x = (0:N-1)*L/N;
[x1,x2,x3] = ndgrid(x,x,x);

% k-space vectors
k = -N /2:N/2-1;
[k1,k2,k3] = ndgrid(k,k,k);
ksq = k1.^2+k2.^2+k3.^2;
ksqinv = 1./ksq;
ksqinvsq = 1./(ksq.*ksq);
ksqinv(N/2+1,N/2+1,N/2+1) = 0;
ksqinvsq(N/2+1,N/2+1,N/2+1) = 0;

C1 = ksqinv / nu;
C2 = ksqinvsq / nu;

%% FUNCTION CALCULATIONS

fx = sin(2*x1).*sin(x2).*sin(3*x3);
fy = sin(3*x1).*sin(2*x2).*sin(x3);
fz = sin(x1).*sin(3*x2).*sin(2*x3);

% calculate stokes flow from the particle force

fhx = fftns(fx);
fhy = fftns(fy);
fhz = fftns(fz);

kdotfh = (k1.*fhx + k2.*fhy + k3.*fhz).*C2;

u = ifftns(fhx.*C1 - k1.*kdotfh);
v = ifftns(fhy.*C1 - k2.*kdotfh);
w = ifftns(fhz.*C1 - k3.*kdotfh);

