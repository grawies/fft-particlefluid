function [a,b,c] = stokes3dconvergence1fcn(N,N0,deltas)

%% PRELIMINARIES
% constants
M = 1000; % # of particles
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

% particle positions
rng(239847284);
% particle positions
R0 = 1.0;
phi = unifrnd(0,2*pi, 1, M);
costheta = unifrnd(-1,1, 1, M);
u = unifrnd(0,1,1,M);
theta = acos( costheta ) - pi/2;
r = R0 * u.^(1/3);

particlex1 = .5 * L + r .* sin( theta ) .* cos( phi );
particlex2 = .5 * L + r .* sin( theta ) .* sin( phi );
particlex3 = .5 * L + r .* cos( theta );

%% FUNTIONS, EQUATION PROPERTIES

% gravity force on point particles
g = 9.82;
mrel = 0.1;
frel = mrel * g;
fgx = frel * 0;
fgy = frel * 0;
fgz = frel * -1;
% delta function half-width
epsilon = L/N0; % = h

%% DO (...WHILE)

hdeltasum = 0;

for i = 1:M
    hdelta1 = max(1-abs(x1-particlex1(i))/epsilon,0);
    hdelta2 = max(1-abs(x2-particlex2(i))/epsilon,0);
    hdelta3 = max(1-abs(x3-particlex3(i))/epsilon,0);
    hdelta = hdelta1.*hdelta2.*hdelta3;
    hdeltasum = hdeltasum + hdelta;

end

fx = fgx * hdeltasum;
fy = fgy * hdeltasum;
fz = fgz * hdeltasum;

% calculate stokes flow from the particle force

fhx = fftns(fx);
fhy = fftns(fy);
fhz = fftns(fz);
figure();
%semilogy(sqrt(ksq(:)), abs(fhx(:)),'*');
kdotfh = (k1.*fhx + k2.*fhy + k3.*fhz).*C2;

u = ifftns(fhx.*C1 - k1.*kdotfh);
v = ifftns(fhy.*C1 - k2.*kdotfh);
w = ifftns(fhz.*C1 - k3.*kdotfh);

if deltas == 1
    a = fhx;
    b = fhy;
    c = fhz;
else
    a = u;
    b = v;
    c = w;
end
