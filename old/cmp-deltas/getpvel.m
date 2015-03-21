function [px,py,pz,pu,pv,pw,gridx,gridy,gridz,gu,gv,gw] = getpvel(param)

%% PRELIMINARIES

% constants
%N = 16; dt = 0.01;
%M = 1000; % # of particles
N = param.N;
M = param.M;
dt = param.dt;
seed = param.seed;
nsteps = param.nsteps;
L = param.L;
nu = param.nu;
%nsteps=200;
tmax = dt*nsteps;
%L = 2*pi;
%nu = 0.4;

% space vectors
x = (0:N-1)*L/N;
[gridx,gridy,gridz] = ndgrid(x,x,x);

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
R0 = param.radius;
phi = unifrnd(0,2*pi, 1, M);
costheta = unifrnd(-1,1, 1, M);
u = unifrnd(0,1,1,M);
theta = acos( costheta );
r = R0 * u.^(1/3);

pointx = .5 * L + r .* sin( theta ) .* cos( phi );
pointy = .5 * L + r .* sin( theta ) .* sin( phi );
pointz = .5 * L + r .* cos( theta );

pointu = zeros(1,M);
pointv = zeros(1,M);
pointw = zeros(1,M);

%% FUNTIONS, EQUATION PROPERTIES

% gravity force on point particles
g = 9.82;
mrel = 0.1;
frel = mrel * g;
fgx = frel * 0;
fgy = frel * 0;
fgz = frel * -1;
% delta function half-width
epsilon = param.epsilon; % = h for some N

%% DO (...WHILE)

hdeltasum = deltasmear_tri(gridx,gridy,gridz,pointx,pointy,pointz,epsilon);

fx = fgx * hdeltasum;
fy = fgy * hdeltasum;
fz = fgz * hdeltasum;

% calculate stokes flow from the particle force

fhx = fftns(fx);
fhy = fftns(fy);
fhz = fftns(fz);

kdotfh = (k1.*fhx + k2.*fhy + k3.*fhz).*C2;

u = ifftns(fhx.*C1 - k1.*kdotfh);
v = ifftns(fhy.*C1 - k2.*kdotfh);
w = ifftns(fhz.*C1 - k3.*kdotfh);

if param.timestep == false
    px = pointx;
    py = pointy;
    pz = pointz;
    pu = pointu;
    pv = pointv;
    pw = pointw;
    gu = u;
    gv = v;
    gw = w;
    return;
end

%% SIMULATION TIME STEPPING
% timestep constants
t = 0;
plotnum = param.plotnum;
    nplt = nsteps/plotnum;
j = 1;

% plot specs
animate = param.animate;
if animate
    nlines = 8;
    I = 1:N/nlines:N;

%{
    vidfig = figure(1)
    frames = nsteps;
    F(frames) = struct('cdata',[],'colormap',[]);
%}
    % quiver cross section
    [cross2,cross3] = ndgrid(x(I),x(I));
    cross1 = 0.5*L*ones(nlines,nlines);

end

% data storing

pxdata = zeros(M,plotnum);
pydata = zeros(M,plotnum);
pzdata = zeros(M,plotnum);
pudata = zeros(M,plotnum);
pvdata = zeros(M,plotnum);
pwdata = zeros(M,plotnum);
gudata = zeros(N,N,N,plotnum);
gvdata = zeros(N,N,N,plotnum);
gwdata = zeros(N,N,N,plotnum);

% timing variables
t_point2grid = 0;
t_grid2point = 0;
t_fft = 0;
t_timestep = 0;
t_plot = 0;

GridIntp = griddedInterpolant(gridx,gridy,gridz,u);

% timestep!
for n = 1:nsteps
    t = t + dt;
    
    % calculate forcing function from particles
    tic;
    hdeltasum = deltasmear_tri(gridx,gridy,gridz,pointx,pointy,pointz,epsilon);
    
    fx = fgx * hdeltasum;
    fy = fgy * hdeltasum;
    fz = fgz * hdeltasum;
    
    t_point2grid = t_point2grid + toc;
    
    tic;
    % calculate flow velocity at particle positions
    GridIntp.Values = real(u);
    pointu = GridIntp(pointx,pointy,pointz);
    GridIntp.Values = real(v);
    pointv = GridIntp(pointx,pointy,pointz);
    GridIntp.Values = real(w);
    pointw = GridIntp(pointx,pointy,pointz);
    
    t_grid2point = t_grid2point + toc;
    
    % calculate stokes flow from the particle force
    tic;
    fhx = fftns(fx);
    fhy = fftns(fy);
    fhz = fftns(fz);
    
    kdotfh = (k1.*fhx + k2.*fhy + k3.*fhz).*C2;
    
    u = ifftns(fhx.*C1 - k1.*kdotfh);
    v = ifftns(fhy.*C1 - k2.*kdotfh);
    w = ifftns(fhz.*C1 - k3.*kdotfh);
    t_fft = t_fft + toc;
    tic;
    % timestep particles
    [pointx,pointy,pointz] = timestep_eulerfwd(L,dt,pointx,pointy,pointz,pointu,pointv,pointw);
    t_timestep = t_timestep + toc;
    
    
    % plot
    if animate
        clf;
        %view(-25,50);
        %view(-16,18);
        tic;
        plot3(pointx,pointy,pointz,'k.');
        view(-90,0);
        hold on;
        %quiver3(particlex1,particlex2,particlex3,particleu,particlev,particlew);
        %plot3(cross1,cross2,cross3,'r*');
        hsurfaces = slice(gridy,gridx,gridz,...
            abs(u.*u+v.*v+w.*w),...
            .5*L,L*(N-1)/N,0);
        set(hsurfaces,'FaceColor','interp','EdgeColor','none')
        colormap jet
        %{
        quiver3(cross1,cross2,cross3,...
        squeeze(u(N/2,I,I)),...
        squeeze(v(N/2,I,I)),...
        squeeze(w(N/2,I,I)),'r');
        %}
        t_plot = t_plot + toc;
        axis([0 L 0 L 0 L]);
        title(num2str(n));
        drawnow;
        %F(j) = getframe(vidfig);
    end
    
    % save data
    if mod(n,nplt) == 0
        pxdata(:,j) = pointx;
        pydata(:,j) = pointy;
        pzdata(:,j) = pointz;
        pudata(:,j) = pointu;
        pvdata(:,j) = pointv;
        pwdata(:,j) = pointw;
        gudata(:,:,:,j) = u;
        gvdata(:,:,:,j) = v;
        gwdata(:,:,:,j) = w;
        j = j+1;
    end
end

%% REPORT

t_point2grid = t_point2grid/n
t_grid2point = t_grid2point/n
t_fft = t_fft/n
t_timestep = t_timestep/n
t_plot = t_plot*nplt/n

%{
px = particlex1;
py = particlex2;
pz = particlex3;
pu = particleu;
pv = particlev;
pw = particlew; 
%}
px = pxdata;
py = pydata;
pz = pzdata;
pu = pudata;
pv = pvdata;
pw = pwdata;
gu = gudata;
gv = gvdata;
gw = gwdata;