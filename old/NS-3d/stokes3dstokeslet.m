%% PRELIMINARIES
% constants
N = 16; dt = 0.0008;
M = 1000; % # of particles
nsteps=200;
tmax = dt*nsteps;
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
R0 = 1.0;
phi = unifrnd(0,2*pi, 1, M);
costheta = unifrnd(-1,1, 1, M);
u = unifrnd(0,1,1,M);
theta = acos( costheta );
r = R0 * u.^(1/3);

particlex1 = .5 * L + r .* sin( theta ) .* cos( phi );
particlex2 = .5 * L + r .* sin( theta ) .* sin( phi );
particlex3 = .5 * L + r .* cos( theta );

particleu = zeros(1,M);
particlev = zeros(1,M);
particlew = zeros(1,M);

%% FUNTIONS, EQUATION PROPERTIES

% gravity force on point particles
g = 9.82;
mrel = 0.1;
frel = mrel * g;
fgx = frel * 0;
fgy = frel * 0;
fgz = frel * -1;

%% SIMULATION TIME STEPPING
% timestep constants
nmax = round(tmax/dt); t = 0;
% plot specs
plotnum = 500; nplt = nmax/plotnum;
nlines = 8;
I = 1:N/nlines:N;

vidfig = figure(1);
frames = nsteps;
F(frames) = struct('cdata',[],'colormap',[]);
frameindex=1;

% streamlines
[sx,sy] = ndgrid(x(I),x(I));
sz = L*ones(nlines,nlines);
% quiver cross section
[cross1,cross3] = ndgrid(x(I),x(I));
cross2 = 0.5*L*ones(nlines,nlines);

% timing variables
t_interpol = 0;
t_fft = 0;
t_timestep = 0;
t_plot = 0;

distancex = zeros(M,M);
distancey = zeros(M,M);
distancez = zeros(M,M);
distancei = zeros(M,M);
distancei3 = zeros(M,M);
lead = 14;
% timestep!
for n = 1:nmax
    t = t + dt;
    
    % calculate particle distances
    tic;
    for i = 1:M
        for j = 1:i-1
            Dx = particlex1(i)-particlex1(j);
            distancex(i,j) = Dx;
            distancex(j,i) = -Dx;
            Dy = particlex2(i)-particlex2(j);
            distancey(i,j) = Dy;
            distancey(j,i) = -Dy;
            Dz = particlex3(i)-particlex3(j);
            distancez(i,j) = Dz;
            distancez(j,i) = -Dz;
            if i==j
                D = Inf;
            else
                D = sqrt(Dx^2+Dy^2+Dz^2);
            end
            distancei(i,j) = 1/D;
            distancei(j,i) = 1/D;
            D3 = D^3;
            distancei3(i,j) = 1/D3;
            distancei3(j,i) = 1/D3;
        end
    end
    % calculate particle velocity
    for i = 1:M
        dx = distancex(i,:);
        dy = distancey(i,:);
        dz = distancez(i,:);
        di = distancei(i,:);
        di3 = distancei3(i,:);
        dot = dx.*fgx+dy.*fgy+dz.*fgz;
        dotdiv = dot.*di3;
        particleu(i) = sum(fgx.*di + dx .* dotdiv);
        particlev(i) = sum(fgy.*di + dy .* dotdiv);
        particlew(i) = sum(fgz.*di + dz .* dotdiv);
        
    end
    
    t_fft = t_fft + toc;

    tic;
    % timestep particles
    particlex1 = particlex1 + dt * particleu - dt * particleu(lead);
    particlex2 = particlex2 + dt * particlev - dt * particlev(lead);
    particlex3 = particlex3 + dt * particlew - dt * particlew(lead);
    t_timestep = t_timestep + toc;
    % plot
    if mod(n,nplt) == 0
        
        clf;
        % STREAMLINES
        %streamline(stream3(x1,x2,x3,u,v,w,sx,sy,sz,[0.1]));
        
        %view(-25,50);
        view(0,0);
        tic;
        plot3(particlex1,particlex2,particlex3,'k.');
        hold on;
        %quiver3(particlex1,particlex2,particlex3,particleu,particlev,particlew);
        %plot3(cross1,cross2,cross3,'r*');
        %wind_speed = sqrt(u.^2 + v.^2 + w.^2);
        %hsurfaces = slice(x2,x1,x3,...
        %    abs(u.*u+v.*v+w.*w),...
        %    .5*L,L*(N-1)/N,0);
        %set(hsurfaces,'FaceColor','interp','EdgeColor','none')
        %colormap jet
        t_plot = t_plot + toc;
        axis([0 L 0 L 0 L]);
        %axis([-.01*L L*1.01 -.01*L L*1.01 -.01*L L*1.01]);
        %
        % might have to pause if:
        %UIJ_AreThereWindowShowsPending - timeout waiting for window to show up
        drawnow;
        F(frameindex) = getframe(vidfig);
        frameindex=frameindex+1;
    end
end

%% REPORT

t_interpol/n
t_fft/n
t_timestep/n
t_plot*nplt/n