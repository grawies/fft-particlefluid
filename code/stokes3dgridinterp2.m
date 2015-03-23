%% PRELIMINARIES

addpath('lib');
addpath('lib/fft_shifted');

rng(0);
L=2*pi;
N = 16;
S = SetupWorld(L, N,...
                16, 300,...
                0.01, 200,...
                L/N);

P = SetupParticles(S, 'sphere', 0.5, 0.0);

%% FUNTIONS, EQUATION PROPERTIES

F = SetupForces(S);

% delta function half-width

%% DO (...WHILE)

F.hdeltasum = CalcDeltaSum(S, P, F);

[F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);

F.fhx = fftns(F.fx);
F.fhy = fftns(F.fy);
F.fhz = fftns(F.fz);

% calculate stokes flow from the particle force
[S.u, S.v, S.w] = SolveStokes(S,P,F);

%% SIMULATION TIME STEPPING
% timestep constants
t=0;
% plot specs
plotnum = 500; nplt = S.nmax/plotnum;
nlines = 8;
I = 1:S.N/nlines:S.N;

vidfig = figure(1);
frames = S.nsteps;
Fv(frames) = struct('cdata',[],'colormap',[]);
j=1;

% streamlines
[sx,sy] = ndgrid(S.x(I),S.x(I));
sz = S.L*ones(nlines,nlines);
% quiver cross section
[cross2,cross3] = ndgrid(S.x(I),S.x(I));
cross1 = 0.49*S.L*ones(nlines,nlines);

% timing variables
t_interpol = 0;
t_fft = 0;
t_timestep = 0;
t_plot = 0;

GridIntp = griddedInterpolant(S.x1,S.x2,S.x3,u);

% timestep!
for n = 1:S.nmax
    t = t + S.dt;
    
    % calculate forcing function from particles
    hdeltasum = 0;
    tic;
    for i = 1:S.M
        hdelta1 = max(1-abs(S.x1-P.x1(i))/epsilon,0);
        hdelta2 = max(1-abs(S.x2-P.x2(i))/epsilon,0);
        hdelta3 = max(1-abs(S.x3-P.x3(i))/epsilon,0);
        hdelta = hdelta1.*hdelta2.*hdelta3;
        hdeltasum = hdeltasum + hdelta;
    end
    
    % calculate flow velocity at particle positions
    GridIntp.Values = real(u);
    P.u = GridIntp(P.x1,P.x2,P.x3);
    GridIntp.Values = real(v);
    P.v = GridIntp(P.x1,P.x2,P.x3);
    GridIntp.Values = real(w);
    P.w = GridIntp(P.x1,P.x2,P.x3);
    
    fx = fgx * hdeltasum;
    fy = fgy * hdeltasum;
    fz = fgz * hdeltasum;
   
    t_interpol = t_interpol + toc;
    
    % calculate stokes flow from the particle force
    tic;
    fhx = fftns(fx);
    fhy = fftns(fy);
    fhz = fftns(fz);
    
    kdotfh = (S.k1.*fhx + S.k2.*fhy + S.k3.*fhz).*S.C2;
    
    u = ifftns(fhx.*S.C1 - S.k1.*kdotfh);
    v = ifftns(fhy.*S.C1 - S.k2.*kdotfh);
    w = ifftns(fhz.*S.C1 - S.k3.*kdotfh);
    t_fft = t_fft + toc;
    tic;
    % timestep particles
    P.x1 = P.x1 + S.dt * P.u;
    P.x1 = P.x1 + S.L*(P.x1<0) - S.L*(P.x1>S.L);
    P.x2 = P.x2 + S.dt * P.v;
    P.x2 = P.x2 + S.L*(P.x2<0) - S.L*(P.x2>S.L);
    P.x3 = P.x3 + S.dt * P.w;
    P.x3 = P.x3 + S.L*(P.x3<0) - S.L*(P.x3>S.L);
    t_timestep = t_timestep + toc;
    % plot
    if mod(n,nplt) == 0
        
        clf;
        % STREAMLINES
        %streamline(stream3(x1,x2,x3,u,v,w,sx,sy,sz,[0.1]));
        
        %view(-25,50);
        %view(-16,18);
        tic;
        plot3(P.x1,P.x2,P.x3,'k.');
        view(-90,0);
        hold on;
        %quiver3(P.x1,P.x2,P.x3,P.u,P.v,P.w);
        %plot3(cross1,cross2,cross3,'r*');
        wind_speed = sqrt(u.^2 + v.^2 + w.^2);
        hsurfaces = slice(S.x2,S.x1,S.x3,...
            abs(u.*u+v.*v+w.*w),...
            .5*S.L,S.L*(S.N-1)/S.N,0);
        set(hsurfaces,'FaceColor','interp','EdgeColor','none')
        colormap jet
        
        quiver3(cross1,cross2,cross3,...
        squeeze(u(S.N/2,I,I)),...
        squeeze(v(S.N/2,I,I)),...
        squeeze(w(S.N/2,I,I)),'r');
        
        t_plot = t_plot + toc;
        axis([0 S.L 0 S.L 0 S.L]);
        %axis([-.01*L L*1.01 -.01*L L*1.01 -.01*L L*1.01]);
        %
        % might have to pause if:
        %UIJ_AreThereWindowShowsPending - timeout waiting for window to show up
        drawnow;
        Fv(j) = getframe(vidfig);
        j=j+1;
    end
end

%% REPORT

t_interpol/n
t_fft/n
t_timestep/n
t_plot*nplt/n