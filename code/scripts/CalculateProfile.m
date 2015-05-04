%% PRELIMINARIES

Init(0);

L=2*pi;
N = 32;
epsilon = L/16;

[S,F] = SetupWorld(L, 0.4,... 
                N, 8,... % N,M
                0.01, 1000,... % dt, nsteps
                'spline',[1,4]);

P = SetupParticles(S, 'line', 0.5, 0.0);

% plot specs
Plt = SetupPlotting(S, 8, 50);

%% DO (...WHILE)

F.hdeltasum = CalcDeltaSumLocally_CardinalSplines(S, P, F);
%F.hdeltasum = CalcDeltaSumLocally(S, P, F);

[F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);

% calculate stokes flow from the particle force
[S.u, S.v, S.w] = SolveStokes(S,P,F);


uvals = zeros(Plt.nplots,S.N);
ppos = zeros(Plt.nplots,S.M);

pltidx = 1;
%% ...WHILE: SIMULATION TIME STEPPING 

for n = 1:S.nmax
    S.t = S.t + S.dt;
    
    % calculate forcing function from particles
    F.hdeltasum = CalcDeltaSumLocally_CardinalSplines(S, P, F);
    %F.hdeltasum = CalcDeltaSumLocally(S,P,F);
    
    % calculate flow velocity at particle positions
    P = InterpolateGridToParticles(S,P);
    [F.fx,F.fy,F.fz] = CalcGridForces(S,P,F);
   
    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);
    
    % timestep particles
    [S,P] = TimestepFwdEuler(S,P);
    
    % plot
    if mod(n,Plt.interval) == 0
        uvals(pltidx,:) = S.w(N/2+1, N/2+1, :);
        ppos(pltidx,:) = P.x3(:);
        pltidx = pltidx+1;
    end
end

%% PLOT

w = waterfall(abs(uvals));
axis([]);
set(w,'edgecolor','k');
hold on;
plot3(ppos*S.N/S.L, repmat((1:Plt.nplots)',1,S.M), zeros(Plt.nplots, S.M), 'k*');