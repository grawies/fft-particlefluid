%% PRELIMINARIES

addpath('lib');
addpath('lib/fft_shifted');
addpath('lib/plotting');

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

% calculate stokes flow from the particle force
[S.u, S.v, S.w] = SolveStokes(S,P,F);

%% SIMULATION TIME STEPPING
% plot specs
Plt = SetupPlotting(S, 8, 10);

% timing variables
t_deltainterp = 0;
t_interpol = 0;
t_fft = 0;
t_timestep = 0;
t_plot = 0;

% timestep!
for n = 1:S.nmax
    S.t = S.t + S.dt;
    
    % calculate forcing function from particles
    tic;
    F.hdeltasum = CalcDeltaSum(S,P,F);
    t_deltainterp = t_deltainterp + toc;
    
    tic;
    % calculate flow velocity at particle positions
    P = InterpolateGridToParticles(S,P);
    [F.fx,F.fy,F.fz] = CalcGridForces(S,P,F);
    t_interpol = t_interpol + toc;
   
    % calculate stokes flow from the particle force
    tic;
    [S.u, S.v, S.w] = SolveStokes(S,P,F);
    t_fft = t_fft + toc;
    
    % timestep particles
    tic;
    [S,P] = TimestepFwdEuler(S,P);
    t_timestep = t_timestep + toc;
    
    % plot
    if mod(n,Plt.interval) == 0
        tic;
        
        PlotAll(S,P,Plt,[-90,0],...
            [false,true,true,true,true]...
            );
        CaptureFrame(Plt);
        
        t_plot = t_plot + toc;
    end
end

%% REPORT
times = [t_deltainterp/n,...
            t_interpol/n,...
            t_fft/n,...
            t_timestep/n,...
            t_plot/n
            ]
        