clear all
%% PRELIMINARIES

Init(0);

L=2*pi;
N = 64;
M = 1000;

[S,F] = SetupWorld(L, 1,... 
                N, M,... % N, M
                0.005, 800,... % dt, nsteps
                'triangle', [4 0]); % splinetype, deltaparams[2][1]

%P = SetupParticles(S, 'uniform', 0.5, 0.0);
P = SetupParticles(S, 'sphere', 0.5, 0.0);

% plot specs
Plt = SetupPlotting(S, 8, 1);

figure(1);

%% DO (...WHILE)


F.hdeltasum = CalcDeltaSumLocally(S, P, F);

[F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);

% calculate stokes flow from the particle force
[S.u, S.v, S.w] = SolveStokes(S,P,F);

%% ...WHILE: SIMULATION TIME STEPPING 

% timing variables
t_deltainterp = 0;
t_interpol = 0;
t_fft = 0;
t_timestep = 0;
t_plot = 0;

% timestep!

MOV(Plt.nplots) = struct('cdata',[],'colormap',[]);

for n = 1:S.nmax
    S.t = S.t + S.dt;
    
    % calculate forcing function from particles
    tic;
    F.hdeltasum = CalcDeltaSumLocally(S,P,F);
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
    
    prevpos.x = P.x1;
	prevpos.y = P.x2;
	prevpos.z = P.x3;
    
    % timestep particles
    tic;
    [S,P] = TimestepFwdEuler(S,P);
    t_timestep = t_timestep + toc;
    
    %}
    % plot
    
    if mod(n,20) == 0
        totTime = t_deltainterp+t_interpol+t_fft+t_timestep+t_plot
        expected = totTime * S.nmax / n
    end
    
    if mod(n,Plt.interval) == 0
        tic;
        set(Plt.vidfig, 'Position', [400 150 800 600])
        viewAng = [-38,28];
        %viewAng = [-90,0];
        PlotAll(S,prevpos,P,Plt,viewAng,...
            [false,true,true,false,false]...
            );
        MOV(Plt.framecount) = getframe(Plt.vidfig);
        Plt.framecount = Plt.framecount + 1;
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
        
%% SAVE MOVIE
%cd '~/Desktop';
%movie2avi(MOV,'test.avi','fps',20)
%plot3(a,b,c,'k*');