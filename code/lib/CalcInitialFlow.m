function [u,v,w] = CalcInitialFlow(L,N,epsilon,deltaType,p)

    LocalInit(0);


    [S,F] = SetupWorld(L, 0.4,...
                    N, 300,... % N,M
                    0.01, 1000,... % dt, nsteps
                    deltaType, p);

    P = SetupParticles(S, 'sphere', 0.5, 0.0);
    
    if strcmp(deltaType, 'hat')
        F.hdeltasum = CalcDeltaSumLocally(S, P, F);
    end
    if strcmp(deltaType, 'spline')
        F.hdeltasum = CalcDeltaSumNaive_CardinalSplines(S, P, F);
    end
    
    [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);

    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);

    u = S.u;
    v = S.v;
    w = S.w;

end
