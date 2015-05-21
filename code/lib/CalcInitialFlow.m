function [u,v,w] = CalcInitialFlow(L,N,deltaType,p)

    LocalInit(0);


    [S,F] = SetupWorld(L, 1.0,...
                    N, 300,... % N,M
                    0.01, 1000,... % dt, nsteps
                    deltaType, p);

    P = SetupParticles(S, 'sphere', 0.5, 0.0);
    
    F.hdeltasum = CalcDeltaSumLocally(S, P, F);

    
    [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);

    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);

    u = S.u;
    v = S.v;
    w = S.w;

end
