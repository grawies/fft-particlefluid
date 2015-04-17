function CompareRegularizationFunctions(pmax, N)

    Init(0);

    L = 2*pi;
    epsilon = L/N;
    
    [S,F] = SetupWorld(L, 0.4,...
                    N, 300,... % N,M
                    0.01, 1000,... % dt, nsteps
                    epsilon);

    P = SetupParticles(S, 'sphere', 0.5, 0.0);

    F.hdeltasum = CalcDeltaSum(S, P, F);

    [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);
    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);