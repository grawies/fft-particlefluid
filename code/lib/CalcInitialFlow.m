function [u,v,w] = CalcInitialFlow(L,N,epsilonFix)

    Init(0);

    epsilon = 2*pi/N;
    if epsilonFix
        epsilon = 2*pi/16;
    end

    [S,F] = SetupWorld(L, 0.4,...
                    N, 300,... % N,M
                    0.01, 1000,... % dt, nsteps
                    epsilon);

    P = SetupParticles(S, 'sphere', 0.5, 0.0);

    F.hdeltasum = CalcDeltaSum(S, P, F);

    [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);

    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);

    P = InterpolateGridToParticles(S,P);
    [F.fx,F.fy,F.fz] = CalcGridForces(S,P,F);

    u = real(S.u);
    v = real(S.v);
    w = real(S.w);

end