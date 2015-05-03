%% PRELIMINARIES
function [x,p,f,u] = Get1DProfile(L,N,M,method,methodParam)
    LocalInit(0);
    [S,F] = SetupWorld(L, 1,... 
                    N, M,... % N,M
                    0.01, 1000,... % dt, nsteps
                    method, methodParam); % delta type

    P = SetupParticles(S, 'line', 0.5, 0.0);
    if M<5
        zList = L*.5 + [-.1 .1 .2 -.4];
        P.x3 = zList(1:M);
    end
    
    %% DO (...WHILE)

    F.hdeltasum = CalcDeltaSumLocally(S, P, F);
    (S.L/S.N)^3 * sum(F.hdeltasum(:))
    [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);
    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);

    x = S.x;
    p = P.x3;
    f = squeeze( F.fz(N/2+1, N/2+1, :) );
    u = squeeze( S.w(N/2+1, N/2+1, :) );
end
%[x,p,f,u] = Get1DProfile(1.0, 128, 2, 'spline', 4);