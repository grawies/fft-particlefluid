function [S,P,F] = CalcPreTime(L,N,epsilon, step)

    [S,F] = SetupWorld(L, 0.4,... 
                N, 1000,... % N,M
                0.01, 200,... % dt, nsteps
                epsilon);

    P = SetupParticles(S, 'sphere', 0.5, 0.0);

    % interpolate the point gravities to the grid
    F.hdeltasum = CalcDeltaSumLocally(S, P, F);
    if step==1
        return
    end

    % calculate particle force
    [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);
    % calculate stokes flow from the particle force
    [S.u, S.v, S.w] = SolveStokes(S,P,F);
    if step==2
        return
    end

    % interpolate the flow velocity from the grid to particles
    P = InterpolateGridToParticles(S,P);
    if step==3
        return
    end
end