function [S,F] = SetupWorld(L,nu,N,M,dt,nsteps, deltaType, deltaParam)

    S.L = L;
    S.nu = nu;
    S.g = -9.82;
    S.mrel = 0.1;

    S.N = N;
    S.M = M;

    S.h = S.L/S.N;
    
    S.deltaType = deltaType;
    if strcmp(deltaType, 'hat')
        S.epsilon = deltaParam;
    end
    if strcmp(deltaType, 'spline')
        S.splineOrder = deltaParam;
        S.epsilon = S.splineOrder/2 * S.h;
    end
    
    S.dI = round( S.epsilon * N / L ) + 1; % +1 just in case
    S.dIvals = -S.dI:S.dI;
    
    % time values
    S.t = 0;
    S.dt = dt;
    S.nsteps = nsteps;
    S.tmax = dt*nsteps;
    S.nmax = round(S.tmax/S.dt);

    % space vectors
    S.x = (0:N-1)*L/N;
    [S.x1,S.x2,S.x3] = ndgrid(S.x);

    % space interpolation
    
    S.GridInterp = griddedInterpolant(...
        S.x1,S.x2,S.x3,...
        zeros(N,N,N));
    
    % k-space vectors
    S.k = -N/2 : N/2-1;

    [S.k1,S.k2,S.k3] = ndgrid(S.k);
    S.ksq = S.k1.^2+S.k2.^2+S.k3.^2;

    S.ksqinv = 1./S.ksq;
    S.ksqinvsq = 1./(S.ksq.^2);

    S.ksqinv(N/2+1,N/2+1,N/2+1) = 0;
    S.ksqinvsq(N/2+1,N/2+1,N/2+1) = 0;

    % fourier solver constants
    S.C1 = (S.L /2/pi)^2 * S.ksqinv / S.nu;
    S.C2 = (S.L /2/pi)^2 * S.ksqinvsq / S.nu;

    % setup forces
    F = SetupForces(S);
end
