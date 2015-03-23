function S = SetupWorld(L,nu,N,M,dt,nsteps, epsilon)

    S.L = L;
    S.nu = nu;
    S.g = -9.82;
    S.mrel = 0.1;

    S.N = N;
    S.M = M;

    S.epsilon = epsilon;

    % time constants
    S.dt = dt;
    S.nsteps = nsteps;
    S.tmax = dt*nsteps;
    S.nmax = round(S.tmax/S.dt);

    % space vectors
    S.x = (0:N-1)*L/N;
    [S.x1,S.x2,S.x3] = ndgrid(S.x,S.x,S.x);

    % k-space vectors
    S.k = -N /2:N/2-1;

    [S.k1,S.k2,S.k3] = ndgrid(S.k,S.k,S.k);
    S.ksq = S.k1.^2+S.k2.^2+S.k3.^2;

    S.ksqinv = 1./S.ksq;
    S.ksqinvsq = 1./(S.ksq.^2);

    S.ksqinv(N/2+1,N/2+1,N/2+1) = 0;
    S.ksqinvsq(N/2+1,N/2+1,N/2+1) = 0;

    % fourier solver constants
    S.C1 = S.ksqinv / S.nu;
    S.C2 = S.ksqinvsq / S.nu;

end