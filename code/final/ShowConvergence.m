function [varyVals, error] = ShowConvergence(N0, imax, L, g)
    
%L=2*pi; g = @(x,y,z) sin(L*cos(x*2*pi/L)).*cos(y*2*pi/L) + sin(z*2*pi/L);
%L=2*pi; g = @(x,y,z) sin(x*2*pi/L).*cos(y*2*pi/L).*(2/L^2 * z.^3 - 3/L * z.^2 + z);

    LocalInit(0);
    
    L = 2*pi;
    
    udata = zeros(N0*N0*N0,imax);
    vdata = udata;
    wdata = udata;
    NVals = zeros(1,imax);
    velnorm = zeros(1,imax-1);
    I = 1:1:N0;
    
    %% SETUP VARYING VALUES
    for i = 1:imax
        NVals(i) = N0 * i;
    end
    
    %% DATA COLLECTION
    for i = 1:imax
        
        N = NVals(i);
        
        I = 1:N/N0:N;
        
        [S,F] = SetupWorld(L, 1,... 
                        N, 8,... % N, M
                        0.03, 5,... % dt, nsteps
                        'spline', 4);
        P = SetupParticles(S, 'line', 0.5, 0.0);
        
        F.hdeltasum = g(S.x1,S.x2,S.x3);
        [F.fx, F.fy, F.fz] = CalcGridForces(S, P, F);
        [u,v,w] = SolveStokes(S,P,F);
        
        % report
        %z = u.*u+v.*v+w.*w;
        %sum(z(:)/N^3)
        
        utmp = u(I,I,I);
        vtmp = v(I,I,I);
        wtmp = w(I,I,I);
        
        udata(:,i) = utmp(:);
        vdata(:,i) = vtmp(:);
        wdata(:,i) = wtmp(:);
    end
    
    %% DATA PROCESSING
    
    error = zeros(1,imax-1);
    
    for i = 1:imax-1
        % difference (numerator)
        udiff = udata(:,i+1) - udata(:,i);
        vdiff = vdata(:,i+1) - vdata(:,i);
        wdiff = wdata(:,i+1) - wdata(:,i);
        diffsq = udiff.^2+vdiff.^2+wdiff.^2;
        diffsq = diffsq(:);
        
        % norm (denominator)
        u = udata(:,i+1);
        v = vdata(:,i+1);
        w = wdata(:,i+1);
        normsq = u.^2+v.^2+w.^2;
        normsq = normsq(:);
        
        velnorm(i) = sqrt(sum(normsq));
        error(i) = sqrt(sum(diffsq) / sum(normsq));
    end
    NVals = NVals(2:end);
    
    varyVals = NVals;
    
    %% FANCY PLOT
    fontsize = 12;
    loglog(varyVals, error, 'k--','LineWidth',2);
    hold on;
    loglog(varyVals, error, 'k+','MarkerSize',12);
    grid on;
    set(gca,'FontSize',fontsize);
    ylabel('relative 2-norm error in velocity')
    xlabel('Grid size N');
    title('convergence of solution for smooth forcing function')

    %axis([10^.75 10^2.25 10^(-17) 10^3])
end