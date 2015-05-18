function [varyVals, error] = ShowSplineConvergence(N0, imax, M, dparams,pltL,pltP)

    %LocalInit(0);
    
    L = 2.0;
    
    udata = zeros(N0*N0*N0,imax);
    vdata = udata;
    wdata = udata;
    NVals = zeros(1,imax);
    velnorm = zeros(1,imax-1);
    I = 1:1:N0;
    
    dpar1 = dparams(1);
    dpar2 = dparams(2);
    
    %% SETUP VARYING VALUES
    for i = 1:imax
        NVals(i) = N0 * i;
    end
    
    plt = {'r','g','b','y','k','m','c'};
    %% DATA COLLECTION
    for i = 1:imax
        N = NVals(i);
        
        I = 1:N/N0:N;
        [S,F] = SetupWorld(L, 1,... 
                N, M,... % N, M
                0.01, 100,... % dt, nsteps
                'spline', [dpar1*i,dpar2]); % splinetype, deltaparams[2][1]
           
        P = SetupParticles(S, 'line', 0.5, 0.0);
        zList = L*.5 + [-.1 .1 .2 -.4];
        P.x3 = zList(1:M);
        
        F.hdeltasum = CalcDeltaSumLocally(S, P, F);
        %plot(squeeze(F.hdeltasum(N/2,N/2,I)),plt{i}); hold on;
        
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
    return
    %% FANCY PLOT
    
    figure(3)
    
    fontsize = 12;
    loglog(varyVals, error, pltL,'LineWidth',2);
    hold on;
    loglog(varyVals, error, pltP,'MarkerSize',12);
    grid on;
    set(gca,'FontSize',fontsize);
    ylabel('relative 2-norm error in velocity')
    xlabel('Grid size N');
    title('convergence of solution for smooth forcing function')

    %axis([10^.75 10^2.25 10^(-17) 10^3])
end