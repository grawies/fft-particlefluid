function [error, varyVals, velnorm] = CompareInitialFlows(varyParam, iterations, epsilonFix, deltaType, p)
    
    if ~ strcmp(varyParam, 'N')
        epsilonFix = true;
    end
    
    LocalInit(0);
    
    N0 = 16;
    N = N0;
    L0 = 2*pi;
    L = L0;
    
    imax = iterations;
    
    udata = zeros(N0*N0*N0,imax);
    vdata = udata;
    wdata = udata;
    NVals = zeros(1,imax);
    velnorm = zeros(1,imax-1);
    I = 1:1:N0;
    
    %% SETUP VARYING VALUES
    for i = 1:imax
        NVals(i) = N0 * 2^(i-1);
        NVals(i) = N0 * i;
    end
    
    %% DATA COLLECTION
    for i = 1:imax
        
        N = NVals(i)
        
        if strcmp(varyParam, 'N')
            I = 1:N/N0:N;
        end
        if strcmp(varyParam, 'L')
            L = L0 * N / N0;
        end
        epsilon = 2*pi/N;
        if epsilonFix
            epsilon = 2*pi/16;
        end
        
        [u,v,w] = CalcInitialFlow(L,N,deltaType,[p(1)*i,p(2)]);
        if strcmp(varyParam,'N')
            utmp = u(I,I,I);
            vtmp = v(I,I,I);
            wtmp = w(I,I,I);
        end
        if strcmp(varyParam,'L')
            % identify elements in the box [L/2 +- L0/2]^3
            x = (0:N-1)*L/N;
            [x1,x2,x3] = ndgrid(x,x,x);
            
            b = (abs(x1(:)-L/2+0.001*L0) <= L0/2) .* ...
                (abs(x2(:)-L/2+0.001*L0) <= L0/2) .* ...
                (abs(x3(:)-L/2+0.001*L0) <= L0/2);
            [c,sortorder] = sort(b);
            I = sortorder(end-sum(c)+1:end);
            
            utmp = u(:); utmp = utmp(I);
            vtmp = v(:); vtmp = vtmp(I);
            wtmp = w(:); wtmp = wtmp(I);
        end
        
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
    
    if varyParam=='N'
        varyVals = NVals;
    end
    if varyParam=='L'
        varyVals = L0*NVals/N0;
    end
    
    %% FANCY PLOT
    fontsize = 12;
    loglog(varyVals, error, 'k-');
    hold on;
    loglog(varyVals, error, 'kx');
    grid on;
    set(gca,'FontSize',fontsize);
    ylabel('Normed 2-norm of velocity differences');
    
    if varyParam == 'N'
        xlabel('Grid size N');
        title('Convergence of velocity field solution, epsilon ~ 1/N');
        if epsilonFix
            title('Convergence of velocity field solution, epsilon held fixed');
        end
    end
    if varyParam == 'L'
        xlabel('Cube side length L');
        title('Convergence of velocity field solution');
    end
    
end