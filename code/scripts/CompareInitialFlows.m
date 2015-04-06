function [norms, varyVals] = CompareInitialFlows(varyParam, iterations, epsilonFix)
    
    Init(0);
    
    N0 = 16;
    N = N0;
    L0 = 2*pi;
    L = L0;
    
    imax = iterations; % default is 5, 6 too many for N
    
    udata = zeros(N0*N0*N0,imax);
    vdata = udata;
    wdata = udata;
    varyVals = zeros(1,imax);
    I = 1:1:N0;
    
    %% DATA COLLECTION
    for i = 1:imax
        if varyParam == 'N'
            N = N0 * 2^(i-1);
            varyVals(i) = N
            I = 1:N/N0:N;
        end
        if varyParam == 'L'
            L = L0 * 2^(i-1);
            varyVals(i) = L;
            N = N0 * L/L0
            epsilonFix = true;
        end
        
        [u,v,w] = CalcInitialFlow(L,N,epsilonFix);
        if varyParam=='N' 
            utmp = u(I,I,I);
            vtmp = v(I,I,I);
            wtmp = w(I,I,I);
        end
        if varyParam=='L'
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
    
    norms = zeros(1,imax-1);
    
    for i = 1:imax-1
        udiff = udata(:,i+1) - udata(:,i);
        vdiff = vdata(:,i+1) - vdata(:,i);
        wdiff = wdata(:,i+1) - wdata(:,i);
        diffsq = udiff.^2+vdiff.^2+wdiff.^2;
        diffsq = diffsq(:);
        norms(i) = sqrt(sum(diffsq));
    end
    
    %% FANCY PLOT
    fontsize = 12;
    loglog(varyVals(1:end-1), norms, 'k-');
    hold on;
    loglog(varyVals(1:end-1), norms, 'kx');
    grid on;
    set(gca,'FontSize',fontsize);
    ylabel('2-norm of velocity differences');
    
    if varyParam == 'N'
        xlabel('Grid size N');
        title('Convergence of velocity field solution, epsilon ~ 1/N');
        if epsilonFix
            title('Convergence of velocity field solution, epsilon held fixed');
        end
    end
    if varyParam == 'L'
        xlabel('Grid size L');
        title('Convergence of velocity field solution');
    end
    
end