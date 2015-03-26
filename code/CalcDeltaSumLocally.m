function hdeltasum = CalcDeltaSumLocally(S,P,F)

    hdeltasum = zeros(S.N,S.N,S.N);
    
    
    ivals = S.dIvals;
    [di,dj,dk] = ndgrid(ivals,ivals,ivals);
    N = S.N;
    i0vals = 1 + round(S.N/S.L * P.x1);
    j0vals = 1 + round(S.N/S.L * P.x2);
    k0vals = 1 + round(S.N/S.L * P.x3);
        
    for m = 1:S.M
        % for each particle
        % add delta sum locally
        i0 = i0vals(m);
        j0 = j0vals(m);
        k0 = k0vals(m);
        
        idxi = i0+ivals;
        idxj = j0+ivals;
        idxk = k0+ivals;
        
        i = i0+di;
        j = j0+dj;
        k = k0+dk;
        
        if (i0 + S.dI > N ||...
            j0 + S.dI > N ||...
            k0 + S.dI > N ||...
            i0 - S.dI < 1 ||...
            j0 - S.dI < 1 ||...
            k0 - S.dI < 1 )
            i = i + N * (i < 1) - N * (i > N);
            j = j + N * (j < 1) - N * (j > N);
            k = k + N * (k < 1) - N * (k > N);
            idxi = idxi + N * (idxi < 1) - N * (idxi > N);
            idxj = idxj + N * (idxj < 1) - N * (idxj > N);
            idxk = idxk + N * (idxk < 1) - N * (idxk > N);
        end
        
        x = S.L * (i-1) / S.N;
        y = S.L * (j-1) / S.N;
        z = S.L * (k-1) / S.N;
        
        hdeltasum(idxi,idxj,idxk) = hdeltasum(idxi,idxj,idxk) + ...
            max(1-abs(x-P.x1(m))/S.epsilon,0).*...
            max(1-abs(y-P.x2(m))/S.epsilon,0).*...
            max(1-abs(z-P.x3(m))/S.epsilon,0);
            
    end
    
end