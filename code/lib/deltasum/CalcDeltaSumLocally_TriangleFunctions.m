function hdeltasum = CalcDeltaSumLocally_TriangleFunctions(S,P,F)

    hdeltasum = zeros(S.N,S.N,S.N);
    
    
    ivals = S.dIvals;
    [di,dj,dk] = ndgrid(ivals);
    N = S.N;
    h = S.L/S.N;
    epsiloninvcube = 1 / S.epsilon^3;
    % get the indices corresponding to the grid point closest to each
    % particle
    i0vals = 1 + round(S.N/S.L * P.x1);
    j0vals = 1 + round(S.N/S.L * P.x2);
    k0vals = 1 + round(S.N/S.L * P.x3);
        
    for m = 1:S.M
        % for each particle
        % add its local contribution the the delta sum
        
        i0 = i0vals(m);
        j0 = j0vals(m);
        k0 = k0vals(m);
        
        % this code gets intervals of indices within epsilon of particle m
        % does this even mean anything?
        idxi = i0 + ivals;
        idxj = j0 + ivals;
        idxk = k0 + ivals;
        
        % this code gets a cube of indices within epsilon of particle m
        i = i0+di;
        j = j0+dj;
        k = k0+dk;
        
        % restrict indices to { 1,2,...,N }
        %{
        i = mod(i-1, N) + 1;
        j = mod(j-1, N) + 1;
        k = mod(k-1, N) + 1;
        idxi = mod(idxi-1, N) + 1;
        idxj = mod(idxj-1, N) + 1;
        idxk = mod(idxk-1, N) + 1;
        %}
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
        
        x = h * (i-1);
        y = h * (j-1);
        z = h * (k-1);
        
        hdeltasum(idxi,idxj,idxk) = hdeltasum(idxi,idxj,idxk) + ...
            epsiloninvcube * ...
            max(1-abs(x-P.x1(m))/S.epsilon,0).*...
            max(1-abs(y-P.x2(m))/S.epsilon,0).*...
            max(1-abs(z-P.x3(m))/S.epsilon,0);
            
    end
    
end