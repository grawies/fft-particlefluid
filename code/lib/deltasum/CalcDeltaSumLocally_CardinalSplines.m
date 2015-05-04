function hdeltasum = CalcDeltaSumLocally_CardinalSplines(S,P,F)
    hdeltasum = zeros(S.N,S.N,S.N);
    
    p = S.splineOrder;
    hinv = 1.0/S.h;
    hmultinv = 1.0/S.hmult;
    
    ivals = S.dIvals;
    [di,dj,dk] = ndgrid(ivals);
    N = S.N;
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
        idxi = i0 + ivals;
        idxj = j0 + ivals;
        idxk = k0 + ivals;
        
        % this code gets a cube of indices within epsilon of particle m
        i = i0-1+di;
        j = j0-1+dj;
        k = k0-1+dk;
        
        % restrict indices to { 1,2,...,N }
        % the efficiency makes you cringe, yes.
        % in the deepest level of nested for-loops,
        % I have placed an 'if'.
         % #NoToVectorization #FightThePower
        if (i0 + S.dI > N ||...
            i0 - S.dI < 1 )
            i = i + N * (i < 1) - N * (i > N);
            idxi = idxi + N * (idxi < 1) - N * (idxi > N);
        end
        if (j0 + S.dI > N ||...
            j0 - S.dI < 1 )
            j = j + N * (j < 1) - N * (j > N);
            idxj = idxj + N * (idxj < 1) - N * (idxj > N);
        end
        if (k0 + S.dI > N ||...
            k0 - S.dI < 1 )
            k = k + N * (k < 1) - N * (k > N);
            idxk = idxk + N * (idxk < 1) - N * (idxk > N);
        end
        Mx = CardinalSpline(i-P.x1(m)*hinv, p);
        My = CardinalSpline(j-P.x2(m)*hinv, p);
        Mz = CardinalSpline(k-P.x3(m)*hinv, p);
        %{
        Mx = CardinalSpline(hmultinv*(i-P.x1(m)*hinv), p);
        My = CardinalSpline(hmultinv*(j-P.x2(m)*hinv), p);
        Mz = CardinalSpline(hmultinv*(k-P.x3(m)*hinv), p);
        %}
        
        hdeltasum(idxi,idxj,idxk) = hdeltasum(idxi,idxj,idxk)...
            + Mx.*My.*Mz;
        
    end
    hdeltasum = hdeltasum * hmultinv^3;
end 
