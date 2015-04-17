function hdeltasum = CalcDeltaSumNaive_CardinalSplines(S,P,F)
% this is in fact the naive way

    hdeltasum = zeros(S.N,S.N,S.N);
    
    p = S.splineOrder;
    h = S.L/S.N;
    
    for m = 1:S.M
        % for each particle
        % add delta sum naively
        Mx = CardinalSpline(p, (P.x1(m)-S.x1)/h);
        My = CardinalSpline(p, (P.x2(m)-S.x2)/h);
        Mz = CardinalSpline(p, (P.x3(m)-S.x3)/h);
        hdeltasum = hdeltasum + Mx.*My.*Mz;
    end
    
end 
