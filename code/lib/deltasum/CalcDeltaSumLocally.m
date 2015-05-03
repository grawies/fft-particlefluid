function hdeltasum = CalcDeltaSumLocally(S,P,F)

    if strcmp(S.deltaType,'triangle')
        hdeltasum = CalcDeltaSumLocally_TriangleFunctions(S,P,F);
        return;
    end
    if strcmp(S.deltaType,'spline')
        hdeltasum = CalcDeltaSumLocally_CardinalSplines(S,P,F);
        return;
    end
    
end