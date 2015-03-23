function hdeltasum = CalcDeltaSum(S,P,F)

    hdeltasum = 0;

    for i = 1:S.M
        hdelta1 = max(1-abs(S.x1-P.x1(i))/S.epsilon,0);
        hdelta2 = max(1-abs(S.x2-P.x2(i))/S.epsilon,0);
        hdelta3 = max(1-abs(S.x3-P.x3(i))/S.epsilon,0);
        hdelta = hdelta1.*hdelta2.*hdelta3;
        hdeltasum = hdeltasum + hdelta;
    end
    
end