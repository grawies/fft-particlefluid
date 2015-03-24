function P = InterpolateGridToParticles(S,P)
    
    S.GridInterp.Values = real(S.u);
    P.u = S.GridInterp(P.x1,P.x2,P.x3);
    S.GridInterp.Values = real(S.v);
    P.v = S.GridInterp(P.x1,P.x2,P.x3);
    S.GridInterp.Values = real(S.w);
    P.w = S.GridInterp(P.x1,P.x2,P.x3);

end