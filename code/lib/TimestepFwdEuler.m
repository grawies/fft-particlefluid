function [S,P] = TimestepFwdEuler(S,P)
    
    P.x1 = P.x1 + S.dt * P.u;
    P.x1 = P.x1 + S.L*(P.x1<0) - S.L*(P.x1>S.L);
    P.x2 = P.x2 + S.dt * P.v;
    P.x2 = P.x2 + S.L*(P.x2<0) - S.L*(P.x2>S.L);
    P.x3 = P.x3 + S.dt * P.w;
    P.x3 = P.x3 + S.L*(P.x3<0) - S.L*(P.x3>S.L);

end