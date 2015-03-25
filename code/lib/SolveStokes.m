function [u,v,w] = SolveStokes(S,P,F)

    F.fhx = fftns(F.fx);
    F.fhy = fftns(F.fy);
    F.fhz = fftns(F.fz);
    
    kdotfhdiv = (S.k1.*F.fhx + S.k2.*F.fhy + S.k3.*F.fhz).*S.C2;

    u = ifftns(F.fhx.*S.C1 - S.k1.*kdotfhdiv);
    v = ifftns(F.fhy.*S.C1 - S.k2.*kdotfhdiv);
    w = ifftns(F.fhz.*S.C1 - S.k3.*kdotfhdiv);

end