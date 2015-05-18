function F = SetupForces(S)

% gravity force on point particles
F.frel = S.mrel * S.g;
F.fgx = F.frel * 0;
F.fgy = F.frel *  0;
F.fgz = F.frel * 1;

