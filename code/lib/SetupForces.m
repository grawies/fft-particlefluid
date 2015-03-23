function F = SetupForces(S)

% gravity force on point particles
frel = S.mrel * S.g;
F.fgx = frel * 0;
F.fgy = frel *  0;
F.fgz = frel * -1;

