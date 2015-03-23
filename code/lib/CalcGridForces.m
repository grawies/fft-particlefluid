function [fx,fy,fz] = CalcGridForces(S,P,F)

    fx = F.fgx * F.hdeltasum;
    fy = F.fgy * F.hdeltasum;
    fz = F.fgz * F.hdeltasum;

end