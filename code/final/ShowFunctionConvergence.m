function [varyVals, error] = ShowFunctionConvergence(N0, imax)
    L = 1.0;
    g = @(x,y,z) sin(2*pi*cos(x*2*pi/L)).*cos(y*2*pi/L) + sin(z*2*pi/L);
    [varyVals, error] = ShowConvergence(N0, imax, L, g);
end