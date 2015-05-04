function splinesum = CardinalSpline (u,p)
    
    splinesum = 0*u;
    
    for k = 0:1:p
        relpos = u-k+p*0.5;
        splinesum = splinesum + (-1)^k * nchoosek(p,k) * (0.5*(relpos + abs(relpos))).^(p-1);
    end
    
    splinesum = splinesum / factorial(p-1);
end