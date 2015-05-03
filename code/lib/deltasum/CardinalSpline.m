function splinesum = CardinalSpline (u,p)
    
    splinesum = 0*u;
    
    for k = 0:1:p
        splinesum = splinesum + (-1)^k * nchoosek(p,k) * (0.5*(u-k+p*0.5 + abs(u-k+p*0.5))).^(p-1);
    end
    
    splinesum = splinesum / factorial(p-1);
end