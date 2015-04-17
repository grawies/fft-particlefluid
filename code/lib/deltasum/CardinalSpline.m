function splinesum = CardinalSpline (p,u)
    
    splinesum = 0;
    
    for k = 0:1:p
        splinesum = splinesum + (-1)^k * nchoosek(p,k) * (0.5*(u-k-p*0.5 + abs(u-k-p*0.5))).^(p-1);
    end
    
    %splinesum = (-1).^k * nchoosek
    
    splinesum = splinesum / factorial(p-1);
end