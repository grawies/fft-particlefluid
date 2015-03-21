N = 8;
M = 500;
for i=1:3
    load(strcat('plain_N',int2str(N),'_M',int2str(M),'-notime.mat'));
    coeffs = abs(fftns(S.sim.gw));
    coeffs = coeffs(:);
    
    k = -N /2:N/2-1;
    [k1,k2,k3] = ndgrid(k,k,k);
    kabs = sqrt(k1.^2+k2.^2+k3.^2);
    kabs = kabs(:);
    
    clf;
    loglog(kabs,coeffs,'r*');
    hold on;
    loglog(kabs,500./kabs.^2);
    
    N = N*2;
end

