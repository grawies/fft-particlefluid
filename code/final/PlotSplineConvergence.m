LocalInit(0);

errcount = 8;
piters = 4;
p0 = 3;
N0 = 16;

ndata = zeros(errcount);
errdata = zeros(piters,errcount);

for m = 1:piters
    p = p0 + m-1
    
    [ns,errs] = CalcSplineConvergence(p, errcount+1, N0);
    ndata = ns;
    errdata(m,:) = errs;
end

%% plot nicely

figure(1);
grid on
for m = 1:piters
    nvals = ndata(:);
    errvals = squeeze(errdata(m,:));
    errvals = errvals(:);
    
    loglog(nvals,errvals);
    hold on
    polyfit(log(nvals),log(errvals),1)
end