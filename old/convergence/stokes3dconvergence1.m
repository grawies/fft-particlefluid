%% PARAMETERS

% number of measurements
M = 3;
% initial resolution
N0 = 8;

%% DATA COLLECTION

udata = zeros(M,3*N0,N0,N0);
N = N0;
for m=1:M
    [u,v,w] = stokes3dconvergence1fcn(N,N0,0);
    %[u,v,w] = stokes3dconvergence2fcn(N,N0);
    I = 1:N/N0:N;
    u = u(I,I,I);
    v = v(I,I,I);
    w = w(I,I,I);
    
    udata(m,:,:,:) = [u;v;w];
    
    N = N*2;
    m
end

%% ERROR ANALYSIS
% entrywise matrix norm parameter
p = 2;

errors = zeros(1,M-1);
for m = 2:M
    % entrywise norm
    errors(m-1) = norm(reshape(udata(m,:,:,:)-udata(m-1,:,:,:),[],1), p);
end
Nvals = N0*2.^(1:M-1);
loglog(Nvals,errors,'b*-');
fit = polyfit(log(Nvals),log(errors),1);
title(strcat('p-norm errors, p = ',int2str(p))), grid on;
legend(strcat('slope=',num2str(fit(1))));
xlabel('mesh size'), ylabel('error');