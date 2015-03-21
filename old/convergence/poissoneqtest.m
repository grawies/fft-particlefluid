
M = 8;
N0 = 2;
L = 2*pi;

ndata = zeros(M,1);
udata = zeros(M,N0,N0,N0);

N = N0;
for m = 1:M
    % space vectors
    x = (0:N-1)*L/N;
    [x,y,z] = ndgrid(x,x,x);
    f = exp(sin(x).*sin(y).*sin(z));

    % k-space vectors
    k = -N /2:N/2-1;
    [k1,k2,k3] = ndgrid(k,k,k);
    ksq = k1.^2+k2.^2+k3.^2;
    ksqinv = 1./ksq;
    ksqinv(N/2+1,N/2+1,N/2+1) = 0;
    
    I = 1:N/N0:N;
    u = -ifftns(fftns(f).*ksqinv);
    udata(m,:,:,:) = u(I,I,I);
    ndata(m) = N;
    N = N*2;
    m
end

%% ERRORS
p = 2;

errors = zeros(1,m-1);
for m=2:M
    diff = udata(m,:,:,:)-udata(m-1,:,:,:);
    errors(m-1) = norm(diff(:),p);
end

loglog(ndata(2:end),errors,'b*-');
title(strcat('p-norm errors, p = ',int2str(p))), grid on;
%fit = polyfit(log(Nvals),log(errors),1);
%legend(strcat('slope=',num2str(fit(1))));
xlabel('mesh size'), ylabel('error');