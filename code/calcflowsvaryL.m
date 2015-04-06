%% DATA COLLECTION

n = 6;
n0 = 3;
udata = zeros(2^n,2^n,2^n,n-n0+1);
vdata = zeros(2^n,2^n,2^n,n-n0+1);
wdata = zeros(2^n,2^n,2^n,n-n0+1);
k = 1;
for i = n0:n
    k = k*4;
    i
    r = (2^n - 2^i)/2;
    pad = [r,r,r];
    [u,v,w] = CalcInitialFlow(2^i*pi, 2^i, true);
    up = k*padarray(u,pad,0);
    vp = k*padarray(v,pad,0);
    wp = k*padarray(w,pad,0);
    udata(:,:,:,i+1-n0) = up;
    vdata(:,:,:,i+1-n0) = vp;
    wdata(:,:,:,i+1-n0) = wp;
end

%% ANALYSIS 

f = [];
for i = 2:n+1-n0
    d1 = sort(udata(:,:,:,i))-sort(udata(:,:,:,i-1));
    d2 = sort(vdata(:,:,:,i))-sort(vdata(:,:,:,i-1));
    d3 = sort(wdata(:,:,:,i))-sort(wdata(:,:,:,i-1));
    q = abs(d1.*d1+d2.*d2+d3.*d3);
    f = [f, sum(q(:))];
end

f