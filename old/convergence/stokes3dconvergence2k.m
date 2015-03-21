%% PARAMETERS

% number of measurements
M = 3;
% initial resolution
N0 = 64;

%% DATA COLLECTION

N = N0;
[a,b,c] = stokes3dconvergence1fcn(N,N0,1);
k = -N /2:N/2-1;
[k1,k2,k3] = ndgrid(k,k,k);
ksq = k1.^2+k2.^2+k3.^2;
k = abs(sqrt(ksq(:)));
f = abs(c(:));
hold on;
plot(k,f,'b.');
title(strcat('fourier coefficients of forcing function, 1000 particles, N=',int2str(N))), grid on;
xlabel('wave number modulus |k|'), ylabel('coefficient modulus');
