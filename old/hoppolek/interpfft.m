 function [xi, yi] = interpfft(y, N)

n = length(y);

yhat = ffts(y);
yhat = padarray(yhat',(N-n)/2,0)';
yi = real(iffts(yhat)) * N/n;


xi = 0:2*pi/(N-1):2*pi;

