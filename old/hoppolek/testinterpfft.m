n = 20;
N = 160;

f1 = @(x) 0.1*(x-pi).^2;
f2 = @(x) -0.1* x .* (x-pi).^2 .* (x-2*pi);
f3 = @(x) sin(x);

x = 0:2*pi/(n-1):2*pi;
y = f1(x);

[x2,y2] = interpfft(y,N);

plot(x,y,'b');
hold on;
plot(x2,y2,'r');
grid on;