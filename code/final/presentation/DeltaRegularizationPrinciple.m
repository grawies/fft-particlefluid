
x = 0:0.001:1;

figure(1)
plot([0 0.5 0.5 0.5 1], [0 0 1 0 0], 'k')
axis([0,1,-0.1,1.1])
grid on
y = zeros(max(size(x),1));
y(round(end/2)) = 1;
hold on
plot(x(7:25:end),y(7:25:end),'k*')

figure(2)
y = exp(-512*(x-0.5).^2 );
plot(x,y,'k');
axis([0,1,-0.1,1.1])
grid on
hold on
plot(x(7:25:end),y(7:25:end),'k*')