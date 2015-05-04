LocalInit(0);

%[x,p,f,u] = Get1DProfile(2*pi, 128, 3, 'triangle',8);
[x,p,f,u] = Get1DProfile(2*pi, 128, 3, 'spline', [1 5]);

plot(p,zeros(max(size(p)),1),'k*');
hold on;
[ax,hline1,hline2] = plotyy(x,-f,x,-u);
%originLim = get(
%set(hLine1, 'Ylim', originlim);
%set(hLine1,'Color', 'r');
%set(hLine2,'Color', 'b');

xlabel('z-value','FontSize',12);

axes(ax(1)); ylabel('-f','FontSize',12);
%axis([0 2*pi 0 2]);
axes(ax(2)); ylabel('-u','FontSize',12);
%axis([0 2*pi 0 0.01]);
grid on;
title('triangle regularization','FontSize',12);
legend('solution','particles','forcing function');