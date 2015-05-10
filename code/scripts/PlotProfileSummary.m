LocalInit(0);
 
%[x,p,f,u] = Get1DProfile(2*pi, 128, 3, 'triangle',8);
L = 2.0;
N = 128;
M = 3;
[~,~,ftS,utS] = Get1DProfile(L, N, M, 'triangle', [4 5]);
[~,~,ftL,utL] = Get1DProfile(L, N, M, 'spline', [8 5]);
[~,~,fsSS,usSS] = Get1DProfile(L, N, M, 'spline', [1 2]);
[~,~,fsSM,usSM] = Get1DProfile(L, N, M, 'spline', [1 4]);
[x,p,fsSL,usSL] = Get1DProfile(L, N, M, 'spline', [1 6]);

%% plot the velocity field
figure(1)
set(gca,'YDir','reverse'); hold all
plot(p,zeros(max(size(p)),1),'k*');
plot(x,usSS);
plot(x,usSM);
plot(x,usSL);
plot(x,utS);
plot(x,utL);

%% plot the profiles
figure(2)
plot(p,zeros(max(size(p)),1),'k*');

xlabel('z-value','FontSize',12);

axes(ax(1)); ylabel('-f','FontSize',12);
%axis([0 2*pi 0 2]);
axes(ax(2)); ylabel('-u','FontSize',12);
%axis([0 2*pi 0 0.01]);
grid on;
title('triangle regularization','FontSize',12);
legend('solution','particles','forcing function');