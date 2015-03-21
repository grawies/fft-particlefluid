figure();
title('speed distribution');
grid on;
hold all;
Mvals = [100 200 300 500 1000];

for i = 1:5
    S.N = 32;
    S.M = Mvals(i);
    S.dt = 0.02;
    
    filename = strcat('plain_N',int2str(S.N),'_M',int2str(S.M),'_dt',num2str(S.dt),'.mat');
    load(filename);
    
    %% EVALUATION
    
    [az,el,speed] = cart2sph(S.sim.u,S.sim.v,S.sim.w);
    
    
    [f_speed,xi] = ksdensity(speed(:,2));
    plot(xi,f);
end
legend('M=100','M=200','M=300','M=500','M=1000')
%{
title(strcat('speed distribution, ',int2str(S.M),' particles'));
grid on;
print(strcat('speed-distribution-m',int2str(S.M)),'-dpng');
%}
%{
[f_el,ti] = ksdensity(el(:,2));
ti = ti * 180/pi;
f_el = f_el * pi/180;
plot(ti,f_el);
hold on;
line([-90 -90], get(gca, 'ylim'));
title(strcat('angle distribution, ',int2str(S.M),' particles'));
grid on;
%}