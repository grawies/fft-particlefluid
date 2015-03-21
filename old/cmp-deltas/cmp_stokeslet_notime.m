




load('stokeslet-M2000-notime.mat');
load('plain_N64_M2000-notime.mat');
S = S.sim;

[Taz,Tel,Tspeed] = cart2sph(T.u,T.v,T.w);
[Saz,Sel,Sspeed] = cart2sph(S.u,S.v,S.w);

[Tf,Txi] = ksdensity(Tel);
[Sf,Sxi] = ksdensity(Sel);

Sxi = Sxi * 180/pi;
Txi = Txi * 180/pi;
Sf = Sf * pi/180;
Tf = Sf * pi/180;
plot(Txi,Tf);
hold on;
plot(Sxi,Sf);
line([-90 -90], get(gca, 'ylim'));
grid on;
title('angle distribution, m = 2000');