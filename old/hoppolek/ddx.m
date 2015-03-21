function yx = ddx(y)
n = size(y)
k = 1:1:max(n);
k2 = k'.*k';
yx = iffts(k2.*ffts(y));
