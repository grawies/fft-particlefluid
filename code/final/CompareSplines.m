LocalInit(0);

L = 2.0;
N = 256;
M = 3;

[x,p,fa,ua] = Get1DProfile(L, N, M, 'spline', [2 2]);
[~,~,fb,ub] = Get1DProfile(L, N, M, 'spline', [2 3]);
[~,~,fc,uc] = Get1DProfile(L, N, M, 'spline', [2 4]);
[~,~,fd,ud] = Get1DProfile(L, N, M, 'spline', [2 6]);

f = {fa,fb,fc,fd};
u = {ua,ub,uc,ud};

legendList = {'particle positions',...
    'B-splines, p = 2, scale = 2',...
    'B-splines, p = 3, scale = 2',...
    'B-splines, p = 4, scale = 2',...
    'B-splines, p = 6, scale = 2'};

PlotProfileSummary(x,p,f,u,legendList);
