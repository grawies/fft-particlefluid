LocalInit(0);

L = 2.0;
N = 256;
M = 3;

[x,p,fa,ua] = Get1DProfile(L, N, M, 'spline', [12 2]);
[~,~,fb,ub] = Get1DProfile(L, N, M, 'spline', [8 3]);
[~,~,fc,uc] = Get1DProfile(L, N, M, 'spline', [6 4]);
[~,~,fd,ud] = Get1DProfile(L, N, M, 'spline', [4 6]);
[~,~,fe,ue] = Get1DProfile(L, N, M, 'spline', [3 8]);
[~,~,ff,uf] = Get1DProfile(L, N, M, 'spline', [2 12]);

f = {fa,fb,fc,fd,fe,ff};
u = {ua,ub,uc,ud,ue,uf};

legendList = {'particle positions',...
    'B-splines, p = 12, scale = 1',...
    'B-splines, p = 6, scale = 2',...
    'B-splines, p = 4, scale = 3',...
    'B-splines, p = 3, scale = 4',...
    'B-splines, p = 2, scale = 6'};

PlotProfileSummary(x,p,f,u,legendList);
