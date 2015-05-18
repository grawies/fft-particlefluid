LocalInit(0);

L = 2.0;
N = 256;
M = 3;

[x,p,fa,ua] = Get1DProfile(L, N, M, 'triangle', [2 0]);
[~,~,fb,ub] = Get1DProfile(L, N, M, 'triangle', [3 0]);
[~,~,fc,uc] = Get1DProfile(L, N, M, 'triangle', [4 0]);
[~,~,fd,ud] = Get1DProfile(L, N, M, 'triangle', [6 0]);

f = {fa,fb,fc,fd};
u = {ua,ub,uc,ud};

legendList = {'particle positions',...
    'triangles, \epsilon = 2h',...
    'triangles, \epsilon = 3h',...
    'triangles, \epsilon = 4h',...
    'triangles, \epsilon = 6h'};

PlotProfileSummary(x,p,f,u,legendList);
