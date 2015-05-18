function [nvals, er_rel] = CalcSplineConvergence(p, iters, N0)
    LocalInit(0);

    L = 2.0;
    M = 3;
    hmult0 = 2;

    Ndata = zeros(iters,1);
    fdata = zeros(iters,N0);
    udata = zeros(iters,N0);

    for i = 1:iters

        N = N0 * i;
        hmult = hmult0 * i;

        I = 1:i:N;

        [x,pos,f,u] = Get1DProfile(L, N, M, 'spline', [hmult, p]);

        Ndata(i) = N;
        fdata(i,:) = f(I);
        udata(i,:) = u(I);
    end

    er_rel = zeros(iters-1,1);

    for j=1:iters-1
        uj = udata(j,:);
        udiff = uj - udata(j+1,:);
        er_rel(j) = sum(udiff.^2) / sum(uj.^2);
    end

    nvals = Ndata(1:end-1);

end