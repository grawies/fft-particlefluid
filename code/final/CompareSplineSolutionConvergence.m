%% calculate data

legendNames = {};
plotLineSpec = {'k--','k--','k--','k--','k--','k--','k--','k--'};
plotPointSpec = {'k*','k.','k^','kv','kx','k^','kv','k<','k>'};
pvals = [2,3,4,6,8,10];
icount = 5;

N0 = 16;
iters = 12;

errvals = zeros(icount, iters-1);

for i = 1:icount
    p = pvals(i)
    dparams = [1 p];
    [Nvals,err] = ShowSplineConvergence(N0,iters,3,dparams,plotLineSpec{i},plotPointSpec{i});
    errvals(i,:) = err;
    
   
end

%% calc convergence orders

logVals = log(Nvals);
convOrder = zeros(1,icount);

for i = 1:icount
    logerr = log(errvals(i,:));
    ord = polyfit(logVals(end-5:end),logerr(end-5:end),1);
    convOrder(i) = ord(1);
end

%% plot
figure()
set(gca, 'FontSize',12);
for i = 1:icount
    loglog(Nvals,errvals(i,:), plotPointSpec{i});
    hold on;
end
for i = 1:icount
    loglog(Nvals,errvals(i,:), plotLineSpec{i});
    hold on;
end

%% fancy format

convOrderStrings = {', ~N^{-2.63}',', ~N^{-3.98}',...
    ', ~N^{-5.27}',', ~N^{-7.71}',', ~N^{-9.63}'};
for i=1:icount
    p = pvals(i);
    legendNames{i} = strcat('p = ',num2str(p), convOrderStrings{i});
end
legend(legendNames);
grid on;
xlabel('Grid resolution N','FontSize',12);
ylabel('relative 2-norm error of Stokes flows','FontSize',12);
title('convergence of solutions for different regularizations','FontSize',12);
