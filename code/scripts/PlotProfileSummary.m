function PlotProfileSummary(x,p,f,u,legendList)

Nplots = max(size(u));

%% plot the velocity field
figure(1)
set(gca,'YDir','reverse');
grid on;
hold all
plot(p,zeros(max(size(p)),1),'k*');
for i=1:Nplots
    plot(x,u{i});
end
set(gca,'FontSize',11);
legend(legendList)
xlabel('z-axis position [m]')
ylabel('velocity along central z-axis [m/s]')
title('velocity field profile for different regularization functions')

%% plot the regularized particle distributions
figure(2)
set(gca,'YDir','reverse'); 
grid on;
hold all
plot(p,zeros(max(size(p)),1),'k*');
for i=1:Nplots
    plot(x,f{i});
end
set(gca,'FontSize',11);
legend(legendList)
xlabel('z-axis position [m]')
ylabel('force along central z-axis [N/m^3]')
title('forcing function profiles for different regularization functions')
 
 
 