% generate plots of the first 4 even orders of cardinal B-splines centered
% around 0

LocalInit(0);
symbol = '.ox+';


for i=1:4
    u = -i:i;
    p = 2*i;
    supp = (5-i):(4+i);
    points = CardinalSpline(u,p);
    plot(u,points,strcat('r',symbol(i)),'Markers',12);
    hold on;
end

legend('p=2','p=4','p=6','p=8');

for i=1:4
    ucont = -i:0.01:i;
    p = 2*i;
    supp = (4-i):(4+i);
    graph = CardinalSpline(ucont,p);
    hold on;
    plot(ucont,graph,'b','LineWidth',1.5);
end

for i=1:4
    u = -i:i;
    p = 2*i;
    supp = (5-i):(4+i);
    points = CardinalSpline(u,p);
    plot(u,points,strcat('r',symbol(i)),'Markers',12);
    hold on;
end

grid on;
xlabel('u','FontSize',16,'fontWeight','bold');
ylabel('Mp(u)','FontSize',16,'fontWeight','bold');
set(gca,'FontSize',12,'fontWeight','bold');


% indeed, incredibly ugly code.