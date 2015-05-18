function PlotAll(S,prevpos,P,Plt,viewArgs,toPlot)

    clf;
    axis([0 S.L 0 S.L 0 S.L]);
    grid on;
    set(gca,'FontSize',11);
    title(strcat('t = ',num2str(S.t)));
    xlabel('x-position [m]');
    ylabel('y-position [m]');
    zlabel('z-position [m]');
    
    view(viewArgs(1),viewArgs(2));
    hold all;

    if toPlot(1)
        PlotStreamLines(S,P,Plt);
    end
    
    if toPlot(2)
        PlotParticles(S,prevpos,P,Plt);
    end
    
    if toPlot(3)
        PlotSliceSurfaces(S,P,Plt);
    end
    
    if toPlot(4)
        PlotFluidQuiver(S,P,Plt);
    end
    
    if toPlot(5)
        CaptureFrame(Plt);
    end
    
    drawnow;


end