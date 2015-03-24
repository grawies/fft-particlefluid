function PlotAll(S,P,Plt,viewArgs,toPlot)

    clf;
    axis([0 S.L 0 S.L 0 S.L]);
    view(viewArgs(1),viewArgs(2));
    hold all;

    if toPlot(1)
        PlotStreamLines(S,P,Plt);
    end
    
    if toPlot(2)
        PlotParticles(S,P,Plt);
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