function Plt = SetupPlotting(S, nlines, plotinterval)

    % space
    Plt.nlines = nlines;
    Plt.I = 1:S.N/nlines:S.N;

    % timestepping
    Plt.interval = plotinterval;
    Plt.nplots = S.nmax / plotinterval;

    % init recording
    Plt.vidfig = figure(1);
    Plt.frames = S.nsteps;
    Plt.Vid(Plt.frames) = struct('cdata',[],'colormap',[]);
    Plt.framecount = 1;

    % prep plotting streamlines
    [Plt.sx,Plt.sy] = ndgrid(S.x(Plt.I),S.x(Plt.I));
    Plt.sz = S.L*ones(Plt.nlines,Plt.nlines);

    % prep plotting quiver cross section
    [Plt.cross2,Plt.cross3] = ndgrid(S.x(Plt.I),S.x(Plt.I));
    Plt.cross1 = 0.49*S.L*ones(Plt.nlines,Plt.nlines);

end