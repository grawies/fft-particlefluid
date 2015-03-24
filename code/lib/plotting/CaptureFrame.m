function CaptureFrame(Plt)

    Plt.Vid(Plt.framecount) = getframe(Plt.vidfig);
    Plt.framecount = Plt.framecount + 1;

end