function PlotSliceSurfaces(S,P,Plt)

    wind_speed = sqrt(S.u.^2 + S.v.^2 + S.w.^2);
    hsurfaces = slice(S.x2,S.x1,S.x3,...
            wind_speed,...
            .5*S.L,S.L,0);
    set(hsurfaces,'FaceColor','interp','EdgeColor','none');
    colormap jet;

end