function deltasum = deltasmear_tri(x,y,z,px,py,pz,epsilon)


hdeltasum = 0;
for i = 1:length(px)
    hdelta1 = max(1-abs(x-px(i))/epsilon,0);
    hdelta2 = max(1-abs(y-py(i))/epsilon,0);
    hdelta3 = max(1-abs(z-pz(i))/epsilon,0);
    hdelta = hdelta1.*hdelta2.*hdelta3;
    hdeltasum = hdeltasum + hdelta;
end

deltasum = hdeltasum;