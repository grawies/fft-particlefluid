Mvals = [100 200 300 500 1000 2000];

for i=1:6
    rng(13);
    M = Mvals(i);
    M
    [x,y,z,u,v,w] = stokesletnotime(M);
    T.x = x;
    T.y = y;
    T.z = z;
    T.u = u;
    T.v = v;
    T.w = w;
    filename = strcat('stokeslet-M',int2str(M),'-notime.mat');
    save(filename,'T');
end
