function main(N,M,tstep)
%% SETUP

S.N = N;
S.M = M;
S.dt = 0.02;
S.seed = 13;
S.nsteps = 200;
S.L = 2*pi;
S.nu = 1;
S.radius = 1.0;
S.plotnum = S.nsteps;
S.timestep = tstep;
S.animate = true;
S.epsilon = S.L/16;
%% SIMULATION

[x,y,z,u,v,w,gx,gy,gz,gu,gv,gw] = getpvel(S);

S.sim.x = x;
S.sim.y = y;
S.sim.z = z;
S.sim.u = u;
S.sim.v = v;
S.sim.w = w;
S.sim.gx = gx;
S.sim.gy = gy;
S.sim.gz = gz;
S.sim.gu = gu;
S.sim.gv = gv;
S.sim.gw = gw;

%% SAVE
if S.timestep
filename = strcat('plain_N',int2str(S.N),'_M',int2str(S.M),'_dt',num2str(S.dt),'.mat');
else
    filename = strcat('plain_N',int2str(S.N),'_M',int2str(S.M),'-notime.mat');
end
save(filename,'S');

