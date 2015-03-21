function [px,py,pz] = timestep_eulerfwd(L,dt,x,y,z,u,v,w)

px = x + dt*u;
px = px + L*(px < 0) - L*(px>L);
py = y + dt*v;
py = py + L*(py < 0) - L*(py>L);
pz = z + dt*w;
pz = pz + L*(pz < 0) - L*(pz>L);
