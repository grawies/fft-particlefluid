function [x,y,z,u,v,w] = stokesletnotime(M)
%% PRELIMINARIES
% constants
L = 2*pi;
nu = 1.0;

% particle positions
R0 = 1.0;
phi = unifrnd(0,2*pi, 1, M);
costheta = unifrnd(-1,1, 1, M);
u = unifrnd(0,1,1,M);
theta = acos( costheta );
r = R0 * u.^(1/3);

pointx = .5 * L + r .* sin( theta ) .* cos( phi );
pointy = .5 * L + r .* sin( theta ) .* sin( phi );
pointz = .5 * L + r .* cos( theta );

pointu = zeros(1,M);
pointv = zeros(1,M);
pointw = zeros(1,M);

%% FUNTIONS, EQUATION PROPERTIES

% gravity force on point particles
g = 9.82;
mrel = 0.1;
frel = mrel * g;
fgx = frel * 0;
fgy = frel * 0;
fgz = frel * -1;

%% CALCULATION

distancex = zeros(M,M);
distancey = zeros(M,M);
distancez = zeros(M,M);
distancei = zeros(M,M);
distancei3 = zeros(M,M);
    % calculate particle distances
    tic;
    for i = 1:M
        for j = 1:M
            Dx = pointx(i)-pointx(j);
            distancex(i,j) = Dx;
            distancex(j,i) = -Dx;
            Dy = pointy(i)-pointy(j);
            distancey(i,j) = Dy;
            distancey(j,i) = -Dy;
            Dz = pointz(i)-pointz(j);
            distancez(i,j) = Dz;
            distancez(j,i) = -Dz;
            if i==j
                D = Inf;
            else
                D = sqrt(Dx^2+Dy^2+Dz^2);
            end
            distancei(i,j) = 1/D;
            distancei(j,i) = 1/D;
            D3 = D^3;
            distancei3(i,j) = 1/D3;
            distancei3(j,i) = 1/D3;
        end
    end
    % calculate particle velocity
    for i = 1:M
        dx = distancex(i,:);
        dy = distancey(i,:);
        dz = distancez(i,:);
        di = distancei(i,:);
        di3 = distancei3(i,:);
        dot = dx.*fgx+dy.*fgy+dz.*fgz;
        dotdiv = dot.*di3;
        pointu(i) = sum(fgx.*di + dx .* dotdiv);
        pointv(i) = sum(fgy.*di + dy .* dotdiv);
        pointw(i) = sum(fgz.*di + dz .* dotdiv);
        
    end

x = pointx;
y = pointy;
z = pointz;
u = pointu;
v = pointv;
w = pointw;
