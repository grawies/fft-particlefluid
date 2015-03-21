% Set up grids and tensor product Laplacian and solve for u:
  N = 24;% [D,x] = cheb(N); y = x;
  x = (1:N)*2*pi/N;
  y = x;
  [xx,yy] = meshgrid(x,y);
  f = 10*sin(xx.*(yy-1)/(2*pi));
  klist = -N/2:1:N/2-1;
  llist = -N/2:1:N/2-1;
[k,l] = meshgrid(klist,llist);
k2 = k.^2;
l2 = l.^2;
c2 = 0;
  tic, u = iffts(ffts(f) ./ (-k2-l2 + c2)); toc          % solve problem and watch the clock
uu = u;
% Reshape long 1D results onto 2D grid:
  %uu = zeros(N+1,N+1); uu(2:N,2:N) = reshape(u,N-1,N-1);
  %[xx,yy] = meshgrid(x,y);
  value = uu(N/4+1,N/4+1);
mesh(xx,yy,uu), colormap(1e-6*[1 1 1]);
% Interpolate to finer grid and plot:
  %[xxx,yyy] = meshgrid(-1:.04:1,-1:.04:1);
  %uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
  %figure(2), clf, mesh(xxx,yyy,uuu), colormap(1e-6*[1 1 1]);
  %xlabel x, ylabel y, zlabel u
  %text(.4,-.3,-.3,sprintf('u(2^{-1/2},2^{-1/2}) = %14.11f',value))
