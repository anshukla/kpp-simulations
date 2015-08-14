% Solve: u_t=D*(u_{xx}+u_{yy})+gamma*q(u) where q(u)='u.*(1-u)';
%
% The solution of this PDE is obtained by the finite difference
% method, assuming dx=dy=dt=1 as desribed in _An Introduction to Computational
% Engineering With Matlab by Xin-She Yang.

% ----------------------------------------------------
function kpp(n)
% Default the domain size n by n
if nargin<1, n=200; end

% Initialize parameters
% ---- time=100, D=0.2; gamma=0.5; -------------------

time=100; D=0.2; gamma=0.5;

% ---- Initial values for u --------------
% Non-zero on a compact set
%  u=zeros(n);  grad=u*0;
%  u(n/2, n/2) = 1E-5;

% Random
% u=rand(n);   grad=u*0;

% Uniformly 1
% u=ones(n);     grad=u*0;

% Nonzero on a strip
% u = [zeros(n, 90) 0.9*ones(n, 20) zeros(n, 90)];
% grad = u*0;
u=[ones(n, n/4), zeros(n, 3*n/4)];
grad = u*0;

% Vectorization/index for u(i,j) and the loop --------
I = 2:n-1; J = 2:n-1;

for step=1:time,
% Laplace gradient of the equation
 grad(I,J)= u(I,J-1)+u(I,J+1)+u(I-1,J)+u(I+1,J);
 u =(1-4*D)*u+D*grad+gamma*u.*(1-u);
 pcolor(u);  shading interp;
 colorbar; colormap hsv;
 drawnow;
end

surf(u);
shading interp; colormap jet;
view([-25 70]);
