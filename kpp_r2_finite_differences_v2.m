%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We attempt to simulate u_t = u_xx + e^2*u_yy + f(u) using finite
% difference approximations.
%
% Solve: u_t = u_{xx} + epsilon^2*u_{yy} + u(1-u)
% Approximation:
% u(t_i + t_delta, x) = u(t_i, x) + t_delta / h^2 [ \sum u(t_i, x + h_e) - 2nu(t_i, x) ]
%                       + t_delta * u(1-u)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

% -----------------------------------------------------------------

% Initialize parameters
% We want t_delta \approx h^2 where n is the mesh size
d = 2;              % working in R^2
n = 150;             % simulation grid size
h = 1;            % mesh width
t_delta = h^2;      % time increment
time = 80;         % number of steps to simulate. time * t_delta is the total "unit time"
f = zeros(n / h);   % growth rate
epsilon = 0.1;
gamma = 0.3;        % dampen the gradient

% Initial conditions

N = n / h;          % mesh size
u = zeros(N);       % initial conditions
%v0 = 1E-2;
v0 = 0.9;

% compact set

%u(N/2, N/2) = v0; u(N/2 + 1, N/2) = v0; u(N/2, N/2 + 1) = v0; u(N/2 - 1, N/2 + 1) = v0; % nonzero on compact set

% wedge

% specify the base point of our wedge
p = [100 100];

% specify some vector
v = [1 0];
tolerance = cos(pi / 6);
for i = -N:1:N
  for j = -N:1:N
    cos_ = ([i j] * v') / (norm([i j]) * norm(v));
    if tolerance < cos_ && cos_ < 1                   % make sure our [i j] are within the wedge
      i_p = i + p(1);
      j_p = j + p(2);
      if 1 <= i_p && i_p <= N && 1 <= j_p && j_p <= N
        u(i + p(1),j + p(2)) = v0;
      end
    end
  end
end

% initialize values of u
grad = 0 * u;       % we'll store the gradient in this term
u_avg = 0 * u;

% Run simulation

% Indices using which we can represent the gradient
I = 2:N-1; J = 2:N-1;

for step=1:time
  grad(I, J) = epsilon^2*u(I, J - 1) + epsilon^2*u(I, J + 1) + u(I - 1, J) + u(I + 1, J);
  u(I, J) = u(I, J) + gamma * (t_delta / h^2) * (grad(I, J) - 2*(epsilon^2 + 1)*u(I, J)) + t_delta*(u(I,J) .* (1 - u(I,J)));
  u(1, :) = u(2, :); u(N, :) = u(N-1, :); u(:, 1) = u(:, 2); u(:, N) = u(:, N-1);

  pcolor(u); shading interp;
  colorbar; colormap hsv;
  drawnow;
end

surf(u);
shading interp; colormap jet;
view([-25 70]);
