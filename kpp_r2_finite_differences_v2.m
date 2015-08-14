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
n = 200;            % simulation grid size
t_delta = 1;     % time increment
h = t_delta^2;            % mesh width
time = 1E2;         % number of steps to simulate. time * t_delta is the total "unit time"
f = zeros(n / h);   % growth rate
epsilon = 1E-1;

% Initial conditions

N = n / h;          % mesh size
u = zeros(N);       % initial conditions
u(n/2, n/2) = 1E-2; % nonzero on compact set
grad = 0 * u;       % we'll store the gradient in this term

% Run simulation

% Indices using which we can represent the gradient
I = 2:N-1; J = 2:N-1;

for step=1:time
  f = 0;
% f = u .* (ones(N) - u);
  grad(I, J) = u(I, J - 1) + u(I, J + 1) + u(I - 1, J) + u(I + 1, J);
  grad(1, :) = grad(2, :); grad(N, :) = grad(N-1, :); grad(:, 1) = grad(:, 2); grad(:, N) = grad(:, N-1);
  u = u + (t_delta / h^2) * (grad - 2*d*u) + t_delta*f;
  pcolor(u); shading interp;
  colorbar; colormap hsv;
  drawnow;
end

surf(u);
shading interp; colormap jet;
view([-25 70]);
