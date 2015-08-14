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
n = 20;             % simulation grid size
h = 0.5;            % mesh width
t_delta = h^2;      % time increment
time = 1E2;         % number of steps to simulate. time * t_delta is the total "unit time"
f = zeros(n / h);   % growth rate
epsilon = 1;
gamma = 0.1;        % dampen the gradient

% Initial conditions

N = n / h;          % mesh size
u = zeros(N);       % initial conditions
v0 = 1E-2;
u(N/2, N/2) = v0; u(N/2 + 1, N/2) = v0; u(N/2, N/2 + 1) = v0; u(N/2 + 1, N/2 + 1) = v0; % nonzero on compact set
grad = 0 * u;       % we'll store the gradient in this term
u_avg = 0 * u;

% Run simulation

% Indices using which we can represent the gradient
I = 2:N-1; J = 2:N-1;

for step=1:time
  %f = 0;
  u_avg(I, J) = (1/5)*(u(I, J) + u(I, J-1) + u(I, J+1) + u(I-1, J) + u(I + 1, J));
  u_avg(1, :) = u_avg(2, :); u_avg(N, :) = u_avg(N-1, :); u_avg(:, 1) = u_avg(:, 2); u_avg(:, N) = u_avg(:, N-1);
  f = u_avg .* (1 - u_avg);
  grad(I, J) = 1/gamma*(epsilon^2*u(I, J - 1) + epsilon^2*u(I, J + 1) + u(I - 1, J) + u(I + 1, J));
  grad(1, :) = grad(2, :); grad(N, :) = grad(N-1, :); grad(:, 1) = grad(:, 2); grad(:, N) = grad(:, N-1);
  u = u + (t_delta / h^2) * (grad - 2*(epsilon^2 + 1)*u) + t_delta*f;
  if min(min(u)) < 0
      fprintf('something just became negative');
  end
  if max(max(isnan(u))) == 1
      fprintf('something messed up here')
  end
  %set(gcf,'renderer','painters');
  pcolor(u); %shading interp;
  colorbar; colormap hsv;
  drawnow;
end

surf(u);
shading interp; colormap jet;
view([-25 70]);
