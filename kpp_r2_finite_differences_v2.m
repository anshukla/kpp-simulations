%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We attempt to simulate u_t = u_xx + e^2*u_yy + mu(x)f(u) using finite
% difference approximations.
%
% Solve: u_t = u_{xx} + epsilon^2*u_{yy} + mu(x)f(u)
% Approximation:
% u(t_i + t_delta, x) = u(t_i, x) + t_delta / h^2 [ \sum u(t_i, x + h_e) - 2nu(t_i, x) ]
%                       + t_delta * u(1-u)
%
% NOTES:
%       1. To animate, you may have to increase your JVM heap size
%          by going to Preferences > General > Java Heap Memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

% -----------------------------------------------------------------

COMPACT_SUPPORT = 0;
WEDGE           = 1;

MU_UNIFORM      = 0;
MU_PERIODIC     = 1;
MU_CUSTOM       = 2;

INITIAL_DATA_TYPE = COMPACT_SUPPORT;        % set to the type of initial condition
                                            % we want to simulate
MU_TYPE = MU_UNIFORM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulation size parameters
d = 2;                  % working in R^2
n = 150;                % simulation grid size
h = 1;                  % mesh width

% We want t_delta \approx h^2 where n is the mesh size. This cancels out
% some ugly constant factors in the update equation.
t_delta = h^2;          % time increment
time = 120;              % number of steps to simulate. time * t_delta is the total "unit time"
f = zeros(n / h);       % growth rate
epsilon = 1;
gamma = 0.3;            % dampen the gradient

% Parameters for compactly supported initial condition.
v0 = 9E-1;              % some small nonzero intial value

% Parameters for wedge initial condition direction.
p = [100 100];                % specify the base point of our wedge
v = [1 1];                    % specify the vector which determines wedge direction
theta = 2*pi / 5;              % angle within which you want to set up the wedge

                              % expression for ``f`` in terms of u(I, J)
f_expr = '(mu(I, J) .* (u(I,J) .* (1 - u(I,J)) .* (u(I,J) - 1/4)))';
%f_expr = 'mu(I, J) .* (u(I, J) .* (1 - u(I, J)))';

mu_expr = '';                 % expression for the (i, j)th entry of mu
if MU_TYPE == MU_UNIFORM      % offer some sensible defaults for mu
  mu_expr = '1';
elseif MU_TYPE == MU_PERIODIC
  mu_expr = '1 + sin ( (pi / 4) * (i + j))';
elseif MU_TYPE == MU_CUSTOM
  mu_expr = 'TYPE_CUSTOM_EXPRESSION_HERE'
end

% Animation parameters
output_folder = 'output/';    % folder where animations are placed
rde_type      = 'bsn';
                              % type of RDE we are modeling, determiend by mu and f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = n / h;              % mesh size
u = zeros(N);           % the matrix ``u`` will hold our simulated values for the solution
grad = 0 * u;           % we'll store the neighbor_average in this term
                        % TODO(ansh): rename to ``neighbors``
filename = '';

% Nonzero on a centered square.
if INITIAL_DATA_TYPE == COMPACT_SUPPORT
  filename = strcat(output_folder, sprintf('%s_compact_support_epsilon_%0.2f_v0_%.2f.gif',...
                                            rde_type, epsilon, v0));
  u(N/2, N/2) = v0; u(N/2 + 1, N/2) = v0; u(N/2, N/2 + 1) = v0; u(N/2 - 1, N/2 + 1) = v0;
end

if INITIAL_DATA_TYPE == WEDGE
  filename = strcat(output_folder, sprintf('%s_epsilon_%.2f_v_%d_%d_theta_%.2f_p_%d_%d.gif',...
                                            rde_type, v(1), v(2), theta, p(1), p(2)));
  tolerance = cos(theta);       % sepcify the degree around the vector which we want nonzero
  for i = -N:1:N
    for j = -N:1:N
      cos_ = ([i j] * v') / (norm([i j]) * norm(v));
      if tolerance <= cos_ && cos_ <= 1                   % make sure our [i j] are within the wedge
        i_p = i + p(1);
        j_p = j + p(2);
        if 1 <= i_p && i_p <= N && 1 <= j_p && j_p <= N
          u(i + p(1),j + p(2)) = v0;
        end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Indices using which we can represent the "inner" terms of the matrix
I = 2:N-1; J = 2:N-1;

% TODO(ansh): add comment about how to control mu
mu = zeros(N);
for i = 1:N
  for j = 1:N
    mu(i, j) = eval(mu_expr);
  end
end

for step=1:time
  grad(I, J) = u(I, J - 1) + u(I, J + 1) + epsilon^2*u(I - 1, J) + epsilon^2*u(I + 1, J);
  u(I, J) = u(I, J) + gamma * (t_delta / h^2) * (grad(I, J) - 2*(epsilon^2 + 1)*u(I, J)) + t_delta*eval(f_expr);
  u(1, :) = u(2, :); u(N, :) = u(N-1, :); u(:, 1) = u(:, 2); u(:, N) = u(:, N-1);

  % normalize so that no value in the solution goes beyond 1
  for i = 1:N
    for j = 1:N
      if u(i, j) > 1
        u(i, j) = 1;
      end
    end
  end

  pcolor(u); shading interp;
  colorbar; colormap jet;
  drawnow;

  % Append the current frame to the animation
  frame = getframe(1);
  im = frame2im(frame);
  [A, map] = rgb2ind(im, 256);
  if step == 1
    imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
  else
    imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
  end
end

surf(u);
shading interp; colormap jet;
view([-25 70]);
