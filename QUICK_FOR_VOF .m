%% 1D Convection-Diffusion Problem with QUICK Scheme - 收敛误差分析
clc;clear all; close all;
% Parameters
L = 1.5;        % Domain length [m]
u = 2.0;        % Velocity [m/s]
rho = 1.0;      % Density [kg/m³]
Gamma = 0.03;   % Diffusion coefficient [kg/(m·s)]
a = -200;       % Source term parameter
b = 100;        % Source term parameter
x1 = 0.6;       % Source distribution parameter [m]
x2 = 0.2;       % Source distribution parameter [m]
% Time parameters
dt = 0.001;     % Time step [s]
t_final = 2.0;  % Final time [s]
nt = round(t_final / dt) + 1;  % Number of time steps

% Governing equation: ∂(ρφ)/∂t + ∂(ρuφ)/∂x = ∂/∂x(Γ∂φ/∂x) + S

%% Analytical Solution (仅用于参考对比)
fprintf('=== Analytical Solution ===\n');

% Calculate Peclet number
P = rho * u / Gamma;
fprintf('Peclet number P = %.4f\n', P);

N = 100; % Number of terms in series

% Calculate a0
a0 = ((x1 + x2) * (a * x1 + b) + b * x1) / (2 * L);
fprintf('a0 = %.6f\n', a0);

% Function to calculate an coefficients
calculate_an = @(n) (2 * L) / (n^2 * pi^2) * ...
    ((a * (x1 + x2) + b) / x2 * cos(n * pi * x1 / L) - ...
     (a + (a * x1 + b) / x2 * cos(n * pi * (x1 + x2) / L)));

% Calculate C2
C2 = a0 / (P^2 * exp(P * L));
for n = 1:N-1
    an = calculate_an(n);
    C2 = C2 + an / exp(P * L) * cos(n * pi) / (P^2 + (n * pi / L)^2);
end
fprintf('C2 = %.6f\n', C2);

% Calculate C1
C1 = -C2 + a0 / P^2;
for n = 1:N-1
    an = calculate_an(n);
    C1 = C1 + an / (P^2 + (n * pi / L)^2);
end
fprintf('C1 = %.6f\n', C1);

% Define analytical solution phi(x)
phi_analytical = @(x) -(C1 + C2 * exp(P * x) - (a0 / P^2) * (P * x + 1) - ...
    arrayfun(@(xi) sum(arrayfun(@(n) calculate_an(n) * (L / (n * pi)) * ...
    (P * sin(n * pi * xi / L) + (n * pi / L) * cos(n * pi * xi / L)) / ...
    (P^2 + (n * pi / L)^2), 1:N-1)), x)) / Gamma;

%% Numerical Solution using QUICK Scheme
fprintf('\n=== Numerical Solution (QUICK Scheme) ===\n');

% Grid generation
nx = 50;  % Number of grid points
dx = L / (nx);
x = linspace(0.5*dx, L-0.5*dx, nx);

fprintf('Grid points: %d\n', nx);
fprintf('Grid spacing dx = %.6f m\n', dx);

%% Initialize solution arrays
phi = zeros(nx, nt);  % Solution matrix: phi(x,t)
phi_old = zeros(nx, 1);  % Previous time step solution
residual_history = zeros(nt-1, 1);  % 存储收敛历史

% Source term function
S = zeros(size(x));
for i = 1:nx
    if x(i) >= 0 && x(i) <= (x1)
        S(i) = a * x(i)  + b;
    elseif x(i) >(x1) && x(i) <=(x1+x2)
        S(i) = (x(i)-(x1+x2))/x2 *(-a*x1-b);
    else
        S(i) = 0;
    end
end

% 存储每个时间步的解，用于收敛分析
phi_convergence = [];
time_steps = [];

for n = 2:nt
    
    % Current time
    t = (n-1) * dt;
    
    % Store previous time step
    phi_old = phi(:, n-1);

    % Initialize coefficient matrix and RHS vector
    A = zeros(nx, nx);
    b_vec = zeros(nx, 1);

    % Convection terms with QUICK scheme
    F = rho * u;
    % Diffusion terms
    D_e = Gamma / dx;
    D_w = Gamma / dx;

    % Interior points using QUICK scheme
    for i = 2:nx-1
        % QUICK coefficients
        if i >= 3 && i <= nx-1
            % Full QUICK scheme
        a_e = -D_e + (3/8)*F;
        a_w = -D_w -(7/8)*F;
        a_ww = F/8;
        a_p = -a_e - a_w - a_ww +rho*dx/dt;
        
        % Fill matrix
        A(i, i-2) = a_ww;
        A(i, i-1) = a_w;
        A(i, i) = a_p;
        A(i, i+1) = a_e;
        end
    end
    
    for i = 1:nx
        % Source term
        b_vec(i) = S(i) * dx + rho*dx/dt*phi_old(i);
    end
    
    % Apply boundary conditions
    % Left boundary: φ = 0
    A(1,1) = 2/3*F +rho*dx/dt+3/2*D_e;
    A(1,2) = 4/9*F - D_e;
    
    A(2,1) = -19/24*F-D_w;
    A(2,2) = 11/36*F +rho*dx/dt+2*D_e;
    A(2,3) = (3/8)*F-D_e;

    % Right boundary: ∂φ/∂x = 0
    A(nx, nx-2) = F/8;
    A(nx, nx-1) = -4/8*F -D_w;
    A(nx, nx) = 3/8*F +rho*dx/dt+D_w;

     % Solve linear system for current time step
    phi(:, n) = A \ b_vec;
    
    % Calculate residual for convergence monitoring
    residual = norm(phi(:, n) - phi_old) / (norm(phi(:, n)) + 1e-12);
    residual_history(n-1) = residual;
    
    % 每隔一定步数存储解用于收敛分析
    if mod(n-1, 100) == 0 || n == nt
        phi_convergence = [phi_convergence, phi(:, n)];
        time_steps = [time_steps, n-1];
    end
    
    % Progress reporting
    if mod(n-1, round(nt/10)) == 0 || n == nt
        fprintf('Time step %d/%d, t = %.4f s, Residual = %.2e\n', ...
                n-1, nt-1, t, residual);
    end
    
    % Check for steady state (optional early termination)
    if n > 100 && residual < 1e-8
        fprintf('Steady state reached at time step %d (t = %.4f s)\n', n-1, t);
        phi_numerical = phi(:,n);
        final_step = n-1;
        break;
    end
end

if ~exist('phi_numerical', 'var')
    phi_numerical = phi(:,end);
    final_step = nt-1;
end

fprintf('Numerical solution completed.\n');

%% 收敛误差分析
% 计算与最终收敛值的误差
convergence_errors = zeros(size(phi_convergence, 2)-1, 1);
for i = 1:size(phi_convergence, 2)-1
    convergence_errors(i) = norm(phi_convergence(:, i) - phi_numerical) / norm(phi_numerical);
end

%% Results and Visualization
% Generate fine grid for analytical solution (仅作参考)
x_fine = linspace(0, L, 200);
phi_analytical_fine = phi_analytical(x_fine);

%% 绘图
figure('Position', [100, 100, 1400, 1000]);

% 主解图 - 数值解和解析解对比
subplot(2, 3, [1, 2]);
plot(x_fine, phi_analytical_fine, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical (Reference)');
hold on;
plot(x, phi_numerical, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'QUICK Numerical (Converged)');
xlabel('Distance x [m]');
ylabel('\phi');
title('1D Convection-Diffusion Problem Solution');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 12);

% 收敛历史图
subplot(2, 3, 3);
semilogy(1:length(residual_history), residual_history, 'b-', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Residual');
title('Convergence History');
grid on;
set(gca, 'FontSize', 10);

% 收敛误差图
subplot(2, 3, 4);
if length(convergence_errors) > 1
    semilogy(time_steps(1:end-1)*dt, convergence_errors, 'r-o', 'LineWidth', 2, 'MarkerSize', 4);
    xlabel('Time [s]');
    ylabel('Relative Error to Converged Solution');
    title('Convergence Error Analysis');
    grid on;
    set(gca, 'FontSize', 10);
else
    text(0.5, 0.5, 'Insufficient data for convergence analysis', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('Convergence Error Analysis');
end

% 源项分布图
subplot(2, 3, 5);
plot(x, S, 'g-', 'LineWidth', 2);
xlabel('Distance x [m]');
ylabel('Source Term S');
title('Source Term Distribution');
grid on;
set(gca, 'FontSize', 10);

% 解的演化图
subplot(2, 3, 6);
if size(phi_convergence, 2) > 1
    plot_indices = 1:min(5, size(phi_convergence, 2));
    colors = lines(length(plot_indices));
    for i = plot_indices
        plot(x, phi_convergence(:, i), 'Color', colors(i, :), ...
             'LineWidth', 1.5, 'DisplayName', sprintf('t = %.3f s', time_steps(i)*dt));
        hold on;
    end
    xlabel('Distance x [m]');
    ylabel('\phi');
    title('Solution Evolution');
    legend('Location', 'best');
    grid on;
    set(gca, 'FontSize', 10);
else
    plot(x, phi_numerical, 'r-', 'LineWidth', 2);
    xlabel('Distance x [m]');
    ylabel('\phi');
    title('Final Solution');
    grid on;
    set(gca, 'FontSize', 10);
end

% 添加参数信息
annotation('textbox', [0.02, 0.85, 0.15, 0.12], 'String', ...
    sprintf('Parameters:\nL = %.1f m\nu = %.1f m/s\n\\rho = %.1f kg/m³\n\\Gamma = %.2f kg/(m·s)\nP = %.2f\nFinal Residual = %.2e', ...
    L, u, rho, Gamma, P, residual_history(end)), 'FontSize', 9, 'BackgroundColor', 'white');

%% Summary
fprintf('\n=== Solution Summary ===\n');
fprintf('Problem: 1D Steady Convection-Diffusion with Source Term\n');
fprintf('Numerical Method: QUICK Scheme\n');
fprintf('Domain: [0, %.1f] m\n', L);
fprintf('Grid Points: %d\n', nx);
fprintf('Peclet Number: %.4f\n', P);
fprintf('Final Time Step: %d\n', final_step);
fprintf('Final Residual: %.6e\n', residual_history(min(final_step, length(residual_history))));
fprintf('Maximum φ value: %.6f\n', max(phi_numerical));
fprintf('Minimum φ value: %.6f\n', min(phi_numerical));

if length(convergence_errors) > 0
    fprintf('Final Convergence Error: %.6e\n', convergence_errors(end));
end

fprintf('\nSolution completed successfully!\n');