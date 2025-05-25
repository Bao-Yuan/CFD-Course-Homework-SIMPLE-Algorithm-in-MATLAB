%% One-Dimensional Advection Equation Solver
% Comparison of Centered Difference, Upwind, and QUICK schemes
% Author: [Your Name]
% Date: [Date]

clear all; close all; clc;

% Call the convergence analysis function
analyze_convergence(); % This function will generate its own plot and output rates

%% Parameters
u = 1;          % Advection velocity
dx = 1;         % Spatial step
Nx = 100;       % Number of spatial points
dt = 0.1;       % Time step
Nt = 1000;      % Number of time steps
CFL = u*dt/dx;  % CFL number

% Domain
x = (0:Nx-1)*dx;
t_final = Nt*dt;

fprintf('CFL number: %.2f\n', CFL);
fprintf('Final time: %.1f\n', t_final);

%% Initial Conditions
% Case 1: Step function
f1_init = zeros(1, Nx);
f1_init(x >= 40 & x <= 60) = 1;

% Case 2: Sinusoidal function
f2_init = sin(0.02*pi*x);

%% Analytical Solutions (for comparison)
% Since it's pure advection, the solution translates
x_analytical = mod(x - u*t_final, Nx*dx);
f1_analytical = zeros(1, Nx);
f1_analytical(x_analytical >= 40 & x_analytical <= 60) = 1;
f2_analytical = sin(0.02*pi*x_analytical);

%% Solve using different schemes
fprintf('\nSolving with different numerical schemes...\n');

% Initialize solution arrays
schemes = {'Centered', 'Upwind', 'QUICK'};
colors = {'r', 'b', 'g'};

% Case 1: Step function
fprintf('\n=== Case 1: Step Function ===\n');
f1_solutions = cell(3,1);
f1_solutions{1} = centered_difference(f1_init, u, dx, dt, Nt);
f1_solutions{2} = upwind_scheme(f1_init, u, dx, dt, Nt);
f1_solutions{3} = quick_scheme(f1_init, u, dx, dt, Nt);

% Case 2: Sinusoidal function
fprintf('\n=== Case 2: Sinusoidal Function ===\n');
f2_solutions = cell(3,1);
f2_solutions{1} = centered_difference(f2_init, u, dx, dt, Nt);
f2_solutions{2} = upwind_scheme(f2_init, u, dx, dt, Nt);
f2_solutions{3} = quick_scheme(f2_init, u, dx, dt, Nt);

%% Error Analysis
fprintf('\n=== Error Analysis (for current parameters) ===\n');
fprintf('Scheme\t\tCase 1 L2 Error\tCase 2 L2 Error\n');
fprintf('------\t\t---------------\t---------------\n');

for i = 1:3
    error1 = calculate_l2_error(f1_solutions{i}, f1_analytical);
    error2 = calculate_l2_error(f2_solutions{i}, f2_analytical);
    fprintf('%s\t\t%.6f\t\t%.6f\n', schemes{i}, error1, error2);
end

%% Plotting Results
% Case 1 Results
figure('Position', [100, 100, 1200, 800]);
subplot(2,1,1);
plot(x, f1_init, 'k--', 'LineWidth', 2); hold on;
for i = 1:3
    plot(x, f1_solutions{i}, colors{i}, 'LineWidth', 1.5);
end
plot(x, f1_analytical, 'ko', 'MarkerSize', 4);
xlabel('x'); ylabel('f(x,t)');
title('Case 1: Step Function - Final Solution');
legend('Initial', schemes{:}, 'Analytical', 'Location', 'best');
grid on;

% Case 2 Results
subplot(2,1,2);
plot(x, f2_init, 'k--', 'LineWidth', 2); hold on;
for i = 1:3
    plot(x, f2_solutions{i}, colors{i}, 'LineWidth', 1.5);
end
plot(x, f2_analytical, 'ko', 'MarkerSize', 4);
xlabel('x'); ylabel('f(x,t)');
title('Case 2: Sinusoidal Function - Final Solution');
legend('Initial', schemes{:}, 'Analytical', 'Location', 'best');
grid on;

sgtitle('One-Dimensional Advection Equation: Numerical Scheme Comparison');

%% Functions

function f_new = centered_difference(f_init, u, dx, dt, Nt)
    % Centered difference scheme
    % fprintf('Running Centered Difference scheme...\n'); % Removed for cleaner output
    
    Nx = length(f_init);
    f = f_init;
    
    for n = 1:Nt
        f_old = f;
        for i = 1:Nx
            % Periodic boundary conditions
            i_plus = mod(i, Nx) + 1;
            i_minus = mod(i-2, Nx) + 1;
            
            f(i) = f_old(i) - (u*dt)/(2*dx) * (f_old(i_plus) - f_old(i_minus));
        end
    end
    
    f_new = f;
end

function f_new = upwind_scheme(f_init, u, dx, dt, Nt)
    % First-order upwind scheme
    % fprintf('Running Upwind scheme...\n'); % Removed for cleaner output
    
    Nx = length(f_init);
    f = f_init;
    
    for n = 1:Nt
        f_old = f;
        for i = 1:Nx
            % Periodic boundary conditions
            i_minus = mod(i-2, Nx) + 1;
            
            f(i) = f_old(i) - (u*dt/dx) * (f_old(i) - f_old(i_minus));
        end
    end
    
    f_new = f;
end

function f_new = quick_scheme(f_init, u, dx, dt, Nt)
    % QUICK (Quadratic Upstream Interpolation) scheme
    % fprintf('Running QUICK scheme...\n'); % Removed for cleaner output
    
    Nx = length(f_init);
    f = f_init;
    
    for n = 1:Nt
        f_old = f;
        for i = 1:Nx
            % Periodic boundary conditions
            i_plus = mod(i, Nx) + 1; % Not used in this specific QUICK formulation, but kept for consistency
            i_minus = mod(i-2, Nx) + 1;
            i_minus2 = mod(i-3, Nx) + 1;
            
            % QUICK interpolation for face value (as provided by user)
            f_wface = (6/8)*f_old(i_minus) + (3/8)*f_old(i) - (1/8)*f_old(i_minus2);
            f_eface = (6/8)*f_old(i) + (3/8)*f_old(i_plus) - (1/8)*f_old(i_minus);
            
            % Update rule using the interpolated face value (as provided by user)
            f(i) = f_old(i) - (u*dt/dx) * (-f_wface + f_eface);
        end
    end
    
    f_new = f;
end

function error = calculate_l2_error(numerical, analytical)
    % Calculate L2 norm error
    error = sqrt(mean((numerical - analytical).^2));
end

%% Additional Analysis Functions

function analyze_convergence()
    % Convergence analysis with different grid sizes
    % This function is configured to reveal the *spatial* order of accuracy
    % by making the temporal error negligible.
    
    fprintf('\n=== Spatial Convergence Analysis ===\n');
    fprintf('Note: This analysis uses a fixed final time (t_final=10).\n');
    fprintf('      Temporal step (dt) is scaled as O(dx^2) to suppress temporal error.\n');
    
    % Define the range of dx values for convergence study
    % The domain length (Lx=100) is kept constant.
    dx_values = [2, 1, 0.5, 0.25, 0.125]; % Added more points for better curve
    
    errors_centered = zeros(size(dx_values));
    errors_upwind = zeros(size(dx_values));
    errors_quick = zeros(size(dx_values));
    
    u = 1;
    t_final_conv = 10; % Fixed final time for convergence study
    
    for i = 1:length(dx_values)
        dx_curr = dx_values(i);
        
        % CRITICAL CHANGE: Scale dt with dx^2 to make temporal error negligible
        % This effectively makes CFL = u * 0.01 * dx_curr, which goes to 0 as dx_curr -> 0
        dt_curr = 0.01 * dx_curr^2; 
        
        Nx_curr = round(100/dx_curr); % Calculate Nx to keep domain length constant (Lx=100)
        Nt_curr = round(t_final_conv/dt_curr); % Calculate Nt to reach fixed final time
        
        % Ensure Nt_curr is at least 1, prevent division by zero for very small dt_curr
        if Nt_curr == 0
            Nt_curr = 1; 
            dt_curr = t_final_conv; % If t_final is very short, use it as dt
        end

        x_curr = (0:Nx_curr-1)*dx_curr;
        f_init_curr = sin(0.02*pi*x_curr); % Use sine wave for convergence analysis
        
        % Analytical solution for current dx and t_final
        x_analytical_curr = mod(x_curr - u*t_final_conv, Nx_curr*dx_curr);
        f_analytical_curr = sin(0.02*pi*x_analytical_curr);
        
        fprintf('  Running for dx = %.3f (Nx = %d, Nt = %d, dt = %.4e)...\n', dx_curr, Nx_curr, Nt_curr, dt_curr);

        % Numerical solutions
        f_centered = centered_difference(f_init_curr, u, dx_curr, dt_curr, Nt_curr);
        f_upwind = upwind_scheme(f_init_curr, u, dx_curr, dt_curr, Nt_curr);
        f_quick = quick_scheme(f_init_curr, u, dx_curr, dt_curr, Nt_curr);
        
        errors_centered(i) = calculate_l2_error(f_centered, f_analytical_curr);
        errors_upwind(i) = calculate_l2_error(f_upwind, f_analytical_curr);
        errors_quick(i) = calculate_l2_error(f_quick, f_analytical_curr);
    end
    
    % Plot convergence
    figure;
    loglog(dx_values, errors_centered, 'r-o', 'LineWidth', 2, 'DisplayName', 'Centered Difference'); hold on;
    loglog(dx_values, errors_upwind, 'b-s', 'LineWidth', 2, 'DisplayName', 'First-Order Upwind');
    loglog(dx_values, errors_quick, 'g-^', 'LineWidth', 2, 'DisplayName', 'QUICK');
    
    % Reference lines for theoretical orders
    % Scale reference lines to pass through the first data point for clear visual comparison
    loglog(dx_values, dx_values.^1 * (errors_upwind(1)/dx_values(1)^1), 'k--', 'LineWidth', 1, 'DisplayName', 'Order 1 Ref');
    loglog(dx_values, dx_values.^2 * (errors_centered(1)/dx_values(1)^2), 'k:', 'LineWidth', 1, 'DisplayName', 'Order 2 Ref');
    loglog(dx_values, dx_values.^3 * (errors_quick(1)/dx_values(1)^3), 'k-.', 'LineWidth', 1, 'DisplayName', 'Order 3 Ref'); % Added for QUICK
    
    xlabel('$\Delta x$', 'Interpreter', 'latex'); ylabel('$L_2$ Error', 'Interpreter', 'latex');
    title('Spatial Convergence Analysis for 1D Convection (Sine Wave)');
    legend('Location', 'best');
    grid on;
    
    % Calculate and display convergence rates
    fprintf('\nObserved Spatial Convergence Rates (for Sine Wave):\n');
    fprintf('%-15s %-15s\n', 'Scheme', 'Observed Rate');
    fprintf('--------------------------------\n');

    % Rates are calculated between consecutive points
    for i = 1:length(dx_values)-1
        dx1 = dx_values(i);
        dx2 = dx_values(i+1);

        rate_centered = log(errors_centered(i)/errors_centered(i+1)) / log(dx1/dx2);
        rate_upwind = log(errors_upwind(i)/errors_upwind(i+1)) / log(dx1/dx2);
        rate_quick = log(errors_quick(i)/errors_quick(i+1)) / log(dx1/dx2);
        
        fprintf('  Between dx=%.3f and dx=%.3f:\n', dx1, dx2);
        fprintf('    Centered: %.2f\n', rate_centered);
        fprintf('    Upwind:   %.2f\n', rate_upwind);
        fprintf('    QUICK:    %.2f\n', rate_quick);
    end
    fprintf('--------------------------------\n');
    fprintf('Expected Spatial Rates:\n');
    fprintf('  Centered: ~2\n');
    fprintf('  Upwind:   ~1\n');
    fprintf('  QUICK:    ~3\n');
end
