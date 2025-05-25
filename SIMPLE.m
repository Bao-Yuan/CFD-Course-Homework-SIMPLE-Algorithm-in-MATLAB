%% 2D Lid-Driven Cavity Flow Simulation using SIMPLE Algorithm
% This code simulates the 2D lid-driven cavity flow problem
% using the SIMPLE algorithm on a staggered grid.
close all;
clear all;
clc;
LD_ratio = 1.25;
Re = 2300; % Reynolds number

% Start total timer
tic;


%% Parameters
% Physical parameters
L  = 1.0;             % Length of the cavity
D = L / LD_ratio;     % Depth of the cavity
rho = 1;              % Density
U0 = 1 ;              % Lid velocity (top wall)
mu =  U0*L*rho/Re;    % Dynamic viscosity (calculated from Re)


% Numerical parameters
nx = 128;        % Number of internal grid points (control volumes) in x direction
ny = 128;        % Number of internal grid points (control volumes) in y direction
maxIter = 10000;       % Maximum number of iterations
continuity_tolerance = 7e-9;  % Convergence tolerance for continuity
alpha_p = 0.5;         % Under-relaxation factor for pressure

% Print simulation info
fprintf('Starting simulation for cavity flow:\n');
fprintf('Reynolds number: %d\n', Re);
fprintf('L/D ratio: %.2f\n', LD_ratio);
fprintf('Grid size (internal cells): %d x %d\n', nx, ny);
fprintf('Using continuity-based convergence with tolerance: %.2e\n', continuity_tolerance);

%% Grid setup (Adhering to Section 2.1 of vanilaCavityFlow_EN.md)
dx = L / nx;     % Grid spacing in x direction
dy = D / ny;     % Grid spacing in y direction

% Grid points for plotting (cell centers of internal domain)
x_plot = linspace(dx/2, L - dx/2, nx); % Cell centers for plotting
y_plot = linspace(dy/2, D - dy/2, ny); % Cell centers for plotting
[X, Y] = meshgrid(x_plot, y_plot); % Grid matrices for plotting

u = zeros(ny + 2, nx + 1);   % x-velocity
v = zeros(ny + 1, nx + 2);   % y-velocity
p = zeros(ny, nx);           % pressure (only internal cells)
pc = zeros(ny, nx);          % pressure correction (only internal cells)

% Initialize predicted velocities
u_star = zeros(ny + 2, nx + 1);
v_star = zeros(ny + 1, nx + 2);

% Coefficients for momentum equations
ae_u = zeros(ny + 2, nx + 1);
aw_u = zeros(ny + 2, nx + 1);
an_u = zeros(ny + 2, nx + 1);
as_u = zeros(ny + 2, nx + 1);
ap_u = zeros(ny + 2, nx + 1);

ae_v = zeros(ny + 1, nx + 2);
aw_v = zeros(ny + 1, nx + 2);
an_v = zeros(ny + 1, nx + 2);
as_v = zeros(ny + 1, nx + 2);
ap_v = zeros(ny + 1, nx + 2);

% Coefficients for pressure correction equation
ae_p = zeros(ny, nx);
aw_p = zeros(ny, nx);
an_p = zeros(ny, nx);
as_p = zeros(ny, nx);
ap_p = zeros(ny, nx);

% Pressure source term (mass imbalance)
b = zeros(ny, nx);
% --- MODIFIED END ---

% Arrays for storing results
continuity_residual_hist = zeros(maxIter, 1);
u_residual_hist = zeros(maxIter, 1);
v_residual_hist = zeros(maxIter, 1);

% Initialize maximum continuity residual for first iteration
initial_continuity_residual = 1.0;

% Initialize the velocity field
u(ny + 2, :) = 2*U0 - u(ny + 1, :);

%% Main SIMPLE Algorithm
converged = false;
iteration = 0;

while ~converged && iteration < maxIter
    iteration = iteration + 1;
    
    % Store old velocities for convergence check
    u_old = u;
    v_old = v;
    
    
    % u-velocity boundary conditions
    % Bottom wall (y=0): No-slip, u=0. Ghost cell u(1,j) mirrors u(2,j)
    u(1, :) = -u(2, :);
    % Top wall (y=D): Moving lid, u=U0. Ghost cell u(ny+2,j) (above y=D)
    u(ny + 2, :) = 2*U0 - u(ny + 1, :);
    % Left wall (x=0): No-penetration, u=0. u(i,1) is on the boundary
    u(:, 1) = 0;
    % Right wall (x=L): No-penetration, u=0. u(i,nx+1) is on the boundary
    u(:, nx + 1) = 0;
    
    % v-velocity boundary conditions
    % Bottom wall (y=0): No-penetration, v=0. v(1,j) is on the boundary
    v(1, :) = 0;
    % Top wall (y=D): No-penetration, v=0. v(ny+1,j) is on the boundary
    v(ny + 1, :) = 0;
    % Left wall (x=0): No-slip, v=0. Ghost cell v(i,1) mirrors v(i,2)
    v(:, 1) = -v(:, 2);
    % Right wall (x=L): No-slip, v=0. Ghost cell v(i,nx+2) mirrors v(i,nx+1)
    v(:, nx + 2) = -v(:, nx + 1);

    %% Step 1: Solve momentum equations
    % Coefficients for u-momentum equation (for INTERIOR u-faces)
    for i = 2:ny + 1 % Iterate over y-indices of u-points (from u_y1 to u_y(ny))
        for j = 2:nx % Iterate over x-indices of u-points (from u_x1 to u_x(nx-1))
            % Convection terms
            % Fe is flux through east face of u-CV. Interpolate u at the east face.
            Fe = rho * dy * 0.5*(u_old(i,j+1)+u_old(i,j));
            % Fw is flux through west face of u-CV. The velocity at this face is u_old(i,j).
            Fw = rho * dy * 0.5*(u_old(i,j)+u_old(i,j-1));
            % Fn is flux through north face of u-CV. Interpolate v at the north face.
            Fn = rho * dx * 0.5 * (v_old(i,j) + v_old(i,j+1));
            % Fs is flux through south face of u-CV. Interpolate v at the south face.
            Fs = rho * dx * 0.5 * (v_old(i-1,j) + v_old(i-1,j+1));
            
            % Diffusion terms
            De = mu * dy / dx; Dw = mu * dy / dx;
            Dn = mu * dx / dy; Ds = mu * dx / dy;
            
            % Coefficients (using Upwind Differencing Scheme for convection)
            aw_u(i,j) = Dw + max(0, Fw);
            ae_u(i,j) = De + max(0, -Fe);
            as_u(i,j) = Ds + max(0, Fs);
            an_u(i,j) = Dn + max(0, -Fn);
            
            % Calculate ap coefficient
            ap_u(i,j) = aw_u(i,j) + ae_u(i,j) + as_u(i,j) + an_u(i,j); 
            
            % Source term (pressure gradient)
            Su = -(p(i-1,j) - p(i-1,j-1)) * dy;
            
            % Solve for u_star
            u_star(i,j) = (aw_u(i,j)*u(i,j-1) + ae_u(i,j)*u(i,j+1) + as_u(i,j)*u(i-1,j) + an_u(i,j)*u(i+1,j) + Su) / ap_u(i,j);
        end
    end
    
    
    % Coefficients for v-momentum equation (for INTERIOR v-faces)
    for i = 2:ny % Iterate over y-indices of v-points (from v_y1 to v_y(ny-1))
        for j = 2:nx + 1 % Iterate over x-indices of v-points (from v_x1 to v_x(nx))
            % Convection terms
            % Fe is flux through east face of v-CV. Interpolate u at the east face.
            Fe = rho * dy * 0.5 * (u_old(i,j) + u_old(i+1,j));
            % Fw is flux through west face of v-CV. Interpolate u at the west face.
            Fw = rho * dy * 0.5 * (u_old(i,j-1) + u_old(i+1,j-1));
            % Fn is flux through north face of v-CV.Interpolate v at the north face.
            Fn = rho * dx * (v_old(i+1,j)+v_old(i,j));
            % Fs is flux through south face of v-CV. The velocity at this face is v_old(i,j).
            Fs = rho * dx * (v_old(i,j)+v_old(i-1,j));
            
            % Diffusion terms
            De = mu * dy / dx; Dw = mu * dy / dx;
            Dn = mu * dx / dy; Ds = mu * dx / dy;
            
            % Coefficients (using Upwind Differencing Scheme for convection)
            aw_v(i,j) = Dw + max(0, Fw);
            ae_v(i,j) = De + max(0, -Fe);
            as_v(i,j) = Ds + max(0, Fs);
            an_v(i,j) = Dn + max(0, -Fn);
            
            % Calculate ap coefficient
            ap_v(i,j) = aw_v(i,j) + ae_v(i,j) + as_v(i,j) + an_v(i,j); 
            
            % Source term (pressure gradient)
            Sv = -(p(i,j-1) - p(i-1,j-1)) * dx;
            
            % Solve for v_star
            v_star(i,j) = (aw_v(i,j)*v(i,j-1) + ae_v(i,j)*v(i,j+1) + as_v(i-1,j)*v(i-1,j) + an_v(i,j)*v(i+1,j) + Sv) / ap_v(i,j);
        end
    end
    
    %% Step 2: Calculate mass imbalance and solve pressure correction equation
    % Calculate source term b (mass imbalance)
    for i = 1:ny % y-index for p-cell
        for j = 1:nx % x-index for p-cell
            b(i,j) = -rho * (u_star(i+1,j+1) * dy - u_star(i+1,j) * dy + v_star(i+1,j+1) * dx - v_star(i,j+1) * dx);
        end
    end

    % Track maximum and total mass imbalance for convergence check
    max_mass_imbalance = max(abs(b(:)));
    total_mass_imbalance = sum(abs(b(:)));
    
    % Normalize the total mass imbalance
    num_interior_cells = ny * nx;
    avg_mass_imbalance = total_mass_imbalance / num_interior_cells;
    
    % Store the initial continuity residual on first iteration
    if iteration == 1
        initial_continuity_residual = avg_mass_imbalance;
        if initial_continuity_residual < 1e-8
            initial_continuity_residual = 1.0;
        end
    end
    
    % Set coefficients for pressure correction equation (for internal p-cells) (Equation 2.13 in MD)
    for i = 1:ny % y-index for p-cell
        for j = 1:nx % x-index for p-cell
            if i == 1 && j == 1
                as_p(i,j) = 0;
                aw_p(i,j) = 0;
                ae_p(i,j) = rho * dy^2 / ap_u(i+1,j+1);
                an_p(i,j) = rho * dx^2 / ap_v(i+1,j+1);   
            elseif i == 1 && j== nx
                as_p(i,j) = 0;
                ae_p(i,j) = 0;
                aw_p(i,j) = rho * dy^2 / ap_u(i+1,j);
                an_p(i,j) = rho * dx^2 / ap_v(i+1,j+1);    
            elseif i == ny && j==1
                as_p(i,j) = rho * dx^2 / ap_v(i,j+1);
                an_p(i,j) = 0;
                aw_p(i,j) = 0;
                ae_p(i,j) = rho * dy^2 / ap_u(i+1,j+1);  
             elseif i == nx && j== nx
                an_p(i,j) = 0;
                ae_p(i,j) = 0;
                aw_p(i,j) = rho * dy^2 / ap_u(i+1,j);
                as_p(i,j) = rho * dx^2 / ap_v(i,j+1);  
            elseif i == 1 
                as_p(i,j) = 0;
                ae_p(i,j) = rho * dy^2 / ap_u(i+1,j+1);
                aw_p(i,j) = rho * dy^2 / ap_u(i+1,j);
                an_p(i,j) = rho * dx^2 / ap_v(i+1,j+1);
            elseif i == ny
                an_p(i,j) = 0;
                ae_p(i,j) = rho * dy^2 / ap_u(i+1,j+1);
                aw_p(i,j) = rho * dy^2 / ap_u(i+1,j);
                as_p(i,j) = rho * dx^2 / ap_v(i,j+1);
            elseif j == 1
                aw_p(i,j) = 0;
                ae_p(i,j) = rho * dy^2 / ap_u(i+1,j+1);
                an_p(i,j) = rho * dx^2 / ap_v(i+1,j+1);
                as_p(i,j) = rho * dx^2 / ap_v(i,j+1);
            elseif j== nx
                ae_p(i,j) = 0;
                aw_p(i,j) = rho * dy^2 / ap_u(i+1,j);
                an_p(i,j) = rho * dx^2 / ap_v(i+1,j+1);
                as_p(i,j) = rho * dx^2 / ap_v(i,j+1);
            else   
            ae_p(i,j) = rho * dy^2 / ap_u(i+1,j+1);
            aw_p(i,j) = rho * dy^2 / ap_u(i+1,j);
            an_p(i,j) = rho * dx^2 / ap_v(i+1,j+1);
            as_p(i,j) = rho * dx^2 / ap_v(i,j+1);
            end
            ap_p(i,j) = ae_p(i,j) + aw_p(i,j) + an_p(i,j) + as_p(i,j);
        end
    end
    
    % Solve pressure correction equation using Gauss-Seidel
    pc = ones(ny, nx); % Re-initialize for this iteration
    max_pc_iter = 100;
    pc_tolerance = 1e-8;
    pc_converged = false;

    for pc_iter = 1:max_pc_iter
        pc_old = pc;
        
        for i = 1:ny % y-index for p-cell
            for j = 1:nx % x-index for p-cell
                pc_E = pc(i, min(j+1, nx));
                pc_W = pc(i, max(j-1, 1));
                pc_N = pc(min(i+1, ny), j);
                pc_S = pc(max(i-1, 1), j);
                pc(i,j) = (ae_p(i,j) * pc_E + aw_p(i,j) * pc_W + an_p(i,j) * pc_N + as_p(i,j) * pc_S + b(i,j)) / ap_p(i,j);
            end
        end

       
        
        % Check convergence of pressure correction equation
        pc_residual = max(max(abs(pc - pc_old))) / (max(max(abs(pc))) + 1e-10);
        if pc_residual < pc_tolerance
            pc_converged = true;
            break;
        end
    end
    
    
    %% Step 3: Correct pressure and velocities
    % Correct pressure with under-relaxation (for internal p-cells)
    for i = 1:ny % y-index for p-cell
        for j = 1:nx % x-index for p-cell
            p(i,j) = p(i,j)+alpha_p* pc(i,j);
        end
    end

    % Correct velocities (for INTERIOR faces based on MD's u/v locations)
    % Correct u-velocity (for u(i,j) where i is 2:ny+1 and j is 2:nx)
    for i = 2:ny + 1 % y-index for u
        for j = 2:nx % x-index for u
            u(i,j) = u_star(i,j) + dy * (pc(i-1,j-1) - pc(i-1,j)) / ap_u(i,j);
        end
    end
    
    % Correct v-velocity (for v(i,j) where i is 2:ny and j is 2:nx+1)
    for i = 2:ny % y-index for v
        for j = 2:nx + 1 % x-index for v
            v(i,j) = v_star(i,j) + dx * (pc(i-1,j-1) - pc(i,j-1)) / ap_v(i,j);
        end
    end
    
    %% Check convergence based on continuity equation
    % Calculate momentum residuals (based on change in velocity, similar to MathWorks' relative residual)
    % Sum over all elements (including ghost cells and boundary faces)
    res_u = sum(abs(u(:) - u_old(:))) / sum(abs(u_old(:)) + 1e-10);
    res_v = sum(abs(v(:) - v_old(:))) / sum(abs(v_old(:)) + 1e-10);
    
    % Store residuals for history
    u_residual_hist(iteration) = res_u;
    v_residual_hist(iteration) = res_v;
    continuity_residual_hist(iteration) = avg_mass_imbalance;
    
    % Calculate relative continuity residual
    relative_continuity_residual = avg_mass_imbalance / initial_continuity_residual;
    
    % Print iteration info
    if mod(iteration, 10) == 0 || iteration == 1
        fprintf('Iteration %d: Continuity Residual = %.6e (Relative: %.6e), U Residual = %.6e, V Residual = %.6e\n', iteration, avg_mass_imbalance, relative_continuity_residual, res_u, res_v);
    end
    
    % Check convergence based on continuity equation
    if avg_mass_imbalance < continuity_tolerance
        converged = true;
        fprintf('\nConverged after %d iterations!\n', iteration);
        fprintf('Final continuity residual: %.6e\n', avg_mass_imbalance);
        fprintf('Final u-momentum residual: %.6e\n', res_u);
        fprintf('Final v-momentum residual: %.6e\n', res_v);
    end
    
    % Check for stalled convergence
    if iteration > 200 && mod(iteration, 100) == 0
        last_100_min = min(continuity_residual_hist(max(1, iteration-99):iteration));
        last_100_max = max(continuity_residual_hist(max(1, iteration-99):iteration));
        
        if (last_100_max - last_100_min) / last_100_max < 0.005
            fprintf('\nWarning: Convergence appears to be stalled. Residual variation < 0.5%% over last 100 iterations.\n');
            fprintf('Consider adjusting under-relaxation factors or grid resolution.\n');
        end
    end
end

if ~converged
    fprintf('\nMaximum iterations reached without convergence.\n');
    fprintf('Final continuity residual: %.6e\n', avg_mass_imbalance);
    fprintf('Final u-momentum residual: %.6e\n', res_u);
    fprintf('Final v-momentum residual: %.6e\n', res_v);
end

time = toc;
fprintf('Total Time: ',time);

%% Post-processing
% Interpolate velocities to cell centers for visualization (matching MD Figure 2.1)
u_center = zeros(ny, nx);
v_center = zeros(ny, nx);

for i = 1:ny % Iterate over p-cell rows
    for j = 1:nx % Iterate over p-cell columns
        % u_center(i,j) is average of u at left face u(i+1,j) and right face u(i+1,j+1)
        u_center(i,j) = (u(i+1, j) + u(i+1, j+1)) / 2;
        
        % v_center(i,j) is average of v at bottom face v(i,j+1) and top face v(i+1,j+1)
        v_center(i,j) = (v(i, j+1) + v(i+1, j+1)) / 2;
    end
end

% Calculate stream function (at cell corners/nodes as per MD Figure 2.2)
psi = zeros(ny + 1, nx + 1);

% Set psi(1,1) as reference (bottom-left corner, x=0, y=0)
psi(1,1) = 0;

% Integrate along the bottom boundary (y=0) using v-velocities
for j = 1:nx
    psi(1, j+1) = psi(1, j) + v(1, j+1) * dx;
end

% Integrate upwards along the left boundary (x=0) using u-velocities
for i = 1:ny
    psi(i+1, 1) = psi(i, 1) - u(i+1, 1) * dy;
end

% Integrate for the rest of the internal points (from left to right, using v-fluxes)
for i = 1:ny
    for j = 1:nx
        psi(i+1,j+1) = psi(i+1,j) + v(i+1, j+1) * dx;
    end
end

% Calculate vorticity (at cell centers, similar to standard CFD)
omega = zeros(ny, nx);
for i = 1:ny
    for j = 1:nx
        % Using central differences on u_center and v_center
        % Handle boundary cells by clamping indices to stay within bounds for the difference.
        dv_dx = (v_center(i, min(j+1,nx)) - v_center(i, max(j-1,1))) / (2*dx);
        du_dy = (u_center(min(i+1,ny), j) - u_center(max(i-1,1), j)) / (2*dy);
        
        omega(i,j) = dv_dx - du_dy;
    end
end

% Extract centerline velocities
u_centerline = u_center(:, round(nx/2));
v_centerline = v_center(round(ny/2), :);

%% Plotting Results
% Plot streamlines (using the (ny+1, nx+1) psi grid and x_corners, y_corners)
x_corners = linspace(0, L, nx + 1);
y_corners = linspace(0, D, ny + 1);
[X_psi, Y_psi] = meshgrid(x_corners, y_corners);

figure(1);
contourf(X_psi, Y_psi, psi, 20); % Use the expanded psi grid for plotting
colormap('jet');
colorbar;
title(sprintf('Stream Function (Re = %d, L/D = %.2f)', Re, LD_ratio));
xlabel('x/L');
ylabel('y/D');
axis equal tight;

% Plot velocity vectors (using u_center, v_center on X, Y for cell centers)
figure(2);
quiver(X(1:2:end,1:2:end), Y(1:2:end,1:2:end), ...
        u_center(1:2:end,1:2:end), v_center(1:2:end,1:2:end), 2);
hold on;
contourf(X, Y, sqrt(u_center.^2 + v_center.^2), 20, 'LineStyle', 'none');
colormap('jet');
colorbar;
title(sprintf('Velocity Field (Re = %d, L/D = %.2f)', Re, LD_ratio));
xlabel('x/L');
ylabel('y/D');
axis equal tight;

% Plot horizontal velocity along vertical centerline
figure(3);
plot(u_centerline/U0, y_plot/D, 'b-', 'LineWidth', 2);
title(sprintf('Horizontal Velocity along Vertical Centerline (Re = %d, L/D = %.2f)', Re, LD_ratio));
xlabel('u/U_0');
ylabel('y/D');
grid on;

% Plot vertical velocity along horizontal centerline
figure(4);
plot(x_plot/L, v_centerline/U0, 'r-', 'LineWidth', 2);
title(sprintf('Vertical Velocity along Horizontal Centerline (Re = %d, L/D = %.2f)', Re, LD_ratio));
xlabel('x/L');
ylabel('v/U_0');
grid on;

% Plot residual history
figure(5);
semilogy(1:iteration, continuity_residual_hist(1:iteration), 'b-', 'LineWidth', 2, 'DisplayName', 'Continuity');
hold on;
semilogy(1:iteration, u_residual_hist(1:iteration), 'r--', 'LineWidth', 1.5, 'DisplayName', 'U-momentum');
semilogy(1:iteration, v_residual_hist(1:iteration), 'g-.', 'LineWidth', 1.5, 'DisplayName', 'V-momentum');
title(sprintf('Residual History (Re = %d, L/D = %.2f)', Re, LD_ratio));
xlabel('Iteration');
ylabel('Residual');
legend('Location', 'best');
grid on;


% Calculate center of cavity statistics for verification
center_i = round(ny/2);
center_j = round(nx/2);
center_u = u_center(center_i, center_j) / U0;
center_v = v_center(center_i, center_j) / U0;
% Pressure is directly at cell center indices
center_p = p(center_i, center_j);
fprintf('Normalized velocity at cavity center (u/U0, v/U0): (%.6f, %.6f)\n', center_u, center_v);
fprintf('Pressure at cavity center: %.6f\n', center_p);

% Save results for comparison
results = struct('u_center', u_center, 'v_center', v_center, ...
                'psi', psi, 'omega', omega, ...
                'u_centerline', u_centerline, 'v_centerline', v_centerline, ...
                'x', x_plot, 'y', y_plot, 'L', L, 'D', D, 'U0', U0, 'Re', Re, 'LD_ratio', LD_ratio);