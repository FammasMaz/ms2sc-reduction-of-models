clear; close all; clc;

%% Global parameters
E = 1;
S = 1;
if_force = 0;
nmodes = 15;

%% Grid definition

L = 1.0; % length of the domain
T = 1.0; % final time
nx = 1000; % number of grid points
nt = 100; % number of time steps

lx0 = 0.0; % left boundary
lt0 = 0.0; % initial time

dx = L/(nx-1); % grid spacing
dt = T/(nt-1); % time step

Ix = sparse(nx, nx);
It = sparse(nt, nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time steps
lx = lx0:dx:L; % grid points
lt = lt0:dt:T; % time points

%% Matrix initialization
f = sparse(nx, nt); % force matrix 
K = sparse(nx, nx); % stiffness matrix
u = sparse(nx, nt); % solution matrix

%% discrete initial condition
ud0d = sin(2*pi*lt)/T; % initial condition
udLd = -sin(4*pi*lt)/T; % boundary condition

%% discrete force
if if_force==1
    f = (1e3 * (sin(3*pi*lt)/T)' * (sin(5*pi*lx)/L))';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constructing the Identities for Assembly

mx = dx/6*[2 1; 1 2];
mt = dt/6*[2 1; 1 2];

for i = 1:nx-1
    Ix(i:i+1, i:i+1) = Ix(i:i+1, i:i+1) + mx;
end
for i = 1:nt-1
    It(i:i+1, i:i+1) = It(i:i+1, i:i+1) + mt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assembly of F and K

F = Ix * f;
ke = 1/dx*[1 -1; -1 1];

for i = 1:nx-1
    K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + ke;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refactoring the solution matrix
u(1, :) = ud0d;
u(end, :) = udLd;
dof_b = [1, nx];
dof_u = setdiff(1:nx, dof_b);
Kuu = K(dof_u, dof_u);
Kub = K(dof_u, dof_b);
Kbu = K(dof_b, dof_u);
Kbb = K(dof_b, dof_b);

Fu = F(dof_u, :);
u(dof_u, :) = Kuu\(Fu - Kub*u(dof_b, :));

figure
[lx_mesh, lt_mesh] = meshgrid(lx, lt);
mesh(lx_mesh, lt_mesh, u');
xlabel('x')
ylabel('t')
zlabel('u')
title(['Solution of the 1D wave equation with force = ', num2str(if_force)]);
saveas(gcf, strcat('../Final Report/assets/TP2_ref_solution_', num2str(if_force), '.png'));
saveas(gcf, strcat('assets/TP2_ref_solution_', num2str(if_force), '.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot collectively
figure()
subplot(1,3,1)
surf(lx_mesh,lt_mesh,u')
title("Coarse distribution")
shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reduction of the model PGD

ucl = (1-lx/L)'*ud0d + (lx/L)'*udLd;
G = F - K*ucl;

lambda = lt;
er = 1;
iter = 0;
nb_modes = 15;
Lambda = zeros(nt,nb_modes);
Gamma = zeros(nx,nb_modes);
for i = 1:nb_modes
    G - K*Lambda*lambda;
    while er > 1e-3
        old_lambda = lambda;
        iter = iter + 1;
        intlambda = lambda * It * lambda';
        H = intlambda * K;
        J = lambda * It * G';
        Huu = H(dof_u, dof_u);
        Ju = J(:, dof_u);
        Lambda = zeros(nx, 1);
        Lambda(dof_u) = (Huu\Ju')';
        Lambda = Lambda / sqrt(Lambda'*K*Lambda);
        lambda = (Lambda' * G);
        er = ((lambda - old_lambda)*It*(lambda - old_lambda)')/ intlambda;
    end
    
fprintf('\nNumber of iterations: %d', iter);
fprintf('\n\nFinal error: %d', er);


%% Cost = no.modes * n_iter * cost
u_red = ucl + Lambda * lambda;
figure
[lx_mesh, lt_mesh] = meshgrid(lx, lt);
mesh(lx_mesh, lt_mesh, u_red');
xlabel('x')
ylabel('t')
zlabel('u')
title(['Reduced Solution of the 1D wave equation with force = ', num2str(if_force)]);
saveas(gcf, strcat('../Final Report/assets/TP2_red_solution_', num2str(if_force), '.png'));
saveas(gcf, strcat('assets/TP2_red_solution_', num2str(if_force), '.png'));

