clear all; close all; clc;

%% Parameters
L = 1;
T = 1;
Nx = 5;
Nt = 3;
dx = L/(Nx-1);
dt = T/(Nt-1);
lx = 0:dx:L;
lt = 0:dt:T;

% coarse discretization
lx_c = 0:dx*10:L;
lt_c = 0:dt*10:T;

% meshes 
% Meshes (fine and coarse)
[mesh_x_f,mesh_t_f] = meshgrid(lx,lt);
[mesh_x_g,mesh_t_g] = meshgrid(lx_c,lt_c);

% Boundary Conditions
ud_0 = sin(2*pi*lt/T);
ud_L = -sin(4*pi*lt/T);

 %Internal Force
f = zeros(Nx,Nt);
for i = 1:Nx-1
    for j = 1:Nt-1
        f(i,j)=10^3*sin(3*pi*lx(i)/T)*sin(5*pi*lt(j)/L);
    end
end
%fg = 10*rand(Nt/10,Nx/10);
%f = interp2(mesh_x_g,mesh_t_g,fg,mesh_x_f,mesh_t_f,'spline');

% "Stiffness" Matrix
k = (1/dx)*[1 -1;-1 1]; 
K = zeros(Nx);
for i = 1:Nx-1
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k;
end

%% Integration
% Space Integration 
% Elementary Matrix
M_elX = dx/6 * [2 1;1 2];

% Assembly on Space
M_asX = zeros(Nx);
for i = 1:Nx-1
    M_asX(i:i+1,i:i+1) = M_asX(i:i+1,i:i+1) + M_elX;
end

% Time Integration
M_elT = dt/6 * [2 1;1 2];
% Assembly on time
M_asT = zeros(Nt);
for i = 1:Nt-1
    M_asT(i:i+1,i:i+1) = M_asT(i:i+1,i:i+1) + M_elT;
end

%% Right Hand Side

F = M_asX*f;  

%% K matrix decomposition 
DOF = 1:Nx;
DOFd = [1,Nx];
DOFu = setdiff(DOF,DOFd);
Kuu = K(DOFu,DOFu);
Kud = K(DOFu,DOFd);
Fu = F(DOFu);
Ud = [ud_0',ud_L'];

%% Loop for every time step
U = zeros(Nx,Nt);
for i = 1:Nt
    Ft = F(:,i);
    U(DOFd,i) = Ud(i,:)';
    U(DOFu,i) = Kuu\(Ft(2:end-1)-Kud*Ud(i,:)');
end
mesh(U)  

%% Greedy Algorithm
Ucl = (1-lx/L)'*ud_0 + (lx/L)'*ud_L;
G = F-K*Ucl;
lambda_0 = lt;
it = 0;
W = zeros(Nx,Nt);
Modes = 2;
rec_error = zeros(1,Modes);

for mode = 1:Modes
    % Fixed Point Algorithm
    error = 1;
    while error > 1e-3
        H = (lambda_0*M_asT*lambda_0')*K;
        J = lambda_0*M_asT*G';
        Huu = H(DOFu,DOFu);
        Juu = J(DOFu);
        Lambda_k = zeros(Nx,1);
        Lambda_k(DOFu) = Huu\Juu';
        Lambda_k = Lambda_k./sqrt(Lambda_k'*K*Lambda_k);
        lambda_1 = Lambda_k'*G;
        error = ((lambda_1-lambda_0)*M_asT*(lambda_1-lambda_0)')/(lambda_0*M_asT*lambda_0');
        lambda_0 = lambda_1;
        it = it + 1;
    end
    G = G - K*Lambda_k*lambda_1;
    W = W + Lambda_k*lambda_1;
    U_greedy = Ucl + W;

    % Error Calculation 
    reconstuction_errorNum = zeros(1,Nx)';
    reconstuction_errorDen = zeros(1,Nx)';

    for k = 1:Nx
        reconstuction_errorNum(k) = (U(k,:)-U_greedy(k,:))*M_asT*(U(k,:)-U_greedy(k,:))';
        reconstuction_errorDen(k) = (U(k,:)*M_asT*U(k,:)');
    end
    rec_errorNum = reconstuction_errorNum'*M_asX*reconstuction_errorNum;
    rec_errorDen = reconstuction_errorDen'*M_asX*reconstuction_errorDen;
    rec_error(mode) = sqrt(rec_errorNum)/sqrt(rec_errorDen);
    
    % Visualization
    figure(1)
    subplot(1,3,1)
        mesh(U)
        title("Reference Solution")
    subplot(1,3,2)
        mesh(U_greedy);
        title("PGD ")
    subplot(1,3,3)
        plot(real(rec_error))
        xlabel('number of modes')
        ylabel('error')
    pause(0.5)
end
