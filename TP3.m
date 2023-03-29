clear all; close all; clc;

%% Parameters
L = 1;
T = 1;
Nx = 1000;
Nt = 100;
dx = L/(Nx-1);
dt = T/(Nt-1);
lx = 0:dx:L;
lt = 0:dt:T;
alpha = 5;

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

% Internal Force
f = zeros(Nx,Nt);
% for i = 1:Nx-1
%     for j = 1:Nt-1
%         f(i,j)=10^3*sin(3*pi*lx(i)/T)*sin(5*pi*lt(j)/L);
%     end
% end
fg = 100*rand(Nt/10,Nx/10);
%f = interp2(mesh_x_g,mesh_t_g,fg,mesh_x_f,mesh_t_f,'spline');

%% Matrices

% "Stiffness" Matrix
k = (1/dx)*[1 -1;-1 1]; 
K = zeros(Nx);
for i = 1:Nx-1
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k;
end

% Integration
% Space Integration 
% Elementary Matrix
M_elX = dx/6 * [2 1;1 2];

% Assembly on Space
Ix = zeros(Nx);
for i = 1:Nx-1
    Ix(i:i+1,i:i+1) = Ix(i:i+1,i:i+1) + M_elX;
end

% Time Integration
M_elT = dt/6 * [2 1;1 2];
% Assembly on time
It = zeros(Nt);
for i = 1:Nt-1
    It(i:i+1,i:i+1) = It(i:i+1,i:i+1) + M_elT;
end

% da/dt * b integration 
M_elt = 1/2 * [-1 -1;1 1];
Ft = zeros(Nt);
for i = 1:Nt-1
    Ft(i:i+1,i:i+1) = Ft(i:i+1,i:i+1) + M_elt;
end

% Time Derivative Matrix (dx/dt = D * x) 
D = It\Ft';
%test = sin(2*pi*lt');
%testdev = 2*pi*cos(2*pi*lt');
%testhat = D*test;
%figure;hold on; plot(testdev); plot(testhat,"*"); hold off;
%% Brute Force Quantities --> P*u_k = F + M/dt * u_(k-1)
F = Ix*f;  
% P matrix
P = Ix/dt + K;
%% Matrix decomposition (To account for Boundary Conditions)
DOF = 1:Nx;
DOFd = [1,Nx];
DOFu = setdiff(DOF,DOFd);
Puu = P(DOFu,DOFu);
Pud = P(DOFu,DOFd);
Ud = [ud_0',ud_L'];
Muu = Ix(DOFu,DOFu);
%% Brute Force Solution -->  Puu*U_u + Pud*U_d = Q
U = zeros(Nx,Nt);
for i = 2:Nt
    Ft = F(:,i);
    U(DOFd,i) = Ud(i,:)';
    Q = Ft(2:end-1) + Muu/dt*U(DOFu,i-1);
    U(DOFu,i) = Puu\(Q-Pud*Ud(i,:)');
end
figure("Name","Brute Force Solution");mesh(U);
xlabel("Time(s)");ylabel("Position(m)");zlabel("U(x,t)");title("Brute Force Solution")

%% Greedy Algorithm
% Boundary Conditions
Ucl = (1-lx/L)'*ud_0 + (lx/L)'*ud_L;

% Right Hand side --> M*dw/dt + K*w = G
G = F-K*Ucl-Ix*(D*Ucl')';

% Initialization for Fixed Point 
Modes = 1;
lambda_0 = lt;
W = zeros(Nx,Nt);
rec_error = zeros(1,Modes); 
rec_errorSVD = zeros(1,Modes);
lambda = [];
Lambda = [];

% SVD Reconstruction 
[X,S,V] = svd(U);
U_SVD = zeros(Nt,Nx);
it2 = 0;
for mode = 1:Modes
    % Fixed Point Algorithm
    error = 1;
    it = 0;
    while error > 1e-3 
        % System in the form --> H*Lambda_k = J
        H = (lambda_0*It*lambda_0')*K + (lambda_0*It*(D*lambda_0'))*Ix ;
        J = lambda_0*It*G';

        % We account for the Boundary Conditions --> zeros at the borders
        Huu = H(DOFu,DOFu);
        Juu = J(DOFu);
        Lambda_k = zeros(Nx,1);
        Lambda_k(DOFu) = Huu\Juu';

        % Normalization for unicity
        Lambda_k = Lambda_k./sqrt(Lambda_k'*K*Lambda_k);

        % Solve for lambda --> lambda_1'*(m*D + I) = h
        % Account for initial conditions --> 0 at t_0 --> lambda(1) = 0
        lambda_1 = zeros(Nt,1);
        h = (Lambda_k'*G)';
        m = (Lambda_k'*Ix*Lambda_k);
        lambda_1(2:end) = (m*D(2:end,2:end) + eye(size(D(2:end,2:end))))\h(2:end);
        lambda_1 = lambda_1';

        % Stagnation Criteria 
        error = ((lambda_1-lambda_0)*It*(lambda_1-lambda_0)')/(lambda_0*It*lambda_0');

        % Update lambda for next iteration
        lambda_0 = lambda_1;
        it = it + 1;
        if it > 30
            break 
        end
    end
    
    %Lambda = [Lambda Lambda_k];
     % Orthogonalisation
     % Gram-Schmidt algorithm.

     %for p = 1:mode-1
     %   Lambda_k_orth = Lambda_k - sqrt(Lambda_k'*Ix*Lambda(:,p))*Lambda(:,p);
        %lambda(:,p) = lambda(:,p) + lambda_new*sqrt(Gamma_new'*M_asX*Gamma(:,p));
    % end
    %Lambda_k_orth = 
    %Lambda(:,mode)  = Lambda_k_orth/sqrt(Lambda_k_orth'*K*Lambda_k_orth);
    %lambda(:,mode) = Lambda_new*sqrt(Gamma_new_orth'*M_asX*Gamma_new_orth);

    % Update Right Hand Side 
    G = G - K*Lambda_k*lambda_1-Ix*Lambda_k*(D*lambda_1')';

    % Update PGD Solution
    W = W + Lambda_k*lambda_1;
    U_greedy = Ucl + W;
    
    % SVD Reconstruction
    U_SVD = X(:,1:mode)*S(1:mode,1:mode)*V(:,1:mode)';
    
    % Error Calculation 
    reconstuction_errorNum = zeros(1,Nx)';
    reconstuction_errorDen = zeros(1,Nx)';
    reconstuction_errorNumSVD = zeros(1,Nx)';

    for k = 1:Nx
        reconstuction_errorNum(k) = (U(k,:)-U_greedy(k,:))*It*(U(k,:)-U_greedy(k,:))';
        reconstuction_errorDen(k) = (U(k,:)*It*U(k,:)');
        reconstuction_errorNumSVD(k) = (U(k,:)-U_SVD(k,:))*It*(U(k,:)-U_SVD(k,:))';
    end
    rec_errorNum = reconstuction_errorNum'*Ix*reconstuction_errorNum;
    rec_errorDen = reconstuction_errorDen'*Ix*reconstuction_errorDen;
    rec_error(mode) = sqrt(rec_errorNum)/sqrt(rec_errorDen);

    rec_errorNumSVD = reconstuction_errorNumSVD'*Ix*reconstuction_errorNumSVD;
    rec_errorSVD(mode) = sqrt(rec_errorNumSVD)/sqrt(rec_errorDen);

    % Visualization
    figure(1)
    subplot(1,3,1)
        mesh(U)
        title("Reference Solution")
        xlabel("Time(s)");ylabel("Position(m)");zlabel("U(x,t)")
    subplot(1,3,2)
        mesh(U_greedy);
        title("PGD Greedy Solution")
        xlabel("Time(s)");ylabel("Position(m)");zlabel("U(x,t)")
    subplot(1,3,3)
        semilogy(real(rec_error),"Color","red") 
        hold on 
        semilogy(real(rec_errorSVD),"Color","Blue")
        hold off
        legend("PGD","SVD")
        xlabel('number of modes')
        ylabel('error')
        title("Error Greedy")
    pause(0.5)
    it2 = it2 + 1;
end