clear all; close all; clc;

%% Parameters
L = 1;
T = 1;
Nx = 1000;
Nt = 500;
dx = L/(Nx-1);
dt = T/(Nt-1);
lx = 0:dx:L;
lt = 0:dt:T;
alpha = 1;

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
% f = zeros(Nx,Nt)';
%for i = 1:Nx-1
%    for j = 1:Nt-1
%        f(i,j)=10^3*sin(3*pi*lx(i)/T)*sin(5*pi*lt(j)/L);
%    end
%end
fg = 10*rand(Nt/10,Nx/10);
f = interp2(mesh_x_g,mesh_t_g,fg,mesh_x_f,mesh_t_f,'spline');

% Matrices

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
% Brute Force Quantities --> P*u_k = F + M/dt * u_(k-1)
F = Ix*f';  
% P matrix
P = Ix/dt + K;
% Matrix decomposition (To account for Boundary Conditions)
DOF = 1:Nx;
DOFd = [1,Nx];
DOFu = setdiff(DOF,DOFd);
Puu = P(DOFu,DOFu);
Pud = P(DOFu,DOFd);
Ud = [ud_0',ud_L'];
Muu = Ix(DOFu,DOFu);
% Brute Force Solution -->  Puu*U_u + Pud*U_d = Q
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
G_fix =  F - K*Ucl - Ix*(D*Ucl')';
U_greedy = Ucl;

% Initialization for Fixed Point 
Modes = 10;
lambda_0 = lt;
W = zeros(Nx,Nt);
rec_error = zeros(1,Modes); 
rec_errorSVD = zeros(1,Modes);
lambda = [];
Lambda = [];
Lambda_sine_t = [];
new_Lam = [];
G_var = [];

% SVD Reconstruction 
[X,S,V] = svd(U);
U_SVD = zeros(Nt,Nx);
it2 = 0;
update = 1;    

for mode = 1:Modes
    G_var(mode,:,:) = G;
    % Fixed Point Algorithm
    error = 1;
    it = 0;
    while error > 1e-8
        % System in the form --> H*Lambda_k = J
        H = (lambda_0*It*lambda_0')*K + (lambda_0*It*(D*lambda_0'))*Ix ;
        J = lambda_0*It*G';

        % We account for the Boundary Conditions --> zeros at the borders
        Huu = H(DOFu,DOFu);
        Juu = J(DOFu);
        Lambda_k = zeros(Nx,1);
        Lambda_k(DOFu) = Huu\Juu';

        % Normalization for unicity
%         Lambda_k = Lambda_k./sqrt(Lambda_k'*K*Lambda_k);

        % Solve for lambda --> lambda_1'*(m*D + I) = h
        % Account for initial conditions --> 0 at t_0 --> lambda(1) = 0
        lambda_1 = zeros(Nt,1);
        h = (Lambda_k'*G)';
        m = (Lambda_k'*Ix*Lambda_k);
        i_m = (Lambda_k'*K*Lambda_k);
        lambda_1(2:end) = (m*D(2:end,2:end) + i_m*eye(size(D(2:end,2:end))))\h(2:end);

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

if update == 1
%     G = G - K*Lambda_k*lambda_1-Ix*Lambda_k*(D*lambda_1')';
    Lambda = gram_schmidt([Lambda, Lambda_k], K);

%         % Update Phase
    m_update = Lambda'*Ix*Lambda;
    k_update = Lambda'*K*Lambda;
    F_update = Lambda'*G_fix;
%     disp(size(F_update))
    
% Update Right Hand Side 
    [~, lambda_upd] = ode45(@(t,y) deriv(t, y, m_update, k_update, F_update, dt), lt, zeros(mode,1));


% lambda_upd = zeros(mode, Nt);
% 
% for i = 1:Nt-1
% 
% %     lambda_upd(:,i+1) = m_update\(F_update(:,i) - k_update*lambda_upd(:,i))*dt + lambda_upd(:,i);
%     lambda_upd(:,i+1) = (1/dt*m_update + k_update)\(F_update(:,i) + 1/dt * m_update * lambda_upd(:,i));
% 
% end
% lambda_upd = lambda_upd';

    U_greedy = Ucl + Lambda*lambda_upd';
%     G = G - K*Lambda(:,end)*lambda_upd(:,end)'- Ix*Lambda(:,end)*(D*lambda_upd(:,end))';
    G = F-K*Ucl-Ix*(D*Ucl')'  -  K*Lambda*lambda_upd'- Ix*Lambda*(D*lambda_upd)';
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

% end
else
    Lambda = [Lambda, Lambda_k];
    lambda = [lambda; lambda_1];

    % Update Right Hand Side 
    G =  F - K*Ucl - Ix*(D*Ucl')' - K*(Lambda*lambda) - Ix*(Lambda*(D*lambda')');

    % Update PGD Solution
%     W = W + Lambda_k*lambda_1;
    U_greedy = Ucl + Lambda*lambda;
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
%     rec_error(mode) = norm(U-U_greedy)/norm(U);

    rec_errorNumSVD = reconstuction_errorNumSVD'*Ix*reconstuction_errorNumSVD;
    rec_errorSVD(mode) = sqrt(rec_errorNumSVD)/sqrt(rec_errorDen);
%     rec_errorSVD(mode) = norm(U-U_SVD)/norm(U);

end

    it2 = it2 + 1;
    disp(['Mode n ', num2str(mode) ])
end

%  % Visualization
    figure(1)
    subplot(1,3,1)
        mesh(U)
        title("Fine Solution")
    subplot(1,3,2)
        mesh(U_greedy);
        title("PGD Greedy Solution")
    subplot(1,3,3)
        semilogy(real(rec_error),"Color","red") 
        hold on 
        semilogy(real(rec_errorSVD),"Color","Blue")
        hold off
        if update == 1
            legend("PGD with update","SVD")
        else
            legend("PGD without update","SVD")
        end

        xlabel('number of modes')
        ylabel('error')
        title("Error Greedy")

% Ortho verif

ortho_verif_K = zeros(length(Lambda(1,:)));

for i = 1:length(ortho_verif_K)
    for j = i:length(ortho_verif_K)
        ortho_verif_K(i, j) = Lambda(:, i)'*K*Lambda(:,j);
    end
end
 
ortho_verif_M = zeros(length(Lambda(1,:)));

for i = 1:length(ortho_verif_M)
    for j = i:length(ortho_verif_M)
        ortho_verif_M(i, j) = Lambda(:, i)'*Ix*Lambda(:,j);
    end
end

%semilogy(rec_errorSVD);
hold on
semilogy(rec_error);
eta_test = rec_error;

semilogy(eta_test)
legend("SVD", "PGD: After Update",'PGD: before update')
xlabel('Number of modes')
ylabel('Error')
title(["Error Comparison"])


%% TEST


x = 1:1000;
a = 2;
y = sin(a*pi*x/1000);
figure(1)
plot(x,y)

%% Functions
function Lam = update_phase(Lambda, Lambda_tot, i)
% correction = zeros(size(Lambda));
        for j = 1:i - 1
            Lambda = Lambda - Lambda*(Lambda'*Lambda_tot(:,j))./2;
        end
        Lam = Lambda;
end

function Q = gram_schmidt(A, K)
% A: an mxn matrix
% Q: an mxn matrix with orthonormal columns
[m, n] = size(A);
Q = zeros(m, n);

for j = 1:n
    v = A(:, j);
    for i = 1:j-1
        q = Q(:, i);
        v = v - q*(q'*v);
    end
    Q(:, j) = v/norm(v);
end
Q = Q/sqrtm(Q'*K*Q);
end

function [Q,R] = orthonormalize(A,M,K)
% ORTHONORMALIZE Orthonormalize a matrix A according to two matrices M and K
% using the Gram-Schmidt process.
%
%   [Q,R] = ORTHONORMALIZE(A,M,K) returns an orthonormal matrix Q and an
%   upper triangular matrix R such that A = Q*R and Q'*M*Q = I and Q'*K*Q = I.
%
%   Example:
%      A = randn(5);
%      M = randn(5);
%      K = randn(5);
%      [Q,R] = orthonormalize(A,M,K);
%
%   See also GRAM-SCHMIDT PROCESS.

[~,n] = size(A);
Q = zeros(size(A));
R = zeros(n);

for j = 1:n
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

[~,R1] = qr(Q'*M*Q);
[~,R2] = qr(Q'*K*Q);

if size(R1,2) < size(Q,2)
    Q(:,size(R1,2)+1:end)=[];
end

if size(R2,2) < size(Q,2)
    Q(:,size(R2,2)+1:end)=[];
end

end

function deriv = deriv(t, lam, m, k, ht, dt)
%     h = htim(ht, t, lt);
    h = ht(:,round(t/(dt) + 1));
%     h = ht(:,round(t/(t(2)) + 1));
deriv = (m^(-1))*(h - k*lam);
%     deriv = m\(h - eye(size(k))*lam);
%     deriv = deriv';
% deriv = deriv';
end



