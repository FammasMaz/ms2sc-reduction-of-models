%% PGD Approximation of a random space-time field 
%% - A POSTERIORI-
% _V. Matray_

clear all
close all
clc
%% 
% Time discretisation

T   = 10; 
n_T = 400; 
tf  = linspace(0,T,n_T);     % fine discretisation
tc  = linspace(0,T,n_T/10);  % coarse discretisation
dt  = tf(2); % time step
%% 
% Space discretisation

L = 1; 
n_X = 200; 
xf  = linspace(0,T,n_X);    % fine discretisation
xc  = linspace(0,T,n_X/10); % coarse discretisation
dx  = xf(2);
%% 
% Display of the field smoothed to approximated by PGD 

% Meshes (fine and coarse)
[mesh_x_f,mesh_t_f] = meshgrid(xf,tf);
[mesh_x_g,mesh_t_g] = meshgrid(xc,tc);

% Random function on the coarse mesh:
Ug = rand(n_T/10,n_X/10);

% Interpolation on the fine mesh to get a 
% smooth distribution : (If the field is too discontinuous
% the PGD will have trouble approximating it: 
% you can check this at the end of the session).

Uf = interp2(mesh_x_g,mesh_t_g,Ug,mesh_x_f,mesh_t_f,'spline');

% Display distribution
figure()
subplot(1,3,1)
surf(xc,tc,Ug)
title("Coarse distribution")
shading interp

subplot(1,3,2)
surf(xf,tf,Uf)
title("Fine distribution")
shading interp
%% Elementary matrix and assembly for integration 

%TO DO
% Space Integration 
% Elementary Matrix
M_elX = dx/6 * [2 1;1 2];

% Assembly on Space
M_asX = zeros(n_X);
for i = 1:n_X-1
    M_asX(i:i+1,i:i+1) = M_asX(i:i+1,i:i+1) + M_elX;
end

% Time Integration
M_elT = dt/6 * [2 1;1 2];
% Assembly on time
M_asT = zeros(n_T);
for i = 1:n_T-1
    M_asT(i:i+1,i:i+1) = M_asT(i:i+1,i:i+1) + M_elT;
end

%% PGD Approximation

% Number of modes
n_modes    = 200;

% initialisation
Lambda     = zeros(n_T,n_modes);
Gamma      = zeros(n_X,n_modes);
U_PGD      = zeros(n_T,n_X);
PGD_error = zeros(n_modes,1);
%% 
% Greedy algorithm + fixed point 
modes = 14;
rec_error = zeros(1,modes);
for mode = 1:modes
    
    U_PGD = Lambda*Gamma';
    U_star = Uf - U_PGD;
    % Add new mode
    % Compute Gamma new
    Lambda_new = linspace(0,n_T*dt,n_T)';
    s_error = 1;
    
    % Fixed point algorithm
    fix_point_iter = 0;
    while s_error > 0.01
    %for test = 1:1
        
        % Save old function Lambda
        Lambda_old = Lambda_new;
        
        % Compute new space function Gamma
        Gamma_new = (U_star'*M_asT*Lambda_new) / (Lambda_new'*M_asT*Lambda_new) ;

        
        % Normalisation
        Gamma_new = Gamma_new/sqrt(Gamma_new'*M_asX*Gamma_new);
        
        % Compute new time function Lambda
        Lambda_new = U_star*M_asX*Gamma_new / (Gamma_new'*M_asX*Gamma_new);
        
        % update PGD error for stagnation criteria
        s_error = (Lambda_new-Lambda_old)'*M_asT*(Lambda_new-Lambda_old)/(Lambda_old'*M_asT*Lambda_old);
        fix_point_iter = fix_point_iter + 1;
    end
    

    % Orthogonalisation
     %We'll do this after we've successfully run the greedy
     % algorithm and managed to calculate the error in 
     % space-time. All we need to do here is to implement 
     % the Gram-Schmidt algorithm.
     % Note that mode orthogonalization is not useful here,
     % but it will be in the next TP 
     % (a priori use of the PGD).

     Gamma_new_orth =  Gamma_new;
     for p = 1:mode-1
        Gamma_new_orth = Gamma_new_orth - sqrt(Gamma_new'*M_asX*Gamma(:,p))*Gamma(:,p);
        Lambda(:,p) = Lambda(:,p) + Lambda_new*sqrt(Gamma_new'*M_asX*Gamma(:,p));
     end
    
    
    Gamma(:,mode)  = Gamma_new/sqrt(Gamma_new_orth'*M_asX*Gamma_new_orth);
    Lambda(:,mode) = Lambda_new*sqrt(Gamma_new_orth'*M_asX*Gamma_new_orth);

    % compute the new PGD approximation with the added mode
    U_PGD = Lambda*Gamma';
    
    % Calculation of the reconstruction error with the sptio-temporal norm
    reconstuction_errorNum = zeros(1,n_X)';
    reconstuction_errorDen = zeros(1,n_X)';

    for k = 1:n_X
        reconstuction_errorNum(k) = (Uf(:,k)-U_PGD(:,k))'*M_asT*(Uf(:,k)-U_PGD(:,k));
        reconstuction_errorDen(k) = (Uf(:,k)'*M_asT*Uf(:,k));
    end
   
    rec_errorNum = reconstuction_errorNum'*M_asX*reconstuction_errorNum;
    rec_errorDen = reconstuction_errorDen'*M_asX*reconstuction_errorDen;
    rec_error(mode) = sqrt(rec_errorNum)/sqrt(rec_errorDen);

    subplot(1,3,1)
         surf(mesh_x_f,mesh_t_f,real(U_PGD))
         xlabel('x')
         ylabel('y')
         zlabel('U PGD')
         shading interp
         Title = ["Approximation with "+ mode " mode(s)"];
         title(Title)
      subplot(1,3,2)
         surf(mesh_x_f,mesh_t_f,real(Uf))
         xlabel('x')
         ylabel('y')
         zlabel('U exact')
         shading interp
         Title = ["Real U"];
         title(Title)
     subplot(1,3,3)
        plot(real(rec_error),"+")
        xlabel('number of modes')
        ylabel('error')
     pause
end