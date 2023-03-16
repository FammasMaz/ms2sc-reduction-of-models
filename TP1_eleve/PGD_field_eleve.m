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
subplot(1,4,1)
surf(xc,tc,Ug)
title("Coarse distribution")
shading interp

subplot(1,4,2)
surf(xf,tf,Uf)
title("Fine distribution")
shading interp
%% Elementary matrix and assembly for integration 

%TO DO
mx = dx/6*[2 1; 1 2];
Ix = zeros(n_X, n_X);
It = zeros(n_T, n_T);

for i=1:n_X-1
    Ix(i:i+1, i:i+1) = Ix(i:i+1, i:i+1) + mx;
end

mt = dt/6*[2 1; 1 2];
for i=1:n_T-1
    It(i:i+1, i:i+1) = It(i:i+1, i:i+1) + mt;
end
%% PGD Approximation

% Number of modes
n_modes    = 20;

% initialisation
Lambda     = zeros(n_T,n_modes);
Gamma      = zeros(n_X,n_modes);
U_PGD      = zeros(n_T,n_X);
PGD_error = zeros(n_modes,1);
Ustar = zeros(n_T,n_X);
error = []

%% 
% Greedy algorithm + fixed point 
for mode = 1:n_modes
    
    U_PGD = Lambda*Gamma';


    % Add new mode
    
    Ustar = Uf - U_PGD;
    Lambdai = linspace(0, T, n_T)';
    s = 1;
    % Fixed point algorithm
    while s>(10e-3)        
        
        % Save old function Lambda
        Lambdaold = Lambdai;

        
        % Compute new space function Gamma
        Gammai = Ustar'*It*Lambdai/(Lambdai'*It*Lambdai);
        
        % Normalisation
        Gammai = Gammai/sqrt((Gammai'*Ix*Gammai));
        
        % Compute new space function Lambda
        Lambdai = Ustar*Ix*Gammai/(Gammai'*Ix*Gammai);
        
        % update PGD error for stagnation criteria
        s = (Lambdai - Lambdaold)'*It*(Lambdai - Lambdaold)/(Lambdaold'*It*Lambdaold);
        mode
    end
    


    % Orthogonalisation
     %We'll do this after we've successfully run the greedy
     % algorithm and managed to calculate the error in 
     % space-time. All we need to do here is to implement 
     % the Gram-Schmidt algorithm.
     % Note that mode orthogonalization is not useful here,
     % but it will be in the next TP 
     % (a priori use of the PGD).
    Gammanew = Gammai;
    Gammanew_orth = Gammai;
    Lambdanew = Lambdai;
    for m=1:mode-1
        Gammanew_orth = Gammanew_orth - sqrt(Gammanew'*Ix*Gamma(:, m))*Gamma(:, m);
        Lambda(:, m) = Lambda(:, m) + Lambdanew*sqrt(Gammanew'*Ix*Gamma(:,m));
    end

    Gamma(:, mode) = Gammanew/sqrt(Gammanew_orth'*Ix*Gammanew_orth);
    Lambda(:, mode) = Lambdanew * sqrt(Gammanew_orth'*Ix*Gammanew_orth);

    % compute the new PGD approximation with the added mode
    U_PGD = Lambda*Gamma';
    
    % Calculation of the reconstruction error with the sptio-temporal norm
    num = zeros(n_X, 1);
    den = zeros(n_X, 1);
    for i=1:n_X
        num(i) = (Uf(:, i) - U_PGD(:,i))'*It*(Uf(:, i) - U_PGD(:,i));
        den(i) = (Uf(:, i))'*It*(Uf(:, i));
    end
    num = num'*Ix*num;
    den = den'*Ix*den;
    error = [error num/den];
    subplot(1,4,3)
    surf(mesh_x_f,mesh_t_f,U_PGD)
    xlabel('x')
    ylabel('y')
    zlabel('U PGD')
    shading interp
    Title = ["Approximation with "+ mode " mode(s)"];
    title(Title)
    subplot(1,4,4);
    plot(error);
    xlabel('modes');
    ylabel('error');
    Title2 = ["Error with "+ mode " mode(s)"];
    title(Title2)    
    pause(0.5)
end

saveas(figure(1),'Results.png');

