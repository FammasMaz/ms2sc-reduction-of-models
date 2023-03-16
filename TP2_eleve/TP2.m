%% Grid definition

L = 1.0 % length of the domain
T = 1.0 % final time
nx = 1000; % number of grid points
nt = 1000; % number of time steps

lx0 = 0.0; % left boundary
lt0 = 0.0; % initial time

dellx = L/(nx-1); % grid spacing
dellt = T/(nt-1); % time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time steps
lx = lx0:dellx:L; % grid points
lt = lt0:dellt:T; % time points


f = zeros(nx, nt); % solution matrix 

%% Continuous initial condition
ud0 = @(t) sin(2*pi*t)/T; % initial condition
udL = @(t) -sin(4*pi*t)/T; % boundary condition

%% discrete initial condition
ud0d = sin(2*pi*lt)/T; % initial condition
udLd = -sin(4*pi*lt)/T; % boundary condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


