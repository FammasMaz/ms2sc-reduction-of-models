%% *POD applied to an image*
%% 
% _V. Matray_

clear all
close all
clc
%% POD by SVD
% load original image

UCimage = imread('callas.jpg');
UCimage = im2double(UCimage(:,:,1));
figure
subplot(1,2,1)
colormap(gray); imagesc(UCimage(:,:));axis equal;axis tight; axis off;
title ('Maria Callas')
%% 
% get dimensions of the image

n_T = size(UCimage, 1);
n_M = size(UCimage, 2);
U = UCimage;
%% 
% SVD on the image

[phi, S, V] = svd(U);
%% 
% Order for reconstruction ?

order = 25;
%% 
% Initialize the reconstruction 

Ui = zeros(n_M,n_T)';
%% 
% Reconstruction :

for o = 1:order
    for i = 1:n_T
       Ui(:,i) =  Ui(:,i)+ phi(:,o)*(U(:,i)'*phi(:,o)) ;
    end
Ui = Ui +S(o,o)*phi(:,o)*V(:,o)';
end
%% 
% Optional challenge: build Ui using a single Matlab line 

%TODO
%% 
% plot

subplot(1,2,2)
colormap(gray); imagesc(Ui(:,:));axis equal;axis tight; axis off;
title("Approximation avec "+order+" modes")
%% Calculation of the reconstruction error
% Plot the reconstruction error as a function of the number of modes used. The 
% reconstruction error considered here will be calculated with the norm of your 
% choice (you can use the "norm" function of Matlab it's the fastest/ easiest 
% way ;) ). 

%TODO
M = (S*phi*V');
Mord = zeros(size(M));
for o = 1:order
    Mord = Mord + S(o,o)*phi(:,o)*V(:,o)'
end

e = norm(Mord - M)/norm(M);
%% 
% The evolution of the reconstruction error

order = 25
e = []
sior = n_M*n_T;
sizes = [];
for o=1:50    
    Ui = zeros(size(M));
    Ui(:,:) = phi(:, 1:o)*S(1:o, 1:o)*V(:, 1:o)';
    e = [e , norm(U-Ui)/norm(U)];
    sizes = [sizes sior-(o*(n_T+n_M+1))];
end
hold off
figure
plot(e)