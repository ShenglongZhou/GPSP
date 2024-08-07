clc; close all; clear; warning off
addpath(genpath(pwd));

n          = 2000;          % Signal dimension 
m          = ceil(0.5*n);   % Number of measurements
s          = ceil(0.01*n);  % Sparsity level
r          = 0.01;          % Flipping ratio
type       = 'Ind';         % or 'Cor' 
[A,c,co,xo]= random1bcs(type,m,n,s,r);

k          = ceil(0.01*m);
out        = GPSP(A,c,s,k);
x          = out.x; 
fprintf('Time: %6.3f\n',out.time);
fprintf('SNR:  %6.3f\n',-10*log10(norm(x-xo)^2));
fprintf('HD:   %6.3f\n',nnz(sign(A*x)-c)/m)
fprintf('HE:   %6.3f\n',nnz(sign(A*x)-co)/m)
RecoverShow(xo,x,[1000 450 500 250])
