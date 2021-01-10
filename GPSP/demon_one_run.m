clc; close all; clear; warning off

n          = 1000; 
m          = ceil(0.5*n);
s          = ceil(0.01*n);  % sparsity level
r          = 0.01;          % flipping ratio
type       = 'Ind';         % 'Ind' or 'Cor' 
[A,c,co,xo]= random1bcs(type,m,n,s,r);

k          = ceil(0.01*m);
out        = GPSP(A,c,s,k);
x          = out.x; 

fprintf('Time: %6.3f\n',out.time);
fprintf('SNR:  %6.3f\n',-10*log10(norm(x-xo)^2));
fprintf('HD:   %6.3f\n',nnz(sign(A*x)-c)/m)
fprintf('HE:   %6.3f\n',nnz(sign(A*x)-co)/m)

ReoveryShow(xo,x,[1000, 550, 400 200],1),pause(0.05)
figure('Renderer', 'painters', 'Position', [1000, 200, 400 200])
semilogy(1:out.iter,  out.OBJ,'b:','LineWidth',2), hold on
xlabel('iteration'), ylabel('objective')
grid on, legend('GPSP')


