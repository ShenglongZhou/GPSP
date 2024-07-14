clc, clear all, close all, warning off
addpath(genpath(pwd));

n    = 500;           % Signal dimension 
m    = ceil(0.5*n);   % Number of measurements
s    = ceil(0.01*n);  % Signal sparsity
r    = 0.05;          % Flipping ratio
k    = ceil(r*m);     % Upper bound of sign flips
v    = 0.5;           % Correlation parameter

Type = 'Ind';         % or 'Cor' 
test = 'v';           % change 'test' to see effect of GPSP to 
                      % factors {'s','m','r','v','n'}

switch test
  case 'm',   test0 = linspace(0.2,2,10);   
  case 's',   test0 = 2:10;  
  case 'r',   test0 = 0.02:0.02:0.2;  
  case 'v',   test0 = 0.1:0.1:0.9; Type = 'Cor';
  case 'n',   test0 = (5:5:20)*1e3;
end

f      = isequal(test,'n');  
S      = 200*(1-f)+10*f;
recd   = zeros(nnz(test0),4);
for j  = 1:nnz(test0)
    switch test
      case 'm',   m = ceil(test0(j)*n);
      case 's',   s = test0(j);
      case 'r',   r = test0(j);
      case 'v',   v = test0(j);    
      case 'n',   n = test0(j); s=ceil(0.01*n); m=ceil(n/2);
    end

    for ii = 1 : S    
        % Generate data 
        [A,c,co,xo] = random1bcs(Type,m,n,s,r,0.1,v);
        out         = GPSP(A,c,s,k);      
        recd(j,1)   = recd(j,1) - 20*log10(norm(xo-out.x));      
        recd(j,2)   = recd(j,2) + nnz(sign(A*out.x)-c)/m;
        recd(j,3)   = recd(j,3) + nnz(sign(A*out.x)-co)/m;
        recd(j,4)   = recd(j,4) + out.time;
    end
end

recd = recd/S; 

close all;
ylab = {'SNR','HD','HE','TIME'};
figure('Renderer', 'painters', 'Position', [200, 50, 900 190])
xloc = [ -0.09  -0.045  0.005  0.04];
for j  = 1:4
    sub  = subplot(1,4,j); 
    pos = get(sub, 'Position'); 
    tmp = recd(:,j);
    if j < 4
        plot(test0,tmp,'black.-','LineWidth',0.75), hold on,
        axis([min(test0) max(test0) 0-2*(j==1) 0.35+33*(j==1)]); 
    else
        semilogy(test0,tmp,'black.-','LineWidth',0.75), hold on,
    end
    grid on, xlabel(test), ylabel(ylab{j})   
    set(sub, 'Position',pos+[xloc(j),0,0.04,-0.01] )
end   
