clc, clear all, close all, warning off
addpath(genpath(fileparts(mfilename('fullpath'))));
 
n    = 500;           % Signal dimension 
m    = ceil(0.5*n);   % Number of measurements
s    = ceil(0.01*n);  % Signal sparsity
r    = 0.05;          % Probability of sign flips
k    = ceil(r*m);     % Upper bound of sign flips
v    = 0.5;

Type = 'Ind';         % 'Ind' or 'Cor' 
test = 'm';           % change 'test' to see effect of GPSP to different factors

switch test
  case 'm',   test0 = linspace(0.1,1.5,8);   
  case 's',   test0 = 2:10;  
  case 'r',   test0 = 0.02:0.02:0.2;  
  case 'v',   test0 = 0.1:0.1:0.9; Type = 'Cor';
  case 'n',   test0 = (5:5:20)*1e3;
end

f      = isequal(test,'n');  
S      = 100*(1-f)+10*f;
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
ylab = {'SNR','HD','HE','Time'};
figure('Renderer', 'painters', 'Position', [900, 200, 500 400])
set(0,'DefaultAxesTitleFontWeight','normal');
for j  = 1:4
    subplot(2,2,j)
    tmp = recd(:,j);
    plot(test0,tmp,'r*-','LineWidth',1), hold on,
    grid on, xlabel(test), title(ylab{j})  
    axis([min(test0) max(test0) min(tmp)/1.1 max(tmp)*1.1])   
end   
 