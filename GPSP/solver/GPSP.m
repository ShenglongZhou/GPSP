function out = GPSP(A0,c,s,k,pars)
% One-bit compressed sensing problem is recovering sparse signal x from
%
%                c = sign(A0*x)
%
% This code aims at solving one-bit compressed sensing 
% via the double sparsity constrained optimization
%
%                min  ||Ax+y-eps||^2 + eta||x||^2
%                s.t. ||x||_0<=s, ||y_+||_0<=k
%
% where A = diag(c)*A0, eps>0, eta>0, s\in[1,n], k\in[0,m] are given.
% =========================================================================
% Inputs:
%     A0       : The sensing matrix \in R^{m-by-n},              (required)
%     c        : The binary observation \in R^m, c_i\in{-1,1}    (required)
%     s        : Sparsity level of x, an integer \in[1,n]        (required)      
%     k        : Upper bound of sign flips of (A0*x)             (required) 
%                An integer \in[1,m], e.g., k = ceil(0.01m)         
%     pars     : Parameters are all OPTIONAL
%                pars.eps   --  The parameter in the model         (default,1e-4)
%                pars.eta   --  The penalty parameter              (default,0.01/log(n))
%                pars.acc   --  Acceleration is used if acc=1      (default,1)
%                pars.maxit --  Maximum number of iterations       (default,1000) 
%                pars.tol   --  Tolerance of the halting condition (default,1e-9*sqrt(min(m,n)))
% Outputs:
%     out.obj  : Objective function value 
%     out.x    : The sparse solution in \R^n
%     out.y    : The  solution in \R^m
%     out.time : CPU time
%     out.iter : Number of iterations
% =========================================================================
% This code is written by Shenglong Zhou 
% It was programmed based on the algorithm proposed in 
%
% Shenglong Zhou, Ziyan Luo and Naihua Xiu, Computing one-bit 
% compressed sensing via double sparsity constrained optimization, 2021
%
% Send your comments and suggestions to <<< slzhou2021@163.com >>> 
% Warning: Accuracy may not be guaranteed !!!!! 
% =========================================================================

t0     = tic; 

if nargin<4, disp('Inputs are not enough, stop running ...'), return, end
if nargin<5, pars  = []; end

[m,n]  = size(A0);
if  n  <  1e4
    A = c.*A0;
else
    A = spdiags(c,0,m,m)*A0;    
end
Fnorm = @(var)norm(var,'fro')^2;
[maxit,tol,eta,eps,acc]          = GetParameters(m,n); 
if isfield(pars,'maxit'); maxit  = pars.maxit;   end
if isfield(pars,'tol');   tol    = pars.tol;     end
if isfield(pars,'eta');   eta    = pars.eta;     end
if isfield(pars,'eps');   eps    = pars.eps;     end
if isfield(pars,'acc');   acc    = pars.acc;     end

s0        = s;
if     s >= 0.01*n
       s  = ceil((1+m/n)*s); 
elseif s >= 0.005*n
       s  = ceil(1.2*s);  
end
 
x       = zeros(n,1); 
y       = zeros(m,1);
a       = 1;
barx    = x;
bary    = y;

T       = []; 
I       = [];
Axy     = y-eps;
obj     = Fnorm(Axy);  

HAM     = zeros(maxit,1); 
OBJ     = zeros(maxit,1);
stop0   = zeros(maxit,1);

fprintf('\n Start to run the solver: GPSP\n')
fprintf('--------------------------------------------\n');
fprintf(' Iter     Error       HamDist    Objective\n')
fprintf('--------------------------------------------\n');

for iter = 1:maxit
      
    alpha= 1; 
    T0   = T;
    I0   = I;
    x0   = x;
    y0   = y;
    obj0 = obj;
    v    = Axy; 
    u    = (v'*A)' + eta*barx;
 
    for j   = 1 : 15
        [x,T] = ProS(barx - alpha*u,s); 
        y     = ProK(bary - alpha*v,k); 
        AT    = A(:,T);
        Ax    = AT*x(T);
        Axy   = Ax - eps   + y;
        obj   = Fnorm(Axy) + eta*Fnorm(x(T));  
        gap   = Fnorm(x-barx) + Fnorm(y-bary);
        if obj  < obj0 - 1e-6*gap 
           break; 
        end
        alpha = alpha*0.25;
    end
         
    flag = isempty(setdiff(T,T0));    
    if flag
        I = find(y); 
        if  nnz(I) == nnz(I0)
            flag = isempty(setdiff(I,I0));
        end
    end
 
    if  flag  && nnz(I)<m
       if iter>2 && stop0(iter-2)==1 && stop0(iter-1)==1 
           break;
       end
       x1     = zeros(n,1);
       y1     = zeros(m,1);
       AIT    = AT(I,:);
       AT0    = AT(y==0,:);
       tmp    = (AT0'*AT0+eta*speye(s))\(eps*sum(AT0,1)'); 
       x1(T)  = tmp;  
       y1(I)  = eps-AIT*tmp;  
       Ax1    = AT*tmp;
       Axy1   = Ax1 - eps + y1;
       obj1   = Fnorm(Axy1) + eta*Fnorm(x1(T));
       gap    = Fnorm(x(T)-x1(T)) + Fnorm(y(I)-y1(I)); 
       stop0(iter) = 1; 
       if obj1    <= obj-1e-6*gap && nnz(y1>0)<=k   
          x        = x1;
          y        = y1;  
          Ax       = Ax1;
          Axy      = Axy1;
          obj      = obj1;              
      end           
   end
     
    sb        = sign(-c.*Ax); 
    sb(sb==0) = -1;
    ham       = 1-nnz(sb+c)/m;     
    HAM(iter) = ham; 
    OBJ(iter) = obj;   
    
    fprintf('%4d     %6.2e     %6.3f%%    %6.3e\n',...
             iter,gap,ham*100,obj);
    stop1 = (iter > 5 && gap < tol);     
    stop2 = (iter > 5 && std(HAM(iter-5:iter))<1e-6*log(n));  
    stop3 = (iter > 5 && std(OBJ(iter-5:iter))<1e-5*log(n));  
    stop4 = stop0(iter)*(n<1e4)+(n>=1e4);
    stop5 = (ham==1 && gap < 1e-4);
    if (stop1 && (stop2 || stop3) && stop4) || stop5
       break;
    end
 
    if mod(iter,50)==0, k = ceil(k/2);  end  
     
    if  acc % acceleration
        a0    = a;
        a     = (1 + sqrt(4*a0^2+1))/2; 
        barx  = x + ((a0-1)/a)*(x-x0);
        bary  = y + ((a0-1)/a)*(y-y0);
        if  stop0(iter)
            Ax  = AT*barx(T); 
        else
            T   = find(barx);
            Ax  = A(:,T)*barx(T);
        end
        barAxy  = Ax - eps  + bary;
        barobj  = Fnorm(barAxy) + eta*Fnorm(barx(T));
        
        if  barobj > obj
            barx  = x;
            bary  = y;  
            a     = a0;
        else
            Axy   = barAxy;
        end
    else % non-acceleration
        barx  = x;
        bary  = y; 
    end
    
end

fprintf('--------------------------------------------\n');  
if nnz(x)   > s0 
   [~,T]    = maxk(abs(x),s0);
   xn       = zeros(n,1);
   xn(T)    = x(T);
   x        = xn;
end
out.x    = x/norm(x);
out.y    = y;
out.obj  = obj;
out.OBJ  = OBJ(1:iter);
out.time = toc(t0);
out.iter = iter;
clear A b A0 B0 P
end

%--------------------------------------------------------------------------
function [maxit,tol,eta,eps,acc] = GetParameters(m,n)
    maxit = 1e3;
    tol   = 1e-9*sqrt(min(m,n));
    eta   = 1e-4; 
    eps   = 0.01*( (n<1e4) + (n>=1e4)/log(n) );
    acc   = 1;    
end

%--------------------------------------------------------------------------
function [xs,T] = ProS(x,s)
         [~,T]  = maxk(abs(x),s);  
         T      = sort(T);
         xs     = 0*x;
         xs(T)  = x(T);
end

%--------------------------------------------------------------------------
function y  = ProK(y,k)
    if   k  > 0
         ys = maxk(y,k);  
         if ys(end) > 0
            y  = y.*( y<0 | y >= ys(end) ); 
         end
    else
         y(y>0)=0; 
    end   
end
