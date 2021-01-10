function [X,yf,y,xopt]=random1bcs(type,m,n,s,r,nf,v)
% This file aims at generating data of 2 types of Guassian samples 
% for 1-bit compressed sensing
% Inputs:
%       type  -- can be 'Ind','Cor' 
%       m     -- number of samples
%       n     -- number of features
%       s     -- sparsity of the true singnal, e.g., s = 0.01n
%       r     -- flipping ratio, 0~1,          e.g., r = 0.01
%       nf    -- nosie factor,                 e.g.,nf = 0.05
%       v     -- corrolation factor, 0~1,      e.g., v = 0.5
% Outputs:
%       X     --  samples data, m-by-n matrix
%       xopt  --  n-by-1 vector, i.e., the ture singnal
%       y     --  m-by-1 vector, i.e., sign(X*xopt+noise)
%       yf    --  m-by-1 vector, y after flapping some signs
%
% written by Shenglong Zhou, 19/07/2020

if     nargin <= 3
       error('Inputs are ont enough');
elseif nargin <= 4
       v = 0.5; nf = 0.05;  r = 0; 
elseif nargin <= 5 
       v = 0.5; nf = 0.05; 
elseif nargin <= 6 
       v = 0.5;
end
    
    
switch type
    case 'Ind'
          X = randn(m,n);
          
    case 'Cor'
          S = v.^(abs((1:n)-(1:n)'));
          X = mvnrnd(zeros(n,1),S,m); 
end
 
[xopt,T] = sparse(n,s);
y        = sign(X(:,T)*xopt(T)+nf*randn(m,1));
yf       = flip(y,r);
fprintf(' Done generation of sample data with:\n')
fprintf(' 1) Sample size: %d x %d\n',m,n)
fprintf(' 2) Sparsity level: %d\n',s)
fprintf(' 3) Sign flipping ratio: %4.2f\n',r)
fprintf(' 4) Noise ratio: %4.2f\n',nf)
end

% generate a sparse vector ------------------------------------------------
function [x,T] = sparse(n,s)
          I    = randperm(n);
          x    = zeros(n,1);
          T    = I(1:s);
          x(T) = randn(s,1);  
          x(T) = x(T) + sign(x(T));
          x(T) = x(T)/norm(x(T));
end

% flip the signs of a vector ----------------------------------------------
function yf   = flip(yopt,r)
         yf   = yopt;
         m    = length(yopt);
         I    = randperm(m);
         T    = I(1:ceil(r*m));
         yf(T)= -yf(T);
end
