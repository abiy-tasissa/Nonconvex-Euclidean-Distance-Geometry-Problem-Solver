% -------------------------------------------------------------------------
% This code implements the Barzilai-Borwein(BB) gradient method.
% The original paper of Barzilai and Borwein is referenced below.
% -------------------------------------------------------------------------
% Ref 1. Barzilai, Jonathan, and Jonathan M. Borwein."Two-point step size 
% gradient methods." IMA journal of numerical analysis 8.1 (1988): 141-148.
% -------------------------------------------------------------------------
% Note that in the original BB method, the theoretical analysis for the
% convergence is for the case where the function in consideration is 
% quadratic. For a general non-quadratic case, the BB method is combined
% with nonmonotone line search. The code here implements the linesearch
% idea proposed in the paper below.
% -------------------------------------------------------------------------
% Ref 2. Zhang, Hongchao, and William W. Hager. "A nonmonotone line search 
% technique and its application to unconstrained optimization." 
% SIAM journal on Optimization 14.4 (2004): 1043-1056.
% -------------------------------------------------------------------------
% Rongjie Lai, Abiy Tasissa 
% -------------------------------------------------------------------------
function [x]= BBGradient(x, fun, opts)
% -------------------------------------------------------------------------
% Initialize the function f 
[f,g] = feval(fun, x);   
% norm of the gradient of f
gnorm  = norm(g, 'fro');
% Initialize Q , C = f(X_0) and alpha
Q = 1; C = f; alpha = opts.alpha;
% Initialize xmin, fmin
xmin = x; 
fmin = f;
% number of rows of x
n = size(x,1);
% maximum and minimum bounds for alpha
MAXalpha = 1e16;
MINalpha = 1e-16;
for it = 1: opts.maxit
    % Initialize and set parameters for line search
    xold = x ;
    fold = f;
    gold = g;
    numlinesearch = 1;
    wolfe_factor = opts.rho*gnorm^2;
    % Nonmonotone line search algorithm
    while 1
        x =  xold - alpha*gold;
        [f,g] = feval(fun, x);   
        if f <= C - alpha*wolfe_factor || numlinesearch >= 5
            break
        end
        alpha = opts.sigma*alpha;
        numlinesearch = numlinesearch+1; 
    end
    % Check if better minimum is obtained
     if f < fmin
        xmin = x;
        fmin = f;  
     end
    % Stopping criterion parameters
    gnorm  = norm(g, 'fro');
    s = x - xold;
    xstop = norm(s,'fro')/sqrt(n);
    fstop =  abs(fold-f)/(abs(fold)+1);
    % Stop main algorithm if stopping criterion are satisfied
    if ( xstop < opts.xtol && fstop < opts.ftol ) || (gnorm < opts.gtol)
       break;
    end
    % BB updates(Look up BB method for details)
    y = g - gold;   
    normy = norm(y,'fro');  
    sy = abs(trace(s'*y));
    if mod(it,2)==0 
        alpha = norm(s,'fro')^2/sy;
    else
        alpha = sy/(normy^2); 
    end
    % Ensure that alpha does not get too big or too small
    alpha = min(alpha,MAXalpha);
    alpha=  max(alpha,MINalpha);
    % Update cost C : This update is based on the idea in Ref 2. (See the
    % beginning of the code to for details)
    Qold = Q; 
    Q = opts.eta*Qold + 1; 
    C = (opts.eta*Qold*C + f)/Q;
end
% Update x
x = xmin;
end


