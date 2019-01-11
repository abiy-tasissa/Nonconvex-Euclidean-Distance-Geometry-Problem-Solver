% --------------------------------------------------------------------------
% This script is a numerical solution of the Euclidean distance geometry
% problem (EDG). Formally, consider a set of n points where partial 
% inter-point distance information is provided. The goal of EDG is to find
% the coordinate of the points given this partial information.
% --------------------------------------------------------------------------
% The n points usually lie in a low dimensional space of size r << n, low
% rank. With this, the EDG problem can be set as low-rank completion problem
% which can be solved via nuclear norm minimization. We recover the Gram 
% matrix, the inner product matrix, and follow classical MDS to recover the
% coordinates. For details, see the associated paper below.
% --------------------------------------------------------------------------
% Tasissa, Abiy, and Rongjie Lai."Exact Reconstruction of Euclidean Distance 
% Geometry Problem Using Low-rank Matrix Completion." arXiv preprint 
% arXiv:1804.04310 (2018).
% --------------------------------------------------------------------------
% The gram matrix is psd so nuclear norm minimization equates to trace.
% Our algorithm uses the Augmented Lagrangian framework to find the Gram 
% matrix. 
% --------------------------------------------------------------------------
% We assume the partial inter-distance information comes from a uniform
% random sample. The information is assumed to be exact.
% Dist is the full distance matrix. Weight is a binary matrix informing 
% whether a given entry of Dist is chosen or not.
% --------------------------------------------------------------------------
% Rongjie Lai, Abiy Tasissa
% --------------------------------------------------------------------------
function [Global_Coordinate, IPM_Recon, output]=alternating_completion(Dist,Weight,opts,lsopts)
% aug. lagrangian penalty and estimate of the rank
r = opts.r;
Rk = opts.rank;
% dim = number of points
% Dist = D^{2}(i,j) = {d_{i,j}^{2}} matrix
[dim , ~] = size(Dist);
% calculate the ground truth of inner-product matrix
IPM_Truth = Dist - mean(Dist,2)*ones(1,dim);
IPM_Truth =  - 1/2*(IPM_Truth - ones(dim,1)*mean(IPM_Truth,1));
% indices of the randomly chosen entries of D
[I,J] = find(Weight==1);
diagind = (1:dim + 1:dim*dim);
edgeind = I + (J - 1)*dim;
% the minimization problem has linear constraint R_{Omega}(X) = R_{Omega)(M)
% Represent this constraint as A(X) = b where b = R_{omega}(M)
M = Dist(edgeind);
b = M;
% -------------------------------------------------------------------------
% main algorithm
% -------------------------------------------------------------------------
% initialize P, lagrangian multiplier
P = rand(dim,Rk);
R = P;
D1 = zeros(length(b),1);
% initialize energies
E1 = zeros(opts.maxit,1);
E = zeros(opts.maxit,1);
num_it = 0;
cre = 1 ;
% -------------------------------------------------------------------------
% main iteration: BB method to solve for P
% -------------------------------------------------------------------------
for i = 1:opts.maxit
    num_it = num_it + 1;
    % do line search based gradient descent for P    
    P = BBGradient(P,@(P)gradient(P),lsopts);
    % update multiplier D
    tmperr = A_operator(P)-b;
    D1 = D1 + tmperr;
    % total energy
    E(i) =  sum(sum(P.*P,2)) + 0.5*r*norm(D1,'fro')^2;
    if opts.printenergy==1
    fprintf('Iteration %d, TotalE = %f\n',i,E(i));
    end
    % calculate the energy E1 
    E1(i) = 0.5 * r * norm(tmperr,'fro')^2;
    % stopping condition
    if i > 1
        cre = abs(E(i) - E(i-1))/E(i);
    end
    if(E1(i) < opts.tol && cre < opts.tol)
        break;
    end
end
% -------------------------------------------------------------------------
% plotting the different energies
% -------------------------------------------------------------------------
if opts.printerror==1
fig1=figure(1);
set(fig1,'defaulttextinterpreter','latex');
plot(1:num_it,log(E1(1:num_it)),'LineWidth',2);
xlabel('iteration number = k','FontSize',18)
ylabel('$\log(E_{1})$' ,'FontSize',18)
grid on
saveas(fig1,'results/cities_exact_sample_constraint_1','png')
fig2=figure(2);
set(fig2,'defaulttextinterpreter','latex');
plot(1:num_it,log(E(1:num_it)),'LineWidth',2);
xlabel('iteration number = k','FontSize',18);
ylabel('$\log(E)$','FontSize',18);
grid on
saveas(fig2,'results/cities_exact_total_energy_1','png')
end
% -------------------------------------------------------------------------
% get the gram matrix, inner-product matrix, IPM_Recon and follow canonical 
% Multidimensional scaling (MDS)
% -------------------------------------------------------------------------
IPM = P*P';
IPM_Recon = IPM-(1/dim)*repmat(sum(IPM,2),1,dim)-(1/dim)*repmat(sum(IPM,1),dim,1)+...
(1/(dim*dim))*repmat(sum(sum(IPM)),dim,dim);
IPM_Recon = (IPM_Recon + IPM_Recon')/2;
IPM_err = norm(IPM_Truth - IPM_Recon,'fro')/norm(IPM_Truth,'fro');
[V, D] = eigs(IPM_Recon,opts.rank,'lm');
D = diag(D);
% plot the eigenvalues
% figure;
% plot(D);
% title('Eigenvalues','FontSize',20);
[D,IJ] = sort(D,'descend');
V = V(:,IJ);
Global_Coordinate = real(V(:,1:3)*diag(sqrt(D(1:3))));
% -------------------------------------------------------------------------
% output parameters: constraint energy E1, total energy E, relative error 
% in the gram matrix
% -------------------------------------------------------------------------
output.E1 = E1(1:num_it);
output.E =  E(num_it);
output.ReconError = IPM_err;
output.numit = num_it;
% -------------------------------------------------------------------------
% Constructs the operator A which captures the linear operator: 
% R_{\omega}(X) = R_{\omega}(M)
% -------------------------------------------------------------------------
function [Y] = A_operator(X)
        % diagonals and off diagonals of X*X'
        X_diag = sum(X.*X,2);
        X_offdiag = sum(X(I,:).*X(J,:),2);
        % resttiction operator
        Y = X_diag(I)+X_diag(J)-2*X_offdiag;        
 end
% -------------------------------------------------------------------------
% Constructs the adjoint operator for A (see associated paper for details). 
% -------------------------------------------------------------------------
    function [X] = At_operator(y)
        % first part of A^{*}(y) : \sum_{alpha in omega} y^{1}_{alpha} w_{alpha}
        X = zeros(dim,dim);
        X(edgeind) = -2*y(1:end);
        X(diagind) =  -sum(X,2);
    end
% -------------------------------------------------------------------------
% a function handle for the line search BB algorithm
% -------------------------------------------------------------------------
    function [F,G] = gradient(P)
        tmp = A_operator(P)-b+D1;
        % F is the objective function
        F = sum(sum(P.*P,2))+0.5*r*norm(tmp,'fro')^2;
        % G is the gradient
        G = 2*P+2.0*r*At_operator(tmp)*P;
    end
end




