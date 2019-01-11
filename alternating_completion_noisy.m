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
% random sample with additive Gaussian noise. 
% Dist is the full distance matrix. Weight is a binary matrix informing 
% whether a given entry of Dist is chosen or not.
% --------------------------------------------------------------------------
% Rongjie Lai, Abiy Tasissa
% --------------------------------------------------------------------------
function [Global_Coordinate, IPM_Recon, output]=alternating_completion_noisy...
(Dist,Weight,opts,lsopts)
% aug. lagrangian penalty, noise parameter, and estimate of the rank
lamda = opts.lamda;
Rk = opts.rank;
% dim = number of points
% Dist_squ = D^{2}(i,j) = {d_{i,j}^{2}} 
[dim,~] = size(Dist);
% calculate the ground truth of inner-product matrix 
IPM_Truth = Dist - mean(Dist,2)*ones(1,dim);
IPM_Truth =  - 1/2*(IPM_Truth - ones(dim,1)*mean(IPM_Truth,1));
% indices of the randomly chosen entries of D
% note that, edgeind contains both symmetric indices
% i.e (1,2) and (2,1)
[I,J] = find(Weight==1);
diagind = (1:dim + 1:dim*dim);
edgeind = I + (J - 1)*dim;
% the minimization problem has linear constraint R_{Omega}(X) = R_{Omega)(M)
% now perturbed by noise
M_noisefree = Dist(edgeind);
% minimium distance
mindist = min(M_noisefree(M_noisefree>0));
% add random noise, of order minimunm distance, and ensure
% resulting distances are still positive
M = M_noisefree + 3*mindist + mindist*randn(size(M_noisefree));
M(M_noisefree==0) = 0;
b = M;
% -------------------------------------------------------------------------
% main algorithm
% -------------------------------------------------------------------------
% initialize P, lagrangian multipliers
P = rand(dim,Rk);
% initialize energies
E = zeros(opts.maxit,1);
% initialize iteration counters and relative errors
num_it = 0;
cre = 1 ;
% -------------------------------------------------------------------------
% main iteration: BB method to solve for P 
% -------------------------------------------------------------------------
for i = 1:opts.maxit
    num_it = num_it + 1;
    % do line search based gradient descent for P
    P = BBGradient(P,@(P)gradient(P),lsopts);
    % total energy
    E(i) = sum(sum(P.*P,2))+0.5*lamda*norm(A_operator(P)-b,'fro')^2;
    if opts.printenergy==1
        fprintf('Iteration %d, TotalE = %f\n',i, E(i));
    end
    % stopping condition
    if i > 1
        cre = abs(E(i) - E(i-1))/E(i);
    end
    if(cre < 1e-5)
        break;
    end
end
% -------------------------------------------------------------------------
% plotting the different energies
% -------------------------------------------------------------------------
if opts.printerror==1
    fig1=figure(1);
    set(fig1,'defaulttextinterpreter','latex');
    plot(1:num_it,log(E(1:num_it)),'LineWidth',2);
    xlabel('iteration number = k','FontSize',18)
    ylabel('$\log(E)$','FontSize',18);
    grid on
    saveas(fig1,'results/cities_noisy_total_energy_5','png')
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
Y = X_diag(I)+X_diag(J)-2*X_offdiag;
end
% -------------------------------------------------------------------------
% Constructs the adjoint operator A* (see associated paper for details). 
% -------------------------------------------------------------------------
function [X] = At_operator(y)
% first part of A^{*}(y) : \sum_{i} y^{1}_{alpha_i} w_{alpha_i}
X = zeros(dim,dim);
X(edgeind) = -2*y(1:end);
X(diagind) =  -sum(X,2);
end
% -------------------------------------------------------------------------
% a function handle for the line search BB algorithm
% -------------------------------------------------------------------------
function [F,G] = gradient(P)
tmp1 = A_operator(P)-b;
% F is the objective function 
F = sum(sum(P.*P,2))+ 0.5*lamda*norm(tmp1,'fro')^2;
% G is the gradient
G = 2*P+2.0*lamda*At_operator(tmp1)*P;
end
end




