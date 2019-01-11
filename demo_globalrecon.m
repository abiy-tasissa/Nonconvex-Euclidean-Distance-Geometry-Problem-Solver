% -----------------------------------------------------------------------------------
% A script that loads different datum (e.g sphere,horse,kitten,swiss roll,cities,cow)
% and samples a subset of points uniformly at random. The outputs of this script are 
% used in main script: alternating_completion_*.m
% -----------------------------------------------------------------------------------
% Rongjie Lai, Abiy Tasissa
% -----------------------------------------------------------------------------------
close all
% specify choice of underlying object. sphere = 1 , cow = 2 , swiss roll = 3 
% cities = 4, random set of points in R^3 = 5
object_choice = 1;
if object_choice == 1
% uniform sampled sphere
[pt, trg] = ReadOFF('./data/1k.off','1');
num_pt = size(pt,1);
end
if object_choice == 2
% uniformly sampled cow 
[pt, trg]= ReadObjShape('./data/cow.off','1');
num_pt = size(pt,1);
end
if object_choice == 3
%  load swiss roll data with Euclidean distance
load('./data/ptswiss.mat'); 
num_pt = size(pt,1);
end
if object_choice == 4
% load US cities data
load('./data/UScities.mat');
pt = spt;  
num_pt = size(pt,1);
end
if object_choice==5
% generate random points in R^3
pt = rand(128,3);
num_pt = size(pt,1);
end
% -------------------------------------------------------------------------
% sampling rate
% -------------------------------------------------------------------------
rate = 0.02;
% -------------------------------------------------------------------------
% generate the distance matrix Dist_{i,j} = d_{ij}^{2}
% -------------------------------------------------------------------------
P = pt';
Dist = bsxfun(@plus,dot(P,P,1)',dot(P,P,1))-2*(P'*P);
% set diagonals to zero (for numerical reasons)
for i = 1:num_pt
    Dist(i,i)=0.0;
end
Dist = Dist/max(max(Dist));
%--------------------------------------------------------------------------
% generate weight for random-missing distance 
% (Weight==1 means available distance)
%--------------------------------------------------------------------------
Weight= rand(num_pt,num_pt);
Weight(Weight>1-rate)=1;
Weight(Weight<1)=0;
Weight(Weight>0)=1;
for i=1:num_pt
Weight(i,i)= 1;
for j=i+1:num_pt
    Weight(i,j)= Weight(j,i);
end
end
% -------------------------------------------------------------------------
% global coordinate reconstruction
% -------------------------------------------------------------------------
% choose the sampling setting: 
% choice = 1 assumes exact partial information
% choice = 2 assumes noisy partial information
choice = 1;
% aug. lagrangian parameters(r,r1,r2), noise parameter(lamda)
if choice == 1
    opts.r = 1.0;
end
if choice == 2
    opts.lamda = 5;
end
% options for printing
opts.printenergy = 0;
opts.printerror = 0;
% estimate of rank
opts.rank = 10;
% maximum number of iterations
opts.maxit = 30;
% tolerance
opts.tol = 1e-5;
% options for BB gradient method with nonmonotone line search
% optional for line search algorithm
lsopts.maxit = 20;
lsopts.xtol = 1e-8;
lsopts.gtol = 1e-8; 
lsopts.ftol = 1e-10; 
lsopts.alpha  = 1e-3; 
lsopts.rho  = 1e-4; 
lsopts.sigma  = 0.1; 
lsopts.eta  = 0.8; 
% -------------------------------------------------------------------------
% choose script to run based on the sampling setting i.e choice
% -------------------------------------------------------------------------
if choice==1
[GCor, ipm, output] = alternating_completion(Dist,Weight,opts,lsopts);
end
if choice==2
for i = 1:nruns
tic;
[GCor, ipm, output] = alternating_completion_noisy(Dist,Weight,opts);
end
end
% % -----------------------------------------------------------------------
% % output error results, view and save the reconstructed image in a folder
% % of choice. The default folder is results.
% % -----------------------------------------------------------------------
view_results = 1;
if view_results==1
if object_choice == 3 
[VV, DD]=eigs(ipm,num_pt);
ipm_err= output.ReconError;
fig4 = figure(4);
cmap = jet(num_pt);
scatter3(GCor(:,1),GCor(:,2),GCor(:,3),20,cmap);
axis off;
axis equal;
grid off;
str1 = strcat(num2str(object_choice),'_');
str2 = strcat(num2str(choice),'_');
str3 = strcat(num2str(rate*100));
% saveas(gcf,strcat('results/',str1,str2,str3),'jpg');
% saveas(gcf,strcat('results/',str1,str2,str3),'fig');
elseif object_choice == 4 || object_choice == 5
[VV,DD]=eigs(ipm,num_pt);
ipm_err= output.ReconError;
fig4 = figure(4);
ViewPC(GCor);
axis off;
axis equal;
grid off;   
str1 = strcat(num2str(object_choice),'_');
str2 = strcat(num2str(choice),'_');
str3 = strcat(num2str(rate*100));
% saveas(gcf,strcat('results/',str1,str2,str3),'jpg');
% saveas(gcf,strcat('results/',str1,str2,str3),'fig');   
else
[VV, DD]=eigs(ipm,num_pt);
ipm_err= output.ReconError;
fig4 = figure(4);
ViewMesh(GCor,trg);
str1 = strcat(num2str(object_choice),'_');
str2 = strcat(num2str(choice),'_');
str3 = strcat(num2str(rate*100));
% saveas(gcf,strcat('results/',str1,str2,str3),'jpg');
% saveas(gcf,strcat('results/',str1,str2,str3),'fig');
end
end
