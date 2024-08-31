close all
clear all
rng(10);

%algorithm application and, joint density and marginal distribution figures

% loading the data and preparing the y, x and u vectors/matrices
load ('C:\LMU\1-2\Seminar\Paper\Additional Analysis\Additional_Analysis\Data\FredMDlargeHor1.mat')
data=[X Y]; clear X Y;
data=[ones(size(data,1),1) zscore(data,1)]; 
y=data(:,end);
u=data(:,1);
x=data(:,2:end-1);

% options
abeta=1;
bbeta=1;
Abeta=1;
Bbeta=1;
M=110000;
N=10000;

% posterior inference
[store_B,store_z,store_phi,store_q,store_R2,store_gam,store_s2,y,x,u,T,k]=SpikeSlabGLP_Lap8(y,x,u,abeta,bbeta,Abeta,Bbeta,M,N);


% saving the results
save ('C:\LMU\1-2\Seminar\Paper\Additional Analysis\Additional_Analysis\Laplace\Estimation Results\Macro1_PosteriorDraws')

clear all
close all

%figure('Position', [0, 400, 800, 800]); 

cd ('C:\LMU\1-2\Seminar\Paper\Additional Analysis\Additional_Analysis\Laplace\Estimation Results'); 

% joint posterior of q and log(gamma)

% macro 1
load Macro1_PosteriorDraws.mat


if exist('store_z')==0;
    store_z=Z;
    store_q=Q';
    store_gam=sqrt(Gamma2)';
    intq=qq;
    N=0;
end;
Z=store_z(:,N+1:M)';
meanZ=mean(Z);

gridg = -4:0.005:2;
gridq = 0:.01:1;

[Gridg, Gridq] = meshgrid(gridg, gridq);
Grid                = [Gridg(:) Gridq(:)];

Temp = (2*gridq'*k*exp(2*gridg))./((gridq'*k*exp(2*gridg)+1).^2);
Fprior = reshape(Temp, size(Gridg));


f = ksdensity([log(store_gam(N+1:end)),store_q(N+1:end)], Grid);
        Fpost = reshape(f, size(Gridg));

subplot(1,2,1)        
contour(Gridg, Gridq, Fpost, 10,'-','LineWidth',2);
hold on
contour(Gridg, Gridq, Fprior,10,':','LineWidth',1);
hold off        
colormap(bone)
xlabel('log(\gamma)'); ylabel('q')
title('Macro 1')
legend('Posterior','Prior')


clear all

cd ('C:\LMU\1-2\Seminar\Paper\Additional Analysis\Additional_Analysis\Laplace\Estimation Results'); 

% posterior of q

% macro 1
load Macro1_PosteriorDraws.mat
if exist('store_z')==0;
    store_z=Z;
    store_q=Q';
    store_gam=sqrt(Gamma2)';
    intq=qq;
    N=0;
end;
edges=0:.01:1;
subplot(1,2,2)  
histogram(store_q(N+1:end),edges,'Normalization','pdf'); 
xlabel('q')
title('Macro 1')
xlim([0 1])
clear all

clear all