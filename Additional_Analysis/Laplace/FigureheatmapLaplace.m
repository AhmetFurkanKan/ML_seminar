clear all
close all

cd ('C:\LMU\1-2\Seminar\Paper\Additional Analysis\Additional_Analysis\Laplace\Estimation Results'); % This is the folder containing the posterior draws used in the paper

% heatmap
figure('Position', [0, 400, 800, 800]); 

% macro 1
load Macro1_PosteriorDraws110k.mat
if exist('store_z')==0;
    store_z=Z;
    store_q=Q';
    store_gam=sqrt(Gamma2)';
    intq=qq;
    N=0;
end;
Z=store_z(:,N+1:M)';
meanZ=mean(Z);
ax1=subplot(1,1,[1]);
imagesc(meanZ)
colorbar
caxis([0, 1])
colormap(ax1,flipud(colormap(ax1,'hot')))
set(gca,'YTickLabel','','Ytick',[])
xlabel('coefficients')
title('Macro 1')
clear all
