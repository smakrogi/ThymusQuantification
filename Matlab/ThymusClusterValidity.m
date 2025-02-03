function [result] = ThymusClusterValidity(DM, nClusters)
% Cluster validity experiment.

addpath('/home/makrogianniss/Software/sokrepo/m-files/UsefulFunctions-Matlab/FuzzyClusteringToolbox_m/FUZZCLUST');
colors={'r.' 'gx' 'b+' 'ys' 'md' 'cv' 'k.' 'r*' 'g*' 'b*' 'y*' 'm*' 'c*' 'k*' };

%the data
data.X=double(DM);
%normalization
%data=clust_normalize(data,'range');

%parameters
param.c=nClusters;
param.m=2;
param.e=1e-3;
param.val=1;

%FCM clustering
result = FCMclust(data,param);

%validation
result = validity(result,data,param);

if ( size(DM,2)==2 )
    figure,
    plot(data.X(:,1),data.X(:,2),'b.',result.cluster.v(:,1),result.cluster.v(:,2),'ro');
    hold on
    plot(result.cluster.v(:,1),result.cluster.v(:,2),'ro');
    xlabel('Fat-suppressed'), ylabel('Water-suppressed')
elseif ( size(DM,2)==3 )
    figure,
    subplot(131),
    plot(data.X(:,1),data.X(:,2),'b.',result.cluster.v(:,1),result.cluster.v(:,2),'ro');
    hold on
    plot(result.cluster.v(:,1),result.cluster.v(:,2),'ro');
    axis square;
    xlabel('Non-suppressed'), ylabel('Fat-suppressed')
    subplot(132),
    plot(data.X(:,1),data.X(:,3),'b.',result.cluster.v(:,1),result.cluster.v(:,3),'ro');
    hold on
    plot(result.cluster.v(:,1),result.cluster.v(:,3),'ro');
    axis square;
    xlabel('Non-suppressed'), ylabel('Water-suppressed')
    subplot(133),
    plot(data.X(:,2),data.X(:,3),'b.',result.cluster.v(:,2),result.cluster.v(:,3),'ro');
    hold on
    plot(result.cluster.v(:,2),result.cluster.v(:,3),'ro');
    axis square;
    xlabel('Fat-suppressed'), ylabel('Water-suppressed')
end

%evaluation
new.X=data.X;
eval = clusteval(new,result,param);

param.val=2;
result = validity(result,data,param);
% result.validity

saveas(gcf, ['clustervalidity', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.png']);

end