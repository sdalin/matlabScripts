%1/13/16
%Written by Simona Dalin
%K-means and heierarchical clustering of log2 IC50 data
%inputs: data (columns=drugs,rows=cell lines)
%%
%K-means Clustering

%Determine optimal number of clusters
eva = evalclusters(IC50data,'kmeans','gap','KList',[1:32],'Distance','cityblock');

distanceErrors = exp(eva.LogW).*eva.SE;

figure(1)
ax = errorbar(eva.InspectedK,exp(eva.LogW),distanceErrors,'k','LineWidth',2);
grid on
set(gca,'XTick',eva.InspectedK)
xlim([0 32])
xlabel('Number of Clusters')
ylabel('Euclidian Distance Within Clusters')
title('Evaluation of Number of Clusters, K-means, cityblock distance')

idx9 = kmeans(IC50data,9);


%%
%Using linkage then cluster generates this:
T = clusterdata(IC50data,1);

%%
%Heiraricical Clustering
eucD = pdist(data,'euclidean');
clustTreeEuc = linkage(eucD,'average');

cophenet(clustTreeEuc,eucD);

[h,nodes] = dendrogram(clustTreeEuc,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];

CGobj = clustergram(data,'RowLabels',RowLabels,'ColumnLabels',ColumnLabels,'Colormap',redbluecmap,'Symmetric',true,'ColumnLabelsRotate',0);


%%
%Some dimensional reduction to get a look at these suckers

[wcoeff,score,latent,tsquared,explained] = pca(IC50data);

