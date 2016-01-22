%1/13/16
%Written by Simona Dalin
%K-means and heierarchical clustering of log2 IC50 data
%inputs: data (columns=drugs,rows=cell lines)

% cidx = nan(10,1);
% cmeans = nan(10,1);
% sumd = nan(10,1);
% 
% for numClust=1:10
%     [cidx(numClust),cmeans(numClust),sumd(numClust)] = kmeans(IC50data,numClust,'replicates',5);
% end

eucD = pdist(IC50data,'euclidean');
clustTreeEuc = linkage(eucD,'average');

cophenet(clustTreeEuc,eucD)

[h,nodes] = dendrogram(clustTreeEuc,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];

CGobj = clustergram(IC50data,'RowLabels',RowLabels,'ColumnLabels',ColumnLabels,'Colormap',redbluecmap);