% %1/13/16
% %Written by Simona Dalin
% %K-means and heierarchical clustering of log2 IC50 data
% %inputs: data (columns=drugs,rows=cell lines)
% %%
% %K-means Clustering
% 
% %Determine optimal number of clusters
% eva = evalclusters(IC50data,'kmeans','gap','KList',[1:32],'Distance','cityblock');
% 
% distanceErrors = exp(eva.LogW).*eva.SE;
% 
% figure(1)
% ax = errorbar(eva.InspectedK,exp(eva.LogW),distanceErrors,'k','LineWidth',2);
% grid on
% set(gca,'XTick',eva.InspectedK)
% xlim([0 32])
% xlabel('Number of Clusters')
% ylabel('Euclidian Distance Within Clusters')
% title('Evaluation of Number of Clusters, K-means, cityblock distance')
% 
% idx9 = kmeans(IC50data,9);
% 
% 
% %%
% %Using linkage then cluster generates this:
% T = clusterdata(IC50data,1);
% 
% %%
%Heiraricical Clustering
eucD = pdist(IC50data,'euclidean');
clustTreeEuc = linkage(eucD,'average');

cophenet(clustTreeEuc,eucD);

[h,nodes] = dendrogram(clustTreeEuc,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];

%CGobj = clustergram(IC50data,'RowLabels',RowLabels,'ColumnLabels',ColumnLabels,'Colormap',redbluecmap,'Symmetric',true,'ColumnLabelsRotate',0,'Cluster','row');
% 
% 
%%
%Some dimensional reduction to get a look at these suckers

[wcoeff,score,latent,tsquared,explained] = pca(IC50data);
figure(1)
biplot(wcoeff(:,1:3),'Scores',score(:,1:3),'VarLabels',ColumnLabels,'ObsLabels',RowLabels,'MarkerSize',15)

%%
%Testing classical MDS
%Make distance matrix of IC50data
IC50dist = squareform(pdist(IC50data,'cityblock'));
%Find matrix Q:
Q = -0.5*IC50dist.^2;
%Find the centering matrix H:
n = 32;
H = eye(n) - ones(n)/n;
%find the matrix B:
B = H*Q*H;
%Get the eigenvalues of B.
[A,L] = eig(B);
%Find the ordering largest to smallest.
[vals, inds] = sort(diag(L));
inds = flipud(inds);
vals = flipud(vals);
%re-sort based on these.
A = A(:,inds);
L = diag(vals);

%First plot scree-type plot to look for the Elbow.
%The following uses a log scale on the y-axis.
figure(4)
semilogy(vals(1:5),'o');
%Using 2-D for visualization purposes,
%find the coordinates in teh lower-dimensional space.
X = A(:,1:3)*diag(sqrt(vals(1:3)));
%Now plot in a 2-D scatterplot
%nice colors
%c = linspace(1,20,length(X));
%offset labels with index of each vector
dx = 0.1;
dy = 0.1;
dz = 0.1;
a = SelectedRowLabels;
labels = cellstr(a);
figure(5)
scatter3(X(:,1),X(:,2),X(:,3),chosenSizes,putativeClusters,'filled')
text(X(:,1)+dx,X(:,2)+dy,X(:,3)+dz,labels)

%%
%Calculate differences via euclidian distance, correlation and spearman
%correlation.  Find pair with smallest difference -> this is one pair to
%investigate.  Then find other cell line with largest difference from the
%pair.  Look for cell line with smallest difference to this new cell line.
%This is the second pair to investigate.  Etc.

distanceColumnLabels ={'First most similar cell line','Second most similar cell line','Cell line most distant to first cell line','Cell line most distant from second cell line','Cell line most similar to first distant cell line','Cell line most similar to second distant cell line'};

distanceMetric = {'euclidean','seuclidean','correlation','cityblock'};

for metric = 1:length(distanceMetric)
    dist = pdist(IC50data,sprintf('%s',distanceMetric{metric}));
    squareDist = squareform(dist);
    
    %Find most similar two cell lines
    minDist = min(dist);
    idxMinDist = find(squareDist==minDist);
    [rowMinDist,colMinDist]=ind2sub(size(squareDist),min(idxMinDist));

    distantCellLines{metric,1} = RowLabels{rowMinDist};
    distantCellLines{metric,2} = RowLabels{colMinDist};

    %Find most distant cell lines to each of those two
    mostDistant1 = find(squareDist(rowMinDist,:) == max(squareDist(rowMinDist,:)));
    mostDistant2 = find(squareDist(colMinDist,:) == max(squareDist(colMinDist,:))); 

    distantCellLines{metric,3} = RowLabels{mostDistant1};
    distantCellLines{metric,4} = RowLabels{mostDistant2};

    %Find most similar cell lines to each of the most distant cell lines
    secondSimilar1 = find(squareDist(mostDistant1,:) == min(nonzeros((squareDist(mostDistant1,:)))));
    secondSimilar2 = find(squareDist(mostDistant2,:) == min(nonzeros((squareDist(mostDistant2,:)))));

    distantCellLines{metric,5} = RowLabels{secondSimilar1};
    distantCellLines{metric,6} = RowLabels{secondSimilar2};
    
    %Find cell line most distant to both of the above cell line pairs
    combinedDistance1 = squareDist(rowMinDist,:) + squareDist(mostDistant1,:);
    combinedDistance2 = squareDist(rowMinDist,:) + squareDist(mostDistant2,:);
    combinedDistance3 = squareDist(colMinDist,:) + squareDist(mostDistant1,:);
    combinedDistance4 = squareDist(colMinDist,:) + squareDist(mostDistant2,:);
    
    bothDist1 = find((combinedDistance1) == max(combinedDistance1));
    bothDist2 = find((combinedDistance2) == max(combinedDistance2));
    bothDist3 = find((combinedDistance3) == max(combinedDistance3));
    bothDist4 = find((combinedDistance4) == max(combinedDistance4));
    
    distantCellLines{metric,7} = RowLabels{bothDist1};
    distantCellLines{metric,8} = RowLabels{bothDist2};
    distantCellLines{metric,9} = RowLabels{bothDist3};
    distantCellLines{metric,10} = RowLabels{bothDist4};
    
    %Find cell line closest to above four 'both distant' cell lines
    thirdSimilar1 = find(squareDist(bothDist1,:) == min(nonzeros((squareDist(bothDist1,:)))));
    thirdSimilar2 = find(squareDist(bothDist2,:) == min(nonzeros((squareDist(bothDist2,:)))));
    thirdSimilar3 = find(squareDist(bothDist3,:) == min(nonzeros((squareDist(bothDist3,:)))));
    thirdSimilar4 = find(squareDist(bothDist4,:) == min(nonzeros((squareDist(bothDist4,:)))));

    distantCellLines{metric,11} = RowLabels{thirdSimilar1};
    distantCellLines{metric,12} = RowLabels{thirdSimilar2};
    distantCellLines{metric,13} = RowLabels{thirdSimilar3};
    distantCellLines{metric,14} = RowLabels{thirdSimilar4};

end
    

%%
%SOM

D = som_normalize(IC50data,'var');



