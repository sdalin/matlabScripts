%Find optimal k for kmeans clustering, do PCA, then re-do optimal k, for
%all data and then for each drug-resistance group individually

function [eva] = kmeansAnalysis(drugsCleaned,folder,drugType,varargin)
    %% Evaluate number of clusters
    eva = evalclusters(drugsCleaned,'kmeans','gap','KList',[1:size(drugsCleaned,1)],'Distance','cityblock');

    distanceErrors = exp(eva.LogW).*eva.SE;

    figure(1)
    ax = errorbar(eva.InspectedK,exp(eva.LogW),distanceErrors,'k','LineWidth',2);
    grid on
    xlim([0 50])
    xtick = [1,5,10,15,20,25,30,35,40,45,50];
    xlabel('Number of Clusters','FontSize',14)
    ylabel('Cityblock Distance Within Clusters','FontSize',14)
    title('Evaluation of Number of Clusters, cityblock distance','FontSize',14)
    
    %% Save the Plot
    cd(sprintf('%smatlabOutput/PCA',folder))
    %set(gcf,'paperpositionmode','manual')
    %set(gcf,'PaperUnits','inches')
    %set(gcf,'paperposition',[-0.30 -0.25 9 12])
    set(gcf,'Renderer','OpenGL')
    
    if nargin < 4
        if ~isempty(drugType)
            print(sprintf('ElbowPlot_%s',drugType),'-dpdf');
        else
            print('ElbowPlotAllData','-dpdf');
        end
    else
        if ~isempty(drugType)
            print(sprintf('ElbowPlot_%s_%s',varargin{1},drugType),'-dpdf');
        else
            print('ElbowPlot_%s_AllData',varargin{1},'-dpdf');
        end
    end
    
    cd(sprintf('%s../',folder))
    
    clf
    close all hidden
    
    

end