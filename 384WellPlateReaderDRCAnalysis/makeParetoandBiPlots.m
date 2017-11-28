%Make pareto plot and biplots for PCA

function makeParetoandBiPlots(explained,dataType,PCAType,folder,wcoeff,score,drugsCleaned,cellLines,idx)
    %Pareto plot of variability explained
    clear xLabel yLabel Title
    figure(1)
    pareto(explained)
    xLabel = xlabel('Principal Component');
    xLabel.FontSize = 15;
    yLabel = ylabel('Percent Variability Explained');
    yLabel.FontSize = 15;
    
    if ~strcmp(dataType,'AllData')
        Title = title(sprintf('%s (%s)',PCAType,dataType));
    else
        Title = title(sprintf('%s (All Data)',PCAType));
    end
    Title.FontSize = 15;
    
    cd(sprintf('%smatlabOutput/PCA',folder))
    set(gcf,'paperpositionmode','manual')
    set(gcf,'PaperUnits','inches')
    set(gcf,'paperposition',[-0.30 -0.25 9 12])
    set(gcf,'Renderer','OpenGL')
    
    if ~strcmp(dataType,'AllData')
        print(sprintf('%S ParetoPlot_%s',PCAType,dataType),'-dpdf');
    else
        print(sprintf('%s ParetoPlot',PCAType),'-dpdf');
    end
    cd(sprintf('%s../',folder))
    
    clf
    close all hidden
    

    %Biplot of scores and loadings
    components = {[1,2],[1,3],[1:3]};
    for i = 1:size(components,2)
        if i < 3
            yLabel = i + 1;
        else 
            yLabel = [];
        end
        
        makeBiplot(wcoeff(:,components{i}),score(:,components{i}),drugsCleaned,cellLines,folder,dataType,PCAType,idx,yLabel)
    end
end

