%Make a nice biplot

function makeBiplot(wcoeff,score,drugsFixed,cellLines,folder,drug,PCAType,idx,varargin)
    %% Make the plot
    clear xLabel yLabel Title
    figure(2)
    hbi = biplot(wcoeff,'Scores',score,'VarLabels',drugsFixed,'ObsLabels',cellLines,'MarkerSize',15);
    xLabel = get(gca,'xlabel');
    xAxis = get(gca,'xaxis');
    yLabel = get(gca,'ylabel');
    yAxis = get(gca,'yaxis');
    Title = get(gca,'title');
    Title.String = 'PCA Biplot';
    Title.FontSize = 20;
    xAxis.FontSize = 15;
    yAxis.FontSize = 15;
    zAxis.FontSize = 15;
    xLabel.FontSize = 20;
    yLabel.FontSize = 20;
    zLabel.FontSize = 20;
    
    if ~isempty(varargin{1})
        yLabel.String = sprintf('Component %d',varargin{1});
    end
    
    
    
    %Change the colors and markers of the cellLines to correspond to which drug its
    %resistant to
    markerSpec = {'+','o','*','s','.','x','d','^','v','>','<','p','h'};
    cellLineTested = 0;
    for ii = 1:length(hbi)
        if strcmp(get(hbi(ii),'Tag'),'obsmarker')
            cellLineTested = cellLineTested + 1;
            %Set color based on drug type
            if strncmpi(cellLines(cellLineTested),'Dox',3)
                set(hbi(ii),'Color','r');
                set(hbi(ii),'MarkerSize',7);
            elseif strncmpi(cellLines(cellLineTested),'Vin',3)
                set(hbi(ii),'Color','m');
                set(hbi(ii),'MarkerSize',7);
            elseif strncmpi(cellLines(cellLineTested),'Pac',3)
                set(hbi(ii),'Color','g');
                set(hbi(ii),'MarkerSize',7);
            elseif strncmpi(cellLines(cellLineTested),'Cis',3)
                set(hbi(ii),'Color','c');
                set(hbi(ii),'MarkerSize',7);
            end
            
            %Set marker type based on k-means group if # groups less than
            %13
            if ~isnan(idx(cellLineTested))
                set(hbi(ii),'Marker',markerSpec{idx(cellLineTested)});
                set(hbi(ii),'MarkerSize',7);
            end
            
        end
    end
    
    
    %% Save the Plot
    cd(sprintf('%smatlabOutput/PCA',folder))
    %set(gcf,'paperpositionmode','manual')
    %set(gcf,'PaperUnits','inches')
    %set(gcf,'paperposition',[-0.30 -0.25 9 12])
    set(gcf,'Renderer','OpenGL')
    
    if ~isempty(drug)
        drug = drug{1};
    else
        drug = 'All Data';
    end
    
    if ~isempty(varargin{1})
        print(sprintf('%s_BiPlot_1%d_%s',PCAType,varargin{1},drug),'-dpdf');
    else
        print(sprintf('%s_BiPlot_%s',PCAType,drug),'-dpdf');
    end
    cd(sprintf('%s../',folder))
    
    clf
    close all hidden
end