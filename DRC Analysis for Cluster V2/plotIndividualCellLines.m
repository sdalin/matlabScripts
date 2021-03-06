%This function will make a hill plot with input from my hillFitv2
%script.

%INPUT: a fit object, the concentrations in a vertical matrix, the
%viabilities in a vertical matrix, the drug name, the cell line name, and
%the experiment name, and the folder to save the graphs into

%OUTPUT: none, although it saves pdfs of graphs into a folder.

function plotIndividualCellLines(fittedHill,fittedHillNoOutliers,concs,viability,outliers,rsquare,drug,cellLine,experiment,folder,quality)
    clf
    
    %Make the actual plot
    
    %This next line will hide the figures which will probably speeds up
    %running, but since the loop in hillFitV2 depends on the user seeing
    %the plot, jits important to comment it out for now.
    p = figure('visible','off');
    
    if ~isempty(fittedHill)
        p = plot(fittedHill,'r-',concs,viability,'k.',outliers,'m*'); 
        hold on
        plot(fittedHillNoOutliers,'c--')
        set(p,'MarkerSize',25)
    else
        p = scatter(concs,viability);
    end
    

    %Change the way it looks so its easily readable
    set(p,'LineWidth',2);  
    axis([-inf,inf,-0.2,1.2])
    set(gca,'XScale','log');
    xlabel('Drug Concentration (uM)','FontSize',13);
    ylabel('Viability %','FontSize',13);
    legend('off');
    title(sprintf('%s %s %s',drug,cellLine,experiment),'FontSize',15);
    
    %Label data points with their concentration
    textConcs = cellstr(num2str(concs));
    dy = 0.01;
    text(concs, viability+dy, textConcs);
    
    %Add Rsquare value to the plot
    text(10e-5,0,sprintf('individual rsquare is %f',rsquare));
    
    %Save into folder with raw data
    
    if quality == 0
        cd(sprintf('%s/matlabOutput/rawDRCs',folder))
    elseif quality == 1
        cd(sprintf('%s/matlabOutput/rejectedFits',folder))
    else
        'there is a problem in plotIndividualCellLines line 50'
    end
    
    set(gcf,'paperpositionmode','auto')
    %set(gcf,'Renderer','OpenGL')
    print(sprintf('%s %s %s.pdf',drug,cellLine,experiment),'-dpdf');%,'-bestfit');
    cd(sprintf('%s../',folder))
end