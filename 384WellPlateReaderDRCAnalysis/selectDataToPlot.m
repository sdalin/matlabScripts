function [selectedDrugs,selectedCellLines] = selectDataToPlot(dataAfterFit)
% SIMPLE_GUI2 Select drugs and cell lines from list menus. Clicking a
% button returns the drugs/cell lines selected (which will then be plotted
% with a different function
 
   clf

   drugs = fieldnames(dataAfterFit.rawData);
   %strip "drug_" from drug names
   drugsStripHead = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);
   cellLines = dataAfterFit.fitParams.drug_Doxorubicin.EC50(1,2:end);

   %  Create and then hide the GUI as it is being constructed.
   f = figure('Visible','off','Position',[360,500,450,585]);
 
   %  Construct the components.
   plotCombosButton = uicontrol('Style','pushbutton','String','Plot these drug/cell line combos','FontSize',15,...
          'Position',[108,37,251,26],...
          'Callback',@plotCombosButton_Callback); 
   drugText = uicontrol('Style','text','String','Drugs available for plotting','FontSize',20,...
          'Position',[100,550,271,26]);
   cellLineText = uicontrol('Style','text','String','Cell lines available for plotting','FontSize',20,...
          'Position',[88,276,283,26]);
   drugList = uicontrol('Style','listbox',...
          'max',length(drugsStripHead),...
          'min',1,...
          'String',drugsStripHead,...
          'Position',[17,346,423,197],...
          'Callback',@drugList_Callback);
   cellLineList = uicontrol('Style','listbox',...
          'max',length(drugsStripHead),...
          'min',1,...
          'String',cellLines,...
          'Position',[17,67,423,197],...
          'Callback',@cellLineList_Callback);
   
   % Initialize the GUI.
   % Change units to normalized so components resize 
   % automatically.
   f.Units = 'normalized';
   plotCombos.Units = 'normalized';
   drugText.Units = 'normalized';
   cellLineText.Units = 'normalized';
   drugList.Units = 'normalized';
   cellLineList.Units = 'normalized';

   % Assign the GUI a name to appear in the window title.
   f.Name = 'Simple GUI';
   % Move the GUI to the center of the screen.
   movegui(f,'center')
   % Make the GUI visible.
   f.Visible = 'on';
  
   %List callbacks. Read the list menu Values property
   %to determine which items are currently displayed and make it
   %the current data.
   function drugList_Callback(source,eventdata) 
       % Determine the selected data set.
       selectedDrugs = drugs(source.Value);
   end

   function cellLineList_Callback(source,eventdata)
       % Determine the selected data set.
       if ~sum(ismember(cellLines(source.Value),'Parental')) %User didn't select parental cell line
           selectedCellLines = cellLines(source.Value);
           selectedCellLines(end+1) = {'Parental'};
       else
           selectedCellLines = cellLines(source.Value);
       end
   end
   
   % Push button callback. Each callback plots current_data in
   % the specified plot type.
 
   function plotCombosButton_Callback(source,eventdata) 
       % Return the relevant data and close the GUI
       display('Button pressed');
       uiresume(f);
       close(f);
   end

uiwait



end 