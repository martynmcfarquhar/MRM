%==========================================================================
% MRM Post-estimation Tools
%==========================================================================
% Function description
%--------------------------------------------------------------------------
% Define the GUI and functions for the post-estimation tools window
%
%=========================================================================%
% Copyright 2016 Martyn McFarquhar                                        %
%=========================================================================%
% This file is part of MRM.
%
% MRM is free software: you can redistribute it and/or modify it under the 
% terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or (at your option) any 
% later version.
%
% MRM is distributed in the hope that it will be useful, but WITHOUT ANY 
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details.
%
% You should have received a copy of the GNU General Public License along 
% with MRM.  If not, see <http://www.gnu.org/licenses/>.
%=========================================================================% 
function MRM_postEstimationTools(MRM)

screen = get(0, 'ScreenSize');
gui.postEstToolsWin = figure('Position', [1,1,0.8*screen(3),0.8*screen(4)], ...
                             'DockControls', 'off', 'MenuBar', 'none', 'Visible', 'off', ...
                             'Name', 'MRM Post-estimation tools', 'NumberTitle', 'off',  ...
                             'WindowButtonMotionFcn', @mouse_motion, ...
                             'WindowButtonUpFcn',     @mouse_release, ...
                             'ResizeFcn',             @win_resize,    ...
                             'CloseRequestFcn',       @win_close);
  
gui.MRM                  = MRM;
gui.dragging             = [];
gui.overlay              = 0;
gui.overlayInterp        = 0; % 0 = NN, 1 = Tril
guidata(gcf,gui)

create_widgets();
path = fileparts(which('MRM_launcher.m'));
load_background_img([path filesep 'Utilities' filesep 'Template' filesep 'template1.nii']);
position_widgets();
draw_images();
fill_results_list();
fill_plot_list();

set(gui.postEstToolsWin, 'Visible', 'on');

table_row_highlight();

gui = getappdata(gcf, 'UsedByGUIData_m');
gui.oldWinSize = [gui.postEstToolsWin.Position(3) gui.postEstToolsWin.Position(4)];
gui.oldTabSize = [gui.resultsTable.Position(3) gui.resultsTable.Position(4)];
guidata(gcf,gui);

end


%==========================================================================
% Create the widgets for the GUI
%==========================================================================
function create_widgets()

gui = getappdata(gcf, 'UsedByGUIData_m'); 

gui.coronalPlot  = axes('Color', [0,0,0], 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [], 'Box', 'on');
gui.axialPlot    = axes('Color', [0,0,0], 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [], 'Box', 'on');
gui.saggitalPlot = axes('Color', [0,0,0], 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [], 'Box', 'on');
gui.CBarPlot     = axes('Color', [0,0,0], 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [], 'Box', 'off');
                    
gui.infoTable               = uitable('ColumnName', {}, 'RowName', {});
gui.resultsTable            = uitable('ColumnName', {'Location', 'Cluster', 'Cluster size', 'Statistic','p-value', 'X', 'Y', 'Z'}, 'RowName', {}, 'CellSelectionCallback', @table_click);
gui.xCoordBox               = uicontrol('Style', 'edit', 'Callback', @coords_changed);
gui.yCoordBox               = uicontrol('Style', 'edit', 'Callback', @coords_changed);
gui.zCoordBox               = uicontrol('Style', 'edit', 'Callback', @coords_changed);
gui.xhairCheck              = uicontrol('Style', 'checkbox',   'String', 'Xhairs', 'Value', 1, 'Callback', @xhairs_off);
gui.backgroundMinBox        = uicontrol('Style', 'edit',       'String', '0',    'Callback', @change_background_contrast);
gui.backgroundMaxBox        = uicontrol('Style', 'edit',       'String', '0',    'Callback', @change_background_contrast);
gui.autoBackContrastButton  = uicontrol('Style', 'pushbutton', 'String', 'Auto', 'Callback', @auto_background_contrast);
gui.barPlotAxes             = axes('Box', 'on', 'ButtonDownFcn', @plot_clicked);

% Results tab
gui.resultsList         = uicontrol('Style', 'popupmenu',  'Callback', @selected_result_changed);
gui.nnInterpTick        = uicontrol('Style', 'checkbox',   'String', 'Nearest neighbour',   'Value', 1, 'Callback', @interp_change_callback);
gui.trilInterpTick      = uicontrol('Style', 'checkbox',   'String', 'Trilinear',           'Value', 0, 'Callback', @interp_change_callback);
gui.threshImgTick       = uicontrol('Style', 'checkbox',   'String', 'Thresholded image',   'Value', 1, 'Callback', @thresh_image_callback);
gui.unthreshImgTick     = uicontrol('Style', 'checkbox',   'String', 'Unthresholded image', 'Value', 0, 'Callback', @thresh_image_callback);
gui.threshBox           = uicontrol('Style', 'edit',       'String', '0',                   'Enable', 'off', 'Callback', @change_overlay_thresh);
gui.colourMapList       = uicontrol('Style', 'popupmenu',  'String', {'Hot', 'Cool', 'Jet', 'Spring', 'Summer', 'Autumn', 'Winter', 'Parula', 'HSV'}, 'Callback', @change_overlay_colourmap);
gui.overlayMinBox       = uicontrol('Style', 'edit',       'String', '0',    'Callback', @change_overlay_contrast);
gui.overlayMaxBox       = uicontrol('Style', 'edit',       'String', '0',    'Callback', @change_overlay_contrast);
gui.autoContrastButton  = uicontrol('Style', 'pushbutton', 'String', 'Auto', 'Callback', @auto_overlay_contrast);

% Plot tab
gui.plotConList         = uicontrol('Style', 'listbox',     'Callback', @selected_plot_changed);
gui.newConButton        = uicontrol('Style', 'pushbutton',  'String', 'New', 'Callback', @new_plot_con);
gui.viewConButton       = uicontrol('Style', 'pushbutton',  'String', 'View');
gui.deleteConButton     = uicontrol('Style', 'pushbutton',  'String', 'Delete');
gui.SECheck             = uicontrol('Style', 'checkbox',    'String', 'SE',                   'Value', 0, 'Callback', @CILevel_callback);
gui.CICheck             = uicontrol('Style', 'checkbox',    'String', 'CIs',                  'Value', 1, 'Callback', @CILevel_callback);
gui.CILevel             = uicontrol('Style', 'edit',        'String', '95', 'Callback', @CI_level_callback);
gui.plot2dTick          = uicontrol('Style', 'checkbox',    'String', 'Single',               'Value', 1,                  'Callback', @plotTypeCallback);
gui.plot3dTick          = uicontrol('Style', 'checkbox',    'String', 'Multiple',             'Value', 0,                  'Callback', @plotTypeCallback);
gui.neighbour6Radio     = uicontrol('Style', 'radiobutton', 'String', '6',                    'Value', 1, 'Enable', 'off', 'Callback', @neighbourCallback);
gui.neighbour18Radio    = uicontrol('Style', 'radiobutton', 'String', '18',                   'Value', 0, 'Enable', 'off', 'Callback', @neighbourCallback);
gui.neighbour26Radio    = uicontrol('Style', 'radiobutton', 'String', '26',                   'Value', 0, 'Enable', 'off', 'Callback', @neighbourCallback);
%gui.neighbourClustRadio = uicontrol('Style', 'radiobutton', 'String', 'Cluster',              'Value', 0, 'Enable', 'off', 'Callback', @neighbourCallback);
gui.savePlotValsButton  = uicontrol('Style', 'pushbutton',  'String', 'Save to text', 'Callback', @save_plot_text);
gui.plotValsWorkspace   = uicontrol('Style', 'pushbutton',  'String', 'Return in workspace', 'Callback', @return_plot_vals);

% Model tab
gui.paramEstsTable           = uitable('ColumnName', {}, 'RowName', {}, 'BackgroundColor', [1 1 1]);
gui.vcovTable                = uitable('ColumnName', {}, 'RowName', {}, 'BackgroundColor', [1 1 1]);
gui.drawDesignButton         = uicontrol('Style', 'pushbutton', 'String', 'Design', 'Callback', @draw_design);
gui.drawVcovButton           = uicontrol('Style', 'pushbutton', 'String', 'Covariance structure', 'Callback', @draw_vcov);
gui.saveModelValsMenu        = uicontrol('Style', 'popupmenu',  'String', {'All', 'Raw data', 'Parameter estimates', 'Residuals', 'Variance-covariance matrix'});
gui.saveModelValsTxtButton   = uicontrol('Style', 'pushbutton', 'String', 'Save to text');
gui.saveModelVals.WkspButton = uicontrol('Style', 'pushbutton', 'String', 'Return in workspace');

guidata(gcf,gui);
end


%==========================================================================
% Widgets placement
%==========================================================================
function position_widgets()

gui = getappdata(gcf, 'UsedByGUIData_m'); 
screen = get(0, 'ScreenSize');

gui.mainVBox   = uix.VBoxFlex('Parent', gcf, 'Spacing', 5, 'Padding', 5);
topHBox        = uix.HBox('Parent', gui.mainVBox, 'Spacing', 5);
bottomHBox     = uix.HBox('Parent', gui.mainVBox, 'Spacing', 5);

%-------------------------------------------------------------------------
% Image controls
%-------------------------------------------------------------------------
gui.imgControlVBox = uix.VBox('Parent', topHBox, 'Spacing', 5);

gui.controlPanel = uipanel('Parent', gui.imgControlVBox, 'Title', 'Image controls');
controlsHBox     = uix.HBox('Parent', gui.controlPanel, 'Spacing', 5, 'Padding', 5);

uicontrol('Parent', controlsHBox, 'Style', 'text', 'String', 'X');
gui.xCoordBox.Parent = controlsHBox;
uicontrol('Parent', controlsHBox, 'Style', 'text', 'String', 'Y');
gui.yCoordBox.Parent = controlsHBox;
uicontrol('Parent', controlsHBox, 'Style', 'text', 'String', 'Z');
gui.zCoordBox.Parent = controlsHBox;
gui.xhairCheck.Parent = controlsHBox;
uicontrol('Parent', controlsHBox, 'Style', 'text', 'String', 'Min');
gui.backgroundMinBox.Parent = controlsHBox;
uicontrol('Parent', controlsHBox, 'Style', 'text', 'String', 'Max');
gui.backgroundMaxBox.Parent = controlsHBox;
gui.autoBackContrastButton.Parent = controlsHBox;
uix.Empty('Parent', controlsHBox);

controlsHBox.Widths = [20 40 20 40 20 40 60 25 40 25 40 60 -1];

%-------------------------------------------------------------------------
% Images grid
%-------------------------------------------------------------------------
X = gui.images.backgroundHdr.dim(1);
Y = gui.images.backgroundHdr.dim(2);
Z = gui.images.backgroundHdr.dim(3);

% Aspect ratios
corAspRatio = X/Z;
axiAspRatio = Y/X;
sagAspRatio = Y/Z;

% All sizes derived from the coronal plot being 1/4 of the screen height
corPlotHeight = 1/4*screen(4);

% Image sizes (derived from the height of the coronal slice as 25% screen height)
gui.coronalPlot.Position(3)  = corAspRatio*corPlotHeight;
gui.coronalPlot.Position(4)  = corPlotHeight;
gui.axialPlot.Position(3)    = corAspRatio*corPlotHeight;
gui.axialPlot.Position(4)    = axiAspRatio*corAspRatio*corPlotHeight;
gui.saggitalPlot.Position(3) = sagAspRatio*corPlotHeight;
gui.saggitalPlot.Position(4) = corPlotHeight;

gui.imagesGrid = uix.Grid('Parent', gui.imgControlVBox, 'Spacing', 5);
%gui.imagesGrid.BackgroundColor = [0 0 0];

set(gui.coronalPlot,  'Parent', gui.imagesGrid);
set(gui.axialPlot,    'Parent', gui.imagesGrid);
set(gui.saggitalPlot, 'Parent', gui.imagesGrid);

% Info and Cbar panel
gui.imagesTabPanel           = uix.TabPanel('Parent', gui.imagesGrid, 'SelectionChangedFcn', @images_tab_changed);
infoTab                      = uix.Panel('Parent', gui.imagesTabPanel, 'Padding', 0);
gui.infoTable.Parent         = infoTab;
cbarTab                      = uix.Panel('Parent', gui.imagesTabPanel, 'Padding', 0);
gui.CBarPlot.Parent          = uicontainer('Parent',cbarTab);
gui.CBarPlot.Color           = cbarTab.BackgroundColor;
gui.CBarPlot.XColor          = cbarTab.BackgroundColor;
gui.CBarPlot.YColor          = cbarTab.BackgroundColor;
gui.imagesTabPanel.TabTitles = {'Info', 'CBar'};

gui.imagesGrid.Heights = [corPlotHeight axiAspRatio*corAspRatio*corPlotHeight];
gui.imagesGrid.Widths  = [corAspRatio*corPlotHeight sagAspRatio*corPlotHeight];

gui.imgControlVBox.Heights = [45 -1];

%--------------------------------------------------------------------------
% Results table
%--------------------------------------------------------------------------
gui.resultsTable.Parent = bottomHBox;
gui.resultsTable.ColumnWidth = {(gui.postEstToolsWin.Position(3)-10)/8};

%--------------------------------------------------------------------------
% Controls tab panel
%--------------------------------------------------------------------------
rightVBox            = uix.VBox('Parent', topHBox);
gui.controlsTabPanel = uix.TabPanel('Parent', rightVBox, 'SelectionChangedFcn', @tab_changed);
resultsTab           = uix.Panel('Parent', gui.controlsTabPanel, 'Padding', 10);
plotTab              = uix.Panel('Parent', gui.controlsTabPanel, 'Padding', 10);
assumpsTab           = uix.Panel('Parent', gui.controlsTabPanel, 'Padding', 5);
modelTab             = uix.Panel('Parent', gui.controlsTabPanel, 'Padding', 5);
if size(gui.MRM.Design.Y.Y) > 2
    ldaTab = uix.Panel('Parent', gui.controlsTabPanel, 'Padding', 5);
    gui.controlsTabPanel.TabTitles = {'Results', 'Plots', 'Assumps', 'Model', 'dLDA'};
else
    gui.controlsTabPanel.TabTitles = {'Results', 'Plots', 'Assumps', 'Model'};
end

% Results tab
resultsTabVBox                 = uix.VBox('Parent', resultsTab, 'Spacing', 5);
gui.resultsList.Parent         = resultsTabVBox;
threshPanel                    = uipanel('Parent', resultsTabVBox, 'Title', 'Threshold');
threshHBox                     = uix.HBox('Parent', threshPanel, 'Spacing', 5, 'Padding', 5);
gui.threshImgTick.Parent       = threshHBox;
gui.unthreshImgTick.Parent     = threshHBox;
gui.threshBox.Parent           = threshHBox;
uix.Empty('Parent', threshHBox)
threshHBox.Widths              = [130 130 45 -1];  
interpPanel                    = uipanel('Parent', resultsTabVBox, 'Title', 'Interpolation');
interpPanelVBox                = uix.VBox('Parent', interpPanel, 'Spacing', 5, 'Padding', 5);
interpHBox                     = uix.HBox('Parent', interpPanelVBox, 'Spacing', 5, 'Padding', 5);
gui.nnInterpTick.Parent        = interpHBox;
gui.trilInterpTick.Parent      = interpHBox;
uix.Empty('Parent', interpHBox);
interpHBox.Widths              = [120 100 -1];
colourmapPanel                 = uipanel('Parent', resultsTabVBox, 'Title', 'Colourmapping');
colourmapPanelVBox             = uix.VBox('Parent', colourmapPanel, 'Spacing', 5, 'Padding', 5);
gui.colourMapList.Parent       = colourmapPanelVBox;
%imageContrastPanel             = uipanel('Parent', colourmapPanelVBox, 'Title', 'Contrast');
imageContrastPanelHBox         = uix.HBox('Parent', colourmapPanelVBox, 'Spacing', 5, 'Padding', 5);
uicontrol('Parent', imageContrastPanelHBox, 'Style', 'text', 'String', 'Min');
gui.overlayMinBox.Parent       = imageContrastPanelHBox;
uicontrol('Parent', imageContrastPanelHBox, 'Style', 'text', 'String', 'Max');
gui.overlayMaxBox.Parent       = imageContrastPanelHBox;
gui.autoContrastButton.Parent  = imageContrastPanelHBox;
uix.Empty('Parent', imageContrastPanelHBox);
imageContrastPanelHBox.Widths  = [30 50 30 50 60 -1];
colourmapPanelVBox.Heights     = [25 35];
uix.Empty('Parent', resultsTabVBox);
resultsTabVBox.Heights         = [30 50 50 85 -1];

% Plots tab
plotTabHBox                    = uix.HBox('Parent',plotTab, 'Spacing', 10);
consPanel                      = uipanel('Parent', plotTabHBox, 'Title', 'Plot contrasts');
consPanelVBox                  = uix.VBox('Parent', consPanel, 'Spacing', 5, 'Padding', 10);
gui.plotConList.Parent         = consPanelVBox;
consListButtonHBox             = uix.HBox('Parent', consPanelVBox, 'Spacing', 5);
gui.newConButton.Parent        = consListButtonHBox;
gui.viewConButton.Parent       = consListButtonHBox;
gui.deleteConButton.Parent     = consListButtonHBox;
consPanelVBox.Heights          = [-1 35];
plotRightVBox                  = uix.VBox('Parent', plotTabHBox,'Spacing',5);
errorBarPanel                  = uipanel('Parent', plotRightVBox, 'Title', 'Error bars');
errorBarHBox                   = uix.HBox('Parent', errorBarPanel, 'Spacing', 5, 'Padding', 5);
gui.SECheck.Parent             = errorBarHBox;
gui.CICheck.Parent             = errorBarHBox;
gui.CILevel.Parent             = errorBarHBox;
uix.Empty('Parent', errorBarHBox);
errorBarHBox.Widths            = [60 50 40 -1];
plotTypePanel                  = uipanel('Parent', plotRightVBox, 'Title', 'Plot type');
plotTypeVBox                   = uix.VBox('Parent', plotTypePanel, 'Spacing', 5, 'Padding', 5);
plotTypeHBox                   = uix.HBox('Parent', plotTypeVBox, 'Spacing', 5);
gui.plot2dTick.Parent          = plotTypeHBox;
gui.plot3dTick.Parent          = plotTypeHBox;
uix.Empty('Parent', plotTypeHBox);
plotTypeHBox.Widths            = [60 60 -1];
neighbourhoodPanel             = uipanel('Parent', plotTypeVBox, 'Title', 'Neighbourhood');
neighbourhoodPanelHBox         = uix.HBox('Parent', neighbourhoodPanel, 'Spacing', 5, 'Padding', 5);
gui.neighbour6Radio.Parent     = neighbourhoodPanelHBox;
gui.neighbour18Radio.Parent    = neighbourhoodPanelHBox;
gui.neighbour26Radio.Parent    = neighbourhoodPanelHBox;
%gui.neighbourClustRadio.Parent = neighbourhoodPanelHBox;
uix.Empty('Parent', neighbourhoodPanelHBox);
%neighbourhoodPanelHBox.Widths  = [50 50 50 75 -1];
neighbourhoodPanelHBox.Widths  = [50 50 50 -1];
plotTypeVBox.Heights           = [-0.75 -1];
plotValsPanel                  = uipanel('Parent', plotRightVBox, 'Title', 'Plot values');
plotValsPanelHBox              = uix.HBox('Parent', plotValsPanel, 'Spacing', 5, 'Padding', 5);
uix.Empty('Parent', plotValsPanelHBox);
gui.savePlotValsButton.Parent  = plotValsPanelHBox;
gui.plotValsWorkspace.Parent   = plotValsPanelHBox;
uix.Empty('Parent', plotValsPanelHBox);
plotValsPanelHBox.Widths       = [-1 130 130 -1];
uix.Empty('Parent', plotRightVBox);
plotRightVBox.Heights          = [45 110 60 -1];
plotTabHBox.Widths             = [-1 -1.4];

% Assumptions tab

% Model tab
modelTabVBox                        = uix.VBox('Parent', modelTab, 'Spacing', 5, 'Padding', 10);
gui.modelValsTabPanel               = uix.TabPanel('Parent', modelTabVBox, 'SelectionChangedFcn', @vals_tab_changed);
PEsTab                              = uix.Panel('Parent', gui.modelValsTabPanel, 'Padding', 0);
vcovTab                             = uix.Panel('Parent', gui.modelValsTabPanel, 'Padding', 0);
PEsTabHBox                          = uix.HBox('Parent', PEsTab,  'Spacing', 0, 'Padding', 0);
gui.paramEstsTable.Parent           = PEsTabHBox;
vcovTabHBox                         = uix.HBox('Parent', vcovTab, 'Spacing', 0, 'Padding', 0);
gui.vcovTable.Parent                = vcovTabHBox;
gui.modelValsTabPanel.TabTitles     = {'PEs', 'VCov'};

visualisePanel                      = uipanel('Parent', modelTabVBox, 'Title', 'Visualisations');
visualisePanelHBox                  = uix.HBox('Parent', visualisePanel, 'Spacing', 5, 'Padding', 10);
gui.drawDesignButton.Parent         = visualisePanelHBox;
gui.drawVcovButton.Parent           = visualisePanelHBox;
uix.Empty('Parent', visualisePanelHBox);
visualisePanelHBox.Widths           = [130 130 -1];

saveValsPanel                       = uipanel('Parent', modelTabVBox, 'Title', 'Save model values');
saveValsPanelHBox                   = uix.HBox('Parent', saveValsPanel, 'Spacing', 5, 'Padding', 10);
gui.saveModelValsMenu.Parent        = saveValsPanelHBox;
gui.saveModelValsTxtButton.Parent   = saveValsPanelHBox;
gui.saveModelVals.WkspButton.Parent = saveValsPanelHBox;
uix.Empty('Parent', saveValsPanelHBox);
if size(gui.MRM.Design.Y.Cell,2) > 1
    for i = 1:size(gui.MRM.Design.Y.Cell,2)
        gui.paramEstsTable.ColumnName{i} = gui.MRM.Design.Y.Cell{i}.Label;
        gui.vcovTable.ColumnName{i}      = gui.MRM.Design.Y.Cell{i}.Label;
        gui.vcovTable.RowName{i}         = gui.MRM.Design.Y.Cell{i}.Label;
    end
else
    gui.paramEstsTable.ColumnName = [];
    gui.vcovTable.ColumnName      = [];
    gui.vcovTable.RowName{1}      = 'MSE';
end
for i = 1:size(gui.MRM.Design.X.Cell,2)
    gui.paramEstsTable.RowName{i} = gui.MRM.Design.X.Cell{i}.Label;
end
saveValsPanelHBox.Widths = [-1 130 130 -1];
modelTabVBox.Heights = [-1 65 65]; 

% dLDA tab
if size(gui.MRM.Design.Y.Y,2) > 1
    
end

%--------------------------------------------------------------------------
% Plots axes
%--------------------------------------------------------------------------
gui.barPlotAxes.Parent = uicontainer('Parent', rightVBox);
rightVBox.Heights = [-1.15 -1];

%--------------------------------------------------------------------------
% Main VBox
%--------------------------------------------------------------------------
gui.mainVBox.Heights = [-1*(gui.imagesGrid.Heights(1)+gui.imagesGrid.Heights(2)+gui.imgControlVBox.Heights(1)+10) ...
                        -1*(gui.postEstToolsWin.Position(4)-gui.imagesGrid.Heights(1)-gui.imagesGrid.Heights(2)-gui.imgControlVBox.Heights(1)-25)];

%gui.postEstToolsWin.Position(3) = 2*(corAspRatio*corPlotHeight + sagAspRatio*corPlotHeight + 15);

guidata(gcf,gui);
end


%==========================================================================
% Load template image used for the background
%==========================================================================
function load_background_img(path)
gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.images.backgroundVol    = spm_data_read(spm_vol(path));
gui.images.backgroundHdr    = spm_data_hdr_read(spm_vol(path));
gui.images.backgroundRange  = quantile(gui.images.backgroundVol(gui.images.backgroundVol(:) ~= 0), [.05 .95]);
gui.images.iMat             = (inv(gui.images.backgroundHdr.mat))';
gui.images.coordsMat        = [0 0 0 1] * gui.images.iMat;
gui.images.coordsMNI        = [0 0 0 1];
gui.backgroundMinBox.String = num2str(gui.images.backgroundRange(1));
gui.backgroundMaxBox.String = num2str(gui.images.backgroundRange(2));
guidata(gcf,gui);
end


%==========================================================================
% Draw the background images and crosshairs
%==========================================================================
function draw_images()

colormap([gray(64); get_overlay_colourmap()]);
redraw_axial();
redraw_saggital();
redraw_coronal();
fill_info_table();
end


%==========================================================================
% Coronal mouse click
%==========================================================================
function coronalPlot_mouseClick(~,~)

gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.dragging = 1;

matCoords = get(gui.coronalPlot, 'CurrentPoint');
matCoords = matCoords(1,1:2);

oldX = gui.images.coordsMat(1);
oldY = gui.images.coordsMat(2);
oldZ = gui.images.coordsMat(3);

if gui.xhairCheck.Value == 1
    gui.images.corHLine.YData = [matCoords(2) matCoords(2)];
    gui.images.corVLine.XData = [matCoords(1) matCoords(1)];
end

gui.images.coordsMat = [round(matCoords(1)) oldY round(gui.coronalPlot.YLim(2) - matCoords(2) + 1)];
gui.images.coordsMat(1) = clamp(gui.images.coordsMat(1), 1, gui.images.backgroundHdr.dim(1));
gui.images.coordsMat(2) = clamp(gui.images.coordsMat(2), 1, gui.images.backgroundHdr.dim(2));
gui.images.coordsMat(3) = clamp(gui.images.coordsMat(3), 1, gui.images.backgroundHdr.dim(3));

gui.images.coordsMNI = [gui.images.backgroundHdr.mat * [gui.images.coordsMat 1]']';

guidata(gcf,gui);

if gui.images.coordsMat(1) ~= oldX
    redraw_saggital();
else
    redraw_saggital_xhairs();
end

if gui.images.coordsMat(3) ~= oldZ
    redraw_axial();
else
    redraw_axial_xhairs();
end

update_info_table();
update_plot();
fill_param_ests_table();

end

%==========================================================================
% Saggital mouse click
%==========================================================================
function saggitalPlot_mouseClick(~,~)

gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.dragging = 2;

matCoords = get(gui.saggitalPlot, 'CurrentPoint');
matCoords = matCoords(1,1:2);

oldX = gui.images.coordsMat(1);
oldY = gui.images.coordsMat(2);
oldZ = gui.images.coordsMat(3);

if gui.xhairCheck.Value == 1
    gui.images.sagHLine.YData = [matCoords(2) matCoords(2)];
    gui.images.sagVLine.XData = [matCoords(1) matCoords(1)];
end

gui.images.coordsMat = [oldX round(matCoords(1)) round(gui.saggitalPlot.YLim(2) - matCoords(2) + 1)];
gui.images.coordsMat(1) = clamp(gui.images.coordsMat(1), 1, gui.images.backgroundHdr.dim(1));
gui.images.coordsMat(2) = clamp(gui.images.coordsMat(2), 1, gui.images.backgroundHdr.dim(2));
gui.images.coordsMat(3) = clamp(gui.images.coordsMat(3), 1, gui.images.backgroundHdr.dim(3));

gui.images.coordsMNI = [gui.images.backgroundHdr.mat * [gui.images.coordsMat 1]']';

guidata(gcf,gui);

if gui.images.coordsMat(2) ~= oldY
    redraw_coronal();
else
    redraw_coronal_xhairs();
end

if gui.images.coordsMat(3) ~= oldZ
    redraw_axial();
else
    redraw_axial_xhairs();
end

update_info_table();
update_plot();
fill_param_ests_table();

end

%==========================================================================
% Axial mouse click
%==========================================================================
function axialPlot_mouseClick(~,~)

gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.dragging = 3;

matCoords = get(gui.axialPlot, 'CurrentPoint');
matCoords = matCoords(1,1:2);

oldX = gui.images.coordsMat(1);
oldY = gui.images.coordsMat(2);
oldZ = gui.images.coordsMat(3);

if gui.xhairCheck.Value == 1
    gui.images.axiHLine.YData = [matCoords(2) matCoords(2)];
    gui.images.axiVLine.XData = [matCoords(1) matCoords(1)];
end

gui.images.coordsMat = [round(matCoords(1)) round(gui.axialPlot.YLim(2) - matCoords(2) + 1) oldZ];
gui.images.coordsMat(1) = clamp(gui.images.coordsMat(1), 1, gui.images.backgroundHdr.dim(1));
gui.images.coordsMat(2) = clamp(gui.images.coordsMat(2), 1, gui.images.backgroundHdr.dim(2));
gui.images.coordsMat(3) = clamp(gui.images.coordsMat(3), 1, gui.images.backgroundHdr.dim(3));

gui.images.coordsMNI = [gui.images.backgroundHdr.mat * [gui.images.coordsMat 1]']';

guidata(gcf,gui);

if gui.images.coordsMat(1) ~= oldX
    redraw_saggital();
else
    redraw_saggital_xhairs();
end

if gui.images.coordsMat(2) ~= oldY
    redraw_coronal();
else
    redraw_coronal_xhairs();
end

update_info_table();
update_plot();
fill_param_ests_table();

end


%==========================================================================
% Redraw coronal slice
%==========================================================================
function redraw_coronal()
gui = getappdata(gcf, 'UsedByGUIData_m'); 

% Background
gui.images.backgroundCorSlice = rot90(squeeze(gui.images.backgroundVol(:,gui.images.coordsMat(2),:)));

if ~isfield(gui.images, 'backgroundCorImg')
    %gui.images.backgroundCorImg   = image(gui.images.backgroundCorSlice, 'Parent', gui.coronalPlot);
    gui.images.backgroundCorImg   = pcolor(gui.images.backgroundCorSlice, 'Parent', gui.coronalPlot);
    shading(gui.coronalPlot, 'interp');
    gui.images.backgroundCorImg.CDataMapping = 'direct';
end
gui.images.backgroundCorImg.CData = min(64, round((64-1) * (gui.images.backgroundCorSlice - gui.images.backgroundRange(1))/ ...
                                                           (gui.images.backgroundRange(2) - gui.images.backgroundRange(1))) + 1);
                                                        
xhairsOnTop = 0;
                                                       
% Overlay
if gui.overlay == 1
    gui.images.overlayCorSlice = rot90(squeeze(gui.images.overlayVol(:,gui.images.coordsMat(2),:)));
    gui.images.overlayCorSlice(gui.images.overlayCorSlice < gui.images.overlayThresh) = 0;
    alphamap = ones(size(gui.images.overlayCorSlice));
    alphamap(gui.images.overlayCorSlice == 0) = 0;
    if ~isfield(gui.images, 'overlayCorImg')
        hold(gui.coronalPlot, 'on')
        gui.images.overlayCorImg = image(gui.images.overlayCorSlice, 'AlphaData', alphamap, 'Parent', gui.coronalPlot);
        hold(gui.coronalPlot, 'off')
        % Needed to force xhairs ontop of overlay w/o calling uistack
        xhairsOnTop = 1;
        delete(gui.coronalPlot.Children(2:3));
    else
        gui.images.overlayCorImg.AlphaData = alphamap;
    end
    gui.images.overlayCorImg.CData = 64 + min(64, round((64-1) * (gui.images.overlayCorSlice - gui.images.overlayRange(1))/ ...
                                                                 (gui.images.overlayRange(2) - gui.images.overlayRange(1))) + 1); 
    if max(gui.images.overlayCorImg.CData(:)) <= min(gui.images.backgroundCorImg.CData(:))
        caxis(gui.coronalPlot, [min(gui.images.backgroundCorImg.CData(:)) Inf]); 
    else
        caxis(gui.coronalPlot, [min(gui.images.backgroundCorImg.CData(:)) max(gui.images.overlayCorImg.CData(:))]);
    end
else
    caxis(gui.coronalPlot, [min(gui.images.backgroundCorImg.CData(:)) max(gui.images.backgroundCorImg.CData(:))]);
end

gui.coronalPlot.XTick = []; gui.coronalPlot.XTickLabel = [];
gui.coronalPlot.YTick = []; gui.coronalPlot.YTickLabel = [];
gui.coronalPlot.Box = 'on'; gui.coronalPlot.YDir = 'reverse';

if sign(gui.images.backgroundHdr.mat(1,1)) == -1
    gui.coronalPlot.XDir = 'reverse';
end

% Crosshairs
if gui.xhairCheck.Value == 1
    if ~isfield(gui.images, 'corHLine') || xhairsOnTop == 1
        gui.images.corHLine = line('XData',        [0 gui.coronalPlot.XLim(2)],                             ...
                                   'YData',        [gui.coronalPlot.YLim(2) - gui.images.coordsMat(3) + 1   ...
                                                    gui.coronalPlot.YLim(2) - gui.images.coordsMat(3) + 1], ...
                                   'LineStyle',     '-',                                                    ...
                                   'LineWidth',     1,                                                      ...
                                   'Color',         'blue',                                                 ...        
                                   'Parent',        gui.coronalPlot,                                        ...
                                   'ButtonDownFcn', @coronalPlot_mouseClick);
        gui.images.corVLine = line('XData',         [gui.images.coordsMat(1) gui.images.coordsMat(1)], ...
                                   'YData',         [0 gui.coronalPlot.YLim(2)],                       ...
                                   'LineStyle',     '-',                                               ...
                                   'LineWidth',     1,                                                 ...
                                   'Color',         'blue',                                            ...
                                   'Parent',        gui.coronalPlot,                                   ...
                                   'ButtonDownFcn', @coronalPlot_mouseClick);
    else
        redraw_coronal_xhairs();
    end
else
    if isfield(gui.images, 'corHLine')
        gui.images.corHLine.Visible = 'off';
        gui.images.corVLine.Visible = 'off';
    end
end

if gui.overlay == 1                       
    set(gui.images.overlayCorImg,    'ButtonDownFcn', @coronalPlot_mouseClick); 
else
    set(gui.images.backgroundCorImg, 'ButtonDownFcn', @coronalPlot_mouseClick);
end

guidata(gcf,gui);                       
end

%==========================================================================
% Redraw saggital slice
%==========================================================================
function redraw_saggital()
gui = getappdata(gcf, 'UsedByGUIData_m'); 

% Background image
gui.images.backgroundSagSlice     = rot90(squeeze(gui.images.backgroundVol(gui.images.coordsMat(1),:,:)));
if ~isfield(gui.images, 'backgroundSagImg')
    %gui.images.backgroundSagImg   = image(gui.images.backgroundSagSlice, 'Parent', gui.saggitalPlot);
    gui.images.backgroundSagImg   = pcolor(gui.images.backgroundSagSlice, 'Parent', gui.saggitalPlot);
    shading(gui.saggitalPlot, 'interp');
    gui.images.backgroundSagImg.CDataMapping = 'direct';
end
gui.images.backgroundSagImg.CData = min(64, round((64-1) * (gui.images.backgroundSagSlice - gui.images.backgroundRange(1))/ ...
                                                           (gui.images.backgroundRange(2) - gui.images.backgroundRange(1))) + 1);
xhairsOnTop = 0;

% Overlay
if gui.overlay  == 1
    gui.images.overlaySagSlice = rot90(squeeze(gui.images.overlayVol(gui.images.coordsMat(1),:,:)));
    gui.images.overlaySagSlice(gui.images.overlaySagSlice < gui.images.overlayThresh) = 0;
    alphamap = ones(size(gui.images.overlaySagSlice));
    alphamap(gui.images.overlaySagSlice == 0) = 0;
    if ~isfield(gui.images, 'overlaySagImg')
        hold(gui.saggitalPlot, 'on')
        gui.images.overlaySagImg = image(gui.images.overlaySagSlice, 'AlphaData', alphamap, 'Parent', gui.saggitalPlot);
        hold(gui.saggitalPlot, 'off')
        % Needed to force xhairs ontop of overlay w/o calling uistack
        xhairsOnTop = 1;
        delete(gui.saggitalPlot.Children(2:3));
    else
        gui.images.overlaySagImg.AlphaData = alphamap;
    end
    gui.images.overlaySagImg.CData = 64 + min(64, round((64-1) * (gui.images.overlaySagSlice - gui.images.overlayRange(1))/ ...
                                                                 (gui.images.overlayRange(2) - gui.images.overlayRange(1))) + 1); 
    if max(gui.images.overlaySagImg.CData(:)) <= min(gui.images.backgroundSagImg.CData(:))
        caxis(gui.saggitalPlot, [min(gui.images.backgroundSagImg.CData(:)) Inf]); 
    else
        caxis(gui.saggitalPlot, [min(gui.images.backgroundSagImg.CData(:)) max(gui.images.overlaySagImg.CData(:))]);
    end
else
    caxis(gui.saggitalPlot, [min(gui.images.backgroundSagImg.CData(:)) max(gui.images.backgroundSagImg.CData(:))]);
end

gui.saggitalPlot.XTick = []; gui.saggitalPlot.XTickLabel = [];
gui.saggitalPlot.YTick = []; gui.saggitalPlot.YTickLabel = [];
gui.saggitalPlot.Box = 'on'; gui.saggitalPlot.YDir = 'reverse';

% Crosshairs
if gui.xhairCheck.Value == 1
    if ~isfield(gui.images, 'sagHLine') || xhairsOnTop == 1
        gui.images.sagHLine = line('XData',         [0 gui.saggitalPlot.XLim(2)],                             ...
                                   'YData',         [gui.saggitalPlot.YLim(2) - gui.images.coordsMat(3) + 1   ...
                                                     gui.saggitalPlot.YLim(2) - gui.images.coordsMat(3) + 1], ...
                                   'LineStyle',     '-',                                                      ...
                                   'LineWidth',     1,                                                        ...
                                   'Color',         'blue',                                                   ...
                                   'Parent',        gui.saggitalPlot,                                         ...
                                   'ButtonDownFcn', @saggitalPlot_mouseClick);

        gui.images.sagVLine = line('XData',         [gui.images.coordsMat(2) gui.images.coordsMat(2)],        ...
                                   'YData',         [0 gui.saggitalPlot.YLim(2)],                             ...
                                   'LineStyle',     '-',                                                      ...
                                   'LineWidth',     1,                                                        ...
                                   'Color',         'blue',                                                   ...
                                   'Parent',        gui.saggitalPlot,                                         ...
                                   'ButtonDownFcn', @saggitalPlot_mouseClick);
    else
        redraw_saggital_xhairs();
    end
else
    if isfield(gui.images, 'sagHLine')
        gui.images.sagHLine.Visible = 'off';
        gui.images.sagVLine.Visible = 'off';
    end
end
 
if gui.overlay == 1
    set(gui.images.overlaySagImg, 'ButtonDownFcn', @saggitalPlot_mouseClick);
else
    set(gui.images.backgroundSagImg,    'ButtonDownFcn', @saggitalPlot_mouseClick);
end

guidata(gcf,gui);
end

%==========================================================================
% Redraw axial slice
%==========================================================================
function redraw_axial()
gui = getappdata(gcf, 'UsedByGUIData_m'); 

% Background
gui.images.backgroundAxiSlice     = rot90(squeeze(gui.images.backgroundVol(:,:,gui.images.coordsMat(3))));
if ~isfield(gui.images, 'backgroundAxiImg')
    %gui.images.backgroundAxiImg   = image(gui.images.backgroundAxiSlice, 'Parent', gui.axialPlot);
    gui.images.backgroundAxiImg   = pcolor(gui.images.backgroundAxiSlice, 'Parent', gui.axialPlot);
    shading(gui.axialPlot, 'interp');
    gui.images.backgroundAxiImg.CDataMapping = 'direct';
end
gui.images.backgroundAxiImg.CData = min(64, round((64-1) * (gui.images.backgroundAxiSlice - gui.images.backgroundRange(1))/ ...
                                                           (gui.images.backgroundRange(2) - gui.images.backgroundRange(1))) + 1);        
xhairsOnTop = 0;

% Overlay
if gui.overlay  == 1
    gui.images.overlayAxiSlice = rot90(squeeze(gui.images.overlayVol(:,:,gui.images.coordsMat(3))));
    gui.images.overlayAxiSlice(gui.images.overlayAxiSlice < gui.images.overlayThresh) = 0;
    alphamap = ones(size(gui.images.overlayAxiSlice));
    alphamap(gui.images.overlayAxiSlice == 0) = 0;
    if ~isfield(gui.images, 'overlayAxiImg')
        hold(gui.axialPlot, 'on')
        gui.images.overlayAxiImg = image(gui.images.overlayAxiSlice, 'AlphaData', alphamap, 'Parent', gui.axialPlot);
        hold(gui.axialPlot, 'off')
        % Needed to force xhairs ontop of overlay w/o calling uistack
        xhairsOnTop = 1;
        delete(gui.axialPlot.Children(2:3));
    else
        gui.images.overlayAxiImg.AlphaData = alphamap;
    end
    gui.images.overlayAxiImg.CData = 64 + min(64, round((64-1) * (gui.images.overlayAxiSlice - gui.images.overlayRange(1))/ ...
                                                                 (gui.images.overlayRange(2) - gui.images.overlayRange(1))) + 1); 
    if max(gui.images.overlayAxiImg.CData(:)) <= min(gui.images.backgroundAxiImg.CData(:))
        caxis(gui.axialPlot, [min(gui.images.backgroundAxiImg.CData(:)) Inf]); 
    else
        caxis(gui.axialPlot, [min(gui.images.backgroundAxiImg.CData(:)) max(gui.images.overlayAxiImg.CData(:))]);
    end
else
    caxis(gui.axialPlot, [min(gui.images.backgroundAxiImg.CData(:)) max(gui.images.backgroundAxiImg.CData(:))])
end

gui.axialPlot.XTick = []; gui.axialPlot.XTickLabel = [];
gui.axialPlot.YTick = []; gui.axialPlot.YTickLabel = [];
gui.axialPlot.Box = 'on'; gui.axialPlot.YDir = 'reverse';

if sign(gui.images.backgroundHdr.mat(1,1)) == -1
    gui.axialPlot.XDir = 'reverse';
end

% Crosshairs
if gui.xhairCheck.Value == 1
    if ~isfield(gui.images, 'axiHLine') || xhairsOnTop == 1
        gui.images.axiHLine = line('XData',         [0 gui.axialPlot.XLim(2)],                            ...
                                   'YData',         [gui.axialPlot.YLim(2) - gui.images.coordsMat(2) + 1  ...
                                                    gui.axialPlot.YLim(2) - gui.images.coordsMat(2) + 1], ...
                                   'LineStyle',     '-',                                                  ...
                                   'LineWidth',     1,                                                    ...
                                   'Color',         'blue',                                               ...
                                   'Parent',        gui.axialPlot,                                        ...
                                   'ButtonDownFcn', @axialPlot_mouseClick);
        gui.images.axiVLine = line('XData',         [gui.images.coordsMat(1) gui.images.coordsMat(1)],  ...
                                   'YData',         [0 gui.axialPlot.YLim(2)],                          ...
                                   'LineStyle',     '-',                                                ...
                                   'LineWidth',     1,                                                  ...
                                   'Color',         'blue',                                             ...
                                   'Parent',        gui.axialPlot,                                      ...
                                   'ButtonDownFcn', @axialPlot_mouseClick);
    else
        redraw_axial_xhairs();
    end
else
    if isfield(gui.images, 'axiHLine')
        gui.images.axiHLine.Visible = 'off';
        gui.images.axiVLine.Visible = 'off';
    end
end
                       
if gui.overlay == 1
    set(gui.images.overlayAxiImg,    'ButtonDownFcn', @axialPlot_mouseClick);
else
    set(gui.images.backgroundAxiImg, 'ButtonDownFcn', @axialPlot_mouseClick);
end

guidata(gcf,gui);
end


%==========================================================================
% Redraw axial xhairs
%==========================================================================
function redraw_axial_xhairs()
gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.images.axiHLine.YData = [gui.axialPlot.YLim(2) - gui.images.coordsMat(2) + 1   ...
                             gui.axialPlot.YLim(2) - gui.images.coordsMat(2) + 1];
gui.images.axiVLine.XData = [gui.images.coordsMat(1) gui.images.coordsMat(1)];
if gui.xhairCheck.Value == 1
    if strcmp(gui.images.axiHLine.Visible, 'off');
        gui.images.axiHLine.Visible = 'on';
        gui.images.axiVLine.Visible = 'on';
    end
end
guidata(gcf,gui);
end

%==========================================================================
% Redraw coronal xhairs
%==========================================================================
function redraw_coronal_xhairs()
gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.images.corHLine.YData = [gui.coronalPlot.YLim(2) - gui.images.coordsMat(3) + 1   ...
                             gui.coronalPlot.YLim(2) - gui.images.coordsMat(3) + 1];
gui.images.corVLine.XData = [gui.images.coordsMat(1) gui.images.coordsMat(1)];
if gui.xhairCheck.Value == 1
    if strcmp(gui.images.corHLine.Visible, 'off');
        gui.images.corHLine.Visible = 'on';
        gui.images.corVLine.Visible = 'on';
    end
end
guidata(gcf,gui);
end

%==========================================================================
% Redraw saggital xhairs
%==========================================================================
function redraw_saggital_xhairs()
gui = getappdata(gcf, 'UsedByGUIData_m'); 
gui.images.sagHLine.YData = [gui.saggitalPlot.YLim(2) - gui.images.coordsMat(3) + 1   ...
                             gui.saggitalPlot.YLim(2) - gui.images.coordsMat(3) + 1];
gui.images.sagVLine.XData = [gui.images.coordsMat(2) gui.images.coordsMat(2)];
if gui.xhairCheck.Value == 1
    if strcmp(gui.images.sagHLine.Visible, 'off');
        gui.images.sagHLine.Visible = 'on';
        gui.images.sagVLine.Visible = 'on';
    end
end
guidata(gcf,gui);
end



%==========================================================================
% Fill the info table
%==========================================================================
function fill_info_table()

gui = getappdata(gcf, 'UsedByGUIData_m'); 

data(1,1) = {'Millimeters'}; data(1,2) = {num2str(gui.images.coordsMNI(1:3))};
data(2,1) = {'Voxels'};      data(2,2) = {num2str(gui.images.coordsMat(1:3))};
data(3,1) = {'Value'};       data(3,2) = {num2str(gui.images.backgroundVol(gui.images.coordsMat(1), ...
                                                                           gui.images.coordsMat(2), ...
                                                                           gui.images.coordsMat(3)))};
data(4,1) = {'Location'};    data(4,2) = {'Anterior commissure'};
data(5,1) = {'Affine mat'};  data(5,2) = {num2str(gui.images.backgroundHdr.mat(1,:))};
data(6,1) = {''};            data(6,2) = {num2str(gui.images.backgroundHdr.mat(2,:))};  
data(7,1) = {''};            data(7,2) = {num2str(gui.images.backgroundHdr.mat(3,:))};  
gui.infoTable.Data = data;

gui.infoTable.ColumnWidth = {1/3*gui.saggitalPlot.Position(3) - 1, 2/3*gui.saggitalPlot.Position(3) - 1};

% Coordinate boxes
gui.xCoordBox.String = num2str(gui.images.coordsMNI(1));
gui.yCoordBox.String = num2str(gui.images.coordsMNI(2));
gui.zCoordBox.String = num2str(gui.images.coordsMNI(3));

guidata(gcf, gui);
end

%==========================================================================
% Update the info table
%==========================================================================
function update_info_table()

gui = getappdata(gcf, 'UsedByGUIData_m'); 

gui.infoTable.Data(1,2) = {num2str(gui.images.coordsMNI(1:3))};
gui.infoTable.Data(2,2) = {num2str(gui.images.coordsMat(1:3))};

if gui.overlay ~= 1
    gui.infoTable.Data(3,2) = {num2str(gui.images.backgroundVol(gui.images.coordsMat(1), ...
                                                                gui.images.coordsMat(2), ...
                                                                gui.images.coordsMat(3)))};
else
    gui.infoTable.Data(3,2) = {num2str(gui.images.overlayVol(gui.images.coordsMat(1), ...
                                                             gui.images.coordsMat(2), ...
                                                             gui.images.coordsMat(3)))};
end

gui.infoTable.Data(4,2) = {'Anterior commissure'};

% Coordinate boxes
gui.xCoordBox.String = num2str(gui.images.coordsMNI(1));
gui.yCoordBox.String = num2str(gui.images.coordsMNI(2));
gui.zCoordBox.String = num2str(gui.images.coordsMNI(3));

guidata(gcf, gui);
end


%==========================================================================
% Load overlay
%==========================================================================
function load_overlay(path)

gui = getappdata(gcf, 'UsedByGUIData_m'); 

try
    tempVol = spm_data_read(spm_vol(path));
catch err
    msgbox(['Could not load image. Error: ' err.message])
    return
end

tempHdr = spm_data_hdr_read(spm_vol(path));

if gui.images.backgroundHdr.dim(1) ~= tempHdr.dim(1) || ...
   gui.images.backgroundHdr.dim(2) ~= tempHdr.dim(2) || ...
   gui.images.backgroundHdr.dim(3) ~= tempHdr.dim(3)
    
    gui.images.overlayVol = NaN(gui.images.backgroundHdr.dim(1:3));
    for i = 1:gui.images.backgroundHdr.dim(3)
        M = inv(spm_matrix([0 0 -i]) / (gui.images.backgroundHdr.mat) * tempHdr.mat);
        gui.images.overlayVol(:,:,i) = spm_slice_vol(tempVol, M, gui.images.backgroundHdr.dim(1:2), gui.overlayInterp);
    end
else
    gui.images.overlayVol = tempVol;
end

gui.images.overlayHdr  = tempHdr;
gui.images.overlayiMat = (inv(tempHdr.mat))'; 
clearvars 'tempImg' 'tempHdr';

gui.images.overlayVol(isnan(gui.images.overlayVol(:))) = 0;
gui.images.overlayRange  = [min(gui.images.overlayVol(gui.images.overlayVol(:) ~= 0)) ...
                            quantile(gui.images.overlayVol(gui.images.overlayVol(:) ~= 0), 0.95)];
if isnan(gui.images.overlayRange(1)) || isnan(gui.images.overlayRange(2))
    gui.images.overlayRange = [0 0];
end
gui.overlayMinBox.String = num2str(gui.images.overlayRange(1));
gui.overlayMaxBox.String = num2str(gui.images.overlayRange(2));
gui.images.overlayThresh = str2double(gui.threshBox.String);
gui.overlay  = 1;

guidata(gcf, gui);

redraw_saggital();
redraw_axial();
redraw_coronal();

create_colourbar();

end


%==========================================================================
% Refresh the plot list
%==========================================================================
function fill_plot_list()
gui = getappdata(gcf, 'UsedByGUIData_m'); 
strings = gui.plotConList.String;
strings{1} = 'None';
if isfield(gui.MRM.Contrasts, 'Plots')
    for i = 1:size(gui.MRM.Contrasts.Plots.Con,2)
        strings{i+1} = gui.MRM.Contrasts.Plots.Con{i}.Name;
    end
end
gui.plotConList.String = strings;
guidata(gcf,gui);
end


%==========================================================================
% Item selected in the plot list has changed
%==========================================================================
function selected_plot_changed(~,~)
setup_plot();
end


%==========================================================================
% Create bar plot
%==========================================================================
function setup_plot()

gui = getappdata(gcf, 'UsedByGUIData_m'); 

con_num = gui.plotConList.Value-1;

if con_num < 1 || ~isfield(gui.images, 'overlayVol')
    clear_plot();
    return
end

A      = gui.MRM.Contrasts.Plots.Con{con_num}.A;
C      = gui.MRM.Contrasts.Plots.Con{con_num}.C;
coords = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
coords(coords <= 0) = 1;

if gui.plot3dTick.Value == 1
    neighbourCoords = neighbourhood_coordinates(coords(1),coords(2),coords(3));
    means = zeros(size(neighbourCoords,1), size(A,1) * size(C,1));
    for i = 1:size(neighbourCoords,1)
        beta       = gui.paramEsts(:,:,neighbourCoords(i,1),neighbourCoords(i,2),neighbourCoords(i,3));
        means(i,:) = bar_heights(beta,A,C);
    end
    bar3(gui.barPlotAxes, means);
else
    beta  = gui.paramEsts(:,:,coords(1),coords(2),coords(3));
    means = bar_heights(beta,A,C);
    bar(means, 'facecolor', [0, 0.75, 0.75], 'Parent', gui.barPlotAxes);
    xlim(gui.barPlotAxes, [0.5 size(means,2)+0.5]);
    add_error_bars();
end
gui.barPlotAxes.ButtonDownFcn = @plot_clicked;
guidata(gcf,gui);
end

%==========================================================================
% Update bar plot data
%==========================================================================
function update_plot()

gui = getappdata(gcf, 'UsedByGUIData_m'); 
con_num = gui.plotConList.Value-1;
if con_num < 1
    clear_plot();
    return
end
A      = gui.MRM.Contrasts.Plots.Con{con_num}.A;
C      = gui.MRM.Contrasts.Plots.Con{con_num}.C;
coords = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
coords(coords <= 0) = 1;

if gui.plot3dTick.Value == 1
    neighbourCoords = neighbourhood_coordinates(coords(1),coords(2),coords(3));
    means = zeros(size(neighbourCoords,1), size(A,1) * size(C,1));
    for i = 1:size(neighbourCoords,1)
        beta       = gui.paramEsts(:,:,neighbourCoords(i,1),neighbourCoords(i,2),neighbourCoords(i,3));
        means(i,:) = bar_heights(beta,A,C);
    end
    [~,x,~,xx,yy,~,~,~,zz] = makebars(means,'3');
    for i=1:size(yy,2)/4
         gui.barPlotAxes.Children(i).XData = xx+x(i);
         gui.barPlotAxes.Children(i).YData = yy(:,(i-1)*4+(1:4));
         gui.barPlotAxes.Children(i).ZData = zz(:,(i-1)*4+(1:4));
    end
else
    if length(gui.barPlotAxes.Children) > 1
       delete(gui.barPlotAxes.Children(1:length(gui.barPlotAxes.Children)-1)); 
    end
    beta  = gui.paramEsts(:,:,coords(1),coords(2),coords(3));
    means = bar_heights(beta,A,C);
    gui.barPlotAxes.Children(1).YData = means;
end
guidata(gcf,gui);
end

%==========================================================================
% Bar heights from contrast
%==========================================================================
function heights = bar_heights(beta, A, C)
heights  = zeros(1, size(A,1) * size(C,1));
k = 1;
for i = 1:size(A,1)
    for j = 1:size(C,1)      
        a = A(i, :);
        c = C(j, :);
        heights(k) = a * beta * c';
        k = k+1;
    end
end
end

%==========================================================================
% Add error bars to plot
%==========================================================================
function add_error_bars()

gui = getappdata(gcf, 'UsedByGUIData_m'); 

if gui.plotConList.Value == 1 || gui.plot3dTick.Value == 1
    return
end

coords = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
coords(coords <= 0) = 1;

V = gui.vcov(:,:,coords(1),coords(2),coords(3));
A = gui.MRM.Contrasts.Plots.Con{gui.plotConList.Value-1}.A;
C = gui.MRM.Contrasts.Plots.Con{gui.plotConList.Value-1}.C;
X = gui.MRM.Design.X.X;

alpha = 1 - (str2double(gui.CILevel.String)/100);
if gui.CICheck.Value == 1
    zVal  = spm_invNcdf(1 - (alpha/2));
elseif gui.SECheck.Value == 1
    zVal = 1;
end

CI = error_limits(V,A,C,X,zVal);

for i = 1:length(gui.barPlotAxes.Children.XData)
   line([gui.barPlotAxes.Children(i).XData(i)       gui.barPlotAxes.Children(i).XData(i)], ...
        [gui.barPlotAxes.Children(i).YData(i)-CI(i) gui.barPlotAxes.Children(i).YData(i)+CI(i)], ...
        'Parent',gui.barPlotAxes, 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1); 
end

guidata(gcf,gui);

end

%==========================================================================
% Error bars from contrast
%==========================================================================
function lims = error_limits(V,A,C,X,zVal)
lims = zeros(1, size(A,1) * size(C,1));
invXX = inv(X'*X);
k = 1;
for i = 1:size(A,1)
    for j = 1:size(C,1)       
        a = A(i, :);
        c = C(j, :);
        lims(k)  = (c * V * c') * (a * invXX * a');
        k = k + 1;
    end
end
lims = sqrt(lims);
lims = lims * zVal;
end


%==========================================================================
% Get the coordinates for neighbourbood of voxels centered on X,Y,Z
%==========================================================================
function coords = neighbourhood_coordinates(X,Y,Z)
gui = getappdata(gcf, 'UsedByGUIData_m'); 
if gui.neighbour6Radio.Value == 1
    offsets = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1; 0 0 0];
elseif gui.neighbour18Radio.Value == 1
    offsets = [ 0 0  1;  0  1  0;  1  0 0; -1  0  0;  0 -1  0;  0  0 -1; ...
                0 1  1;  1  0  1;  1  1 0; -1 -1  0; -1  0 -1;  0 -1 -1; ...
               -1 0  1; -1  1  0;  0 -1 1;  0  1 -1;  1 -1  0;  1  0 -1; ...
                0 0  0];
elseif gui.neighbour26Radio.Value == 1
    offsets = [ 0 0  1;  0  1  0;  1  0 0; -1  0  0;  0 -1  0;  0  0 -1; ...
                0 1  1;  1  0  1;  1  1 0; -1 -1  0; -1  0 -1;  0 -1 -1; ...
               -1 0  1; -1  1  0;  0 -1 1;  0  1 -1;  1 -1  0;  1  0 -1; ...
                1 1  1; -1 -1 -1; -1  1 1;  1 -1  1;  1  1 -1; -1 -1  1; ...
               -1 1 -1;  1 -1 -1;  0  0 0]; 
end
coords = offsets + repmat([X, Y, Z], size(offsets,1),1);
end


%==========================================================================
% Plot 2D and 3D tick boxes callbacks
%==========================================================================
function plotTypeCallback(obj,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if obj == gui.plot3dTick
   if gui.plot3dTick.Value == 1
       gui.plot2dTick.Value          = 0;
       gui.neighbour6Radio.Enable    = 'on';
       gui.neighbour18Radio.Enable   = 'on';
       gui.neighbour26Radio.Enable   = 'on';
       gui.savePlotValsButton.Enable = 'off';
       gui.plotValsWorkspace.Enable  = 'off';
   else
       gui.plot2dTick.Value = 1;
       gui.neighbour6Radio.Enable    = 'off';
       gui.neighbour18Radio.Enable   = 'off';
       gui.neighbour26Radio.Enable   = 'off';
       gui.savePlotValsButton.Enable = 'on';
       gui.plotValsWorkspace.Enable  = 'on';
       h = rotate3d(gui.barPlotAxes);
       h.Enable = 'off';
   end
else
   if gui.plot2dTick.Value == 1
       gui.plot3dTick.Value  = 0;
       gui.neighbour6Radio.Enable    = 'off';
       gui.neighbour18Radio.Enable   = 'off';
       gui.neighbour26Radio.Enable   = 'off';
       gui.savePlotValsButton.Enable = 'on';
       gui.plotValsWorkspace.Enable  = 'on';
       h = rotate3d(gui.barPlotAxes);
       h.Enable = 'off';
   else
       gui.plot3dTick.Value = 1;
       gui.neighbour6Radio.Enable    = 'on';
       gui.neighbour18Radio.Enable   = 'on';
       gui.neighbour26Radio.Enable   = 'on';
       gui.savePlotValsButton.Enable = 'off';
       gui.plotValsWorkspace.Enable  = 'off';
   end
end
guidata(gcf,gui);
setup_plot();
end


%==========================================================================
% Load results table
%==========================================================================
function fill_results_table()

gui = getappdata(gcf, 'UsedByGUIData_m');

name = regexprep(gui.MRM.Contrasts.Con{gui.resultsList.Value-1}.Name, ' ', '_');

if strcmp(gui.MRM.Options.Thresh.Pvals, 'Permutation')
    filepath = [gui.MRM.Options.Out filesep 'PermutationContrasts' filesep name filesep ...
                'MRM_' name '_Results.txt'];
elseif strcmp(gui.MRM.Options.Thresh.Pvals, 'Approximate')
    filepath = [gui.MRM.Options.Out filesep 'Contrasts' filesep name filesep            ...
                'MRM_' name '_Results.txt'];
end

tablePath = fileparts(filepath);

if ~exist(tablePath,'dir') || strcmp(gui.MRM.Options.Thresh.Level, 'None')
   return
end

fileID = fopen(filepath);
input  = textscan(fileID,'%s %d %.3f %.3f %.3f %d %d %d','Delimiter','|', ...
                  'HeaderLines',1);

inputNew = cell(size(input{1},1), 3);
    
temp = input{3}(:);
tot  = sum(mod(temp,1) == 0);

names = get(gui.resultsTable, 'ColumnName');

if tot ~= size(temp,1)
    names{3} = 'Mass';
else
    names{3} = 'Extent';
end

for i = 1:size(input{1},1)
    inputNew{i,1} = input{1}{i};   
    inputNew{i,2} = sprintf('%d',input{2}(i));
    if strcmp(names{3}, 'Mass') == 1
        inputNew{i,3} = sprintf('%.3f',input{3}(i));
    else
        inputNew{i,3} = sprintf('%d',input{3}(i));
    end
    inputNew{i,4} = sprintf('%.3f',input{4}(i));
    if strcmp(sprintf('%.3f',input{5}(i)), '0.000') == 1;
        inputNew{i,5} = '< 0.001';
    else
        inputNew{i,5} = ['   ' sprintf('%.3f',input{5}(i))];
    end
    for j = 6:8
        inputNew{i,j} = input{j}(i);
    end
end

switch gui.MRM.Options.Thresh.Level
        case 'None'
              names{5} = 'Uncorrected p-value (voxel)';
        case 'Uncorrected'
              names{5} = 'Uncorrected p-value (voxel)';
        case 'Voxel' 
            switch gui.MRM.Options.Thresh.Method
                case 'FWE'
                    names{5} = 'FWE p-value (voxel)';
                case 'FDR'
                    names{5} = 'FDR q-value (voxel)';     
            end
        case 'Cluster'
            names{4} = 'FWE p-value (cluster)';
            names{5} = 'X';
            names{6} = 'Y';
            names{7} = 'Z';
            names{8} = '';
end

fclose(fileID);

columnformat = {'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char'};
           
gui.resultsTable.ColumnFormat = columnformat;
gui.resultsTable.Data         = inputNew;
gui.resultsTable.ColumnName   = names;

guidata(gcf,gui);

end


%==========================================================================
% Fill estimates tables - PEs and VCov
%==========================================================================
function fill_param_ests_table()

gui = getappdata(gcf, 'UsedByGUIData_m');

if gui.controlsTabPanel.Selection ~= 4
    return
end

% PEs
if gui.modelValsTabPanel.Selection == 1
    if isfield(gui, 'paramEsts')
        if isempty(gui.paramEsts)
            return
        end
    else
        return
    end
    coords                  = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
    coords(coords <= 0)     = 1;
    beta                    = gui.paramEsts(:,:,coords(1),coords(2),coords(3));
    gui.paramEstsTable.Data = num2cell(beta);
    guidata(gcf,gui);
% Vcov    
elseif gui.modelValsTabPanel.Selection == 2
    if isfield(gui, 'vcov')
        if isempty(gui.vcov)
            return
        end
    else
        return
    end
    coords              = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
    coords(coords <= 0) = 1;
    vcov                = gui.vcov(:,:,coords(1),coords(2),coords(3));
    gui.vcovTable.Data  = num2cell(vcov);
    guidata(gcf,gui);
end

end


%==========================================================================
% Load parameter estimates into memory
%==========================================================================
function load_parameter_ests()

gui = getappdata(gcf, 'UsedByGUIData_m');

if isfield(gui, 'paramEsts')
    if ~isempty(gui.paramEsts)
        return
    end
end

gui.paramEsts = zeros(size(gui.MRM.Design.X.X,2), size(gui.MRM.Design.Y.Y,2),     ...
                      gui.images.overlayHdr.dim(1), gui.images.overlayHdr.dim(2), ...
                      gui.images.overlayHdr.dim(3));

for i = 1:size(gui.MRM.Design.X.X,2)
    for j = 1:size(gui.MRM.Design.Y.Y,2)
        gui.paramEsts(i,j,:,:,:) = spm_data_read(spm_vol([gui.MRM.Options.Out filesep 'MRM_PE_' num2str(i) '_' num2str(j) '.nii']));
    end
end

%%%% May need to remove below if mem is limited

if isfield(gui, 'vcov')
    if ~isempty(gui.vcov)
        return
    end
end

gui.vcov = zeros(size(gui.MRM.Design.Y.Y,2), size(gui.MRM.Design.Y.Y,2),          ...
                      gui.images.overlayHdr.dim(1), gui.images.overlayHdr.dim(2), ...
                      gui.images.overlayHdr.dim(3));

for i = 1:size(gui.MRM.Design.Y.Y,2)
    for j = 1:size(gui.MRM.Design.Y.Y,2)
        gui.vcov(i,j,:,:,:) = spm_data_read(spm_vol([gui.MRM.Options.Out filesep 'MRM_Covar_' num2str(i) '_' num2str(j) '.nii']));                       
    end
end

%%%%

guidata(gcf,gui);
end


%==========================================================================
% Neighbourhood radio button callbacks
%==========================================================================
function neighbourCallback(obj,~)
gui = getappdata(gcf, 'UsedByGUIData_m');

if obj == gui.neighbour6Radio
    gui.neighbour6Radio.Value     = 1;
    gui.neighbour18Radio.Value    = 0;
    gui.neighbour26Radio.Value    = 0; 
    %gui.neighbourClustRadio.Value = 0; 
elseif obj == gui.neighbour18Radio
    gui.neighbour6Radio.Value     = 0;
    gui.neighbour18Radio.Value    = 1;
    gui.neighbour26Radio.Value    = 0; 
    %gui.neighbourClustRadio.Value = 0; 
elseif obj == gui.neighbour26Radio
    gui.neighbour6Radio.Value     = 0;
    gui.neighbour18Radio.Value    = 0;
    gui.neighbour26Radio.Value    = 1;
    %gui.neighbourClustRadio.Value = 0; 
elseif obj == gui.neighbourClustRadio.Value
    gui.neighbour6Radio.Value     = 0;
    gui.neighbour18Radio.Value    = 0;
    gui.neighbour26Radio.Value    = 0;
    %gui.neighbourClustRadio.Value = 1; 
end
setup_plot();
guidata(gcf,gui);
end


%==========================================================================
% Change the overlay threshold
%==========================================================================
function thresh_image_callback(obj,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if obj == gui.threshImgTick 
    gui.threshImgTick.Value = 1;
    gui.unthreshImgTick.Value = 0;
    gui.threshBox.Enable = 'off';
    gui.threshBox.String = '0';
    guidata(gcf,gui);
    load_overlay(gui.MRM.Contrasts.Con{gui.resultsList.Value-1}.FileThresh);
elseif obj == gui.unthreshImgTick
    gui.threshImgTick.Value = 0;
    gui.unthreshImgTick.Value = 1;
    gui.threshBox.Enable = 'on';
    guidata(gcf,gui);
    load_overlay(gui.MRM.Contrasts.Con{gui.resultsList.Value-1}.File);
end
end


%==========================================================================
% Change interpolation method
%==========================================================================
function interp_change_callback(obj,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if obj == gui.nnInterpTick
    gui.nnInterpTick.Value   = 1;
    gui.trilInterpTick.Value = 0;
elseif obj == gui.trilInterpTick
    gui.nnInterpTick.Value   = 0;
    gui.trilInterpTick.Value = 1;
end
gui.overlayInterp = gui.trilInterpTick.Value;
guidata(gcf,gui);
selectNum = gui.resultsList.Value-1;
if selectNum > 0
    if gui.threshImgTick.Value == 1
        load_overlay(gui.MRM.Contrasts.Con{selectNum}.FileThresh);
    else
        load_overlay(gui.MRM.Contrasts.Con{selectNum}.File);
    end
end
end


%==========================================================================
% Tab changed
%==========================================================================
function tab_changed(~,event)
if event.NewValue == 4
   fill_param_ests_table(); 
end
end

function vals_tab_changed(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if isfield(gui, 'controlsTabPanel') % stop it firing on creation
    if gui.controlsTabPanel.Selection == 4
        fill_param_ests_table();
    end
end
end

function images_tab_changed(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if isfield(gui, 'imagesTabPanel') % stop it firing on creation
    create_colourbar();
end
end


%==========================================================================
% Fill the list of results
%==========================================================================
function fill_results_list()
gui = getappdata(gcf, 'UsedByGUIData_m');

strings = cell(gui.MRM.Contrasts.Number + 1,1);
strings{1} = 'None';

if gui.MRM.Contrasts.Number > 0
    for i = 1:gui.MRM.Contrasts.Number
        strings{i+1} = gui.MRM.Contrasts.Con{i}.Name; 
    end
end

gui.resultsList.String = strings;
guidata(gcf,gui);
end


%==========================================================================
% Selected result changed
%==========================================================================
function selected_result_changed(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
selectNum = gui.resultsList.Value-1;
if selectNum > 0
    if gui.threshImgTick.Value == 1
        load_overlay(gui.MRM.Contrasts.Con{selectNum}.FileThresh);
    else
        load_overlay(gui.MRM.Contrasts.Con{selectNum}.File);
    end
    fill_results_table();
    load_parameter_ests();
end
end


%==========================================================================
% New plot contrast
%==========================================================================
function new_plot_con(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
MRM = MRM_contrastSpecifier(gui.MRM,'plot');
if ~isempty(MRM)
    gui.MRM = MRM;
end
guidata(gcf, gui);
fill_plot_list();
end


%==========================================================================
% Return plot values to workspace
%==========================================================================
function return_plot_vals(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');

if gui.resultsList.Value == 1 || gui.plotConList.Value == 1
    return
end

varName = inputdlg({'Bar heights', 'Errors'}, 'Names', 1, {'Bars', 'Errors'});

if ~isempty(varName)
    coords   = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
    coords(coords <= 0) = 1;
    A        = gui.MRM.Contrasts.Plots.Con{gui.plotConList.Value-1}.A;
    C        = gui.MRM.Contrasts.Plots.Con{gui.plotConList.Value-1}.C;
    alpha    = 1 - (str2double(gui.CILevel.String)/100);
    if gui.CICheck.Value == 1
        zVal = spm_invNcdf(1 - (alpha/2));
    elseif gui.SECheck.Value == 1
        zVal = 1;
    end
    beta     = gui.paramEsts(:,:,coords(1),coords(2),coords(3));
    heights  = bar_heights(beta,A,C);
    V        = gui.vcov(:,:,coords(1),coords(2),coords(3));
    X        = gui.MRM.Design.X.X;
    CI       = error_limits(V,A,C,X,zVal);
    try
        assignin('base', varName{1}, heights);
    catch err
        msgbox(['Error returning values: ' err.message]);
    end
    
    try
        assignin('base', varName{2}, CI);
    catch err
         msgbox(['Error returning values: ' err.message]);
    end
end

end


%==========================================================================
% Save plot values to text file
%==========================================================================
function save_plot_text(~,~)

gui = getappdata(gcf, 'UsedByGUIData_m');

if gui.resultsList.Value == 1 || gui.plotConList.Value == 1
    return
end

[FileName,PathName] = uiputfile('*.txt', 'Save plot values as');

if ~isempty(FileName) && FileName ~= 0
    coords   = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
    coords(coords <= 0) = 1;
    A        = gui.MRM.Contrasts.Plots.Con{gui.plotConList.Value-1}.A;
    C        = gui.MRM.Contrasts.Plots.Con{gui.plotConList.Value-1}.C;
    alpha    = 1 - (str2double(gui.CILevel.String)/100);
    if gui.CICheck.Value == 1
        zVal = spm_invNcdf(1 - (alpha/2));
    elseif gui.SECheck.Value == 1
        zVal = 1;
    end
    beta     = gui.paramEsts(:,:,coords(1),coords(2),coords(3));
    heights  = bar_heights(beta,A,C);
    V        = gui.vcov(:,:,coords(1),coords(2),coords(3));
    X        = gui.MRM.Design.X.X;
    CI       = error_limits(V,A,C,X,zVal);   
    fileID   = fopen([PathName filesep FileName], 'w');
    outTable = [heights', CI'];   
    displaytable(outTable, {'Bar heights', 'Errors'}, [11 9], {'.4f','.4f'}, [], fileID);
    fclose(fileID);
end

end


%==========================================================================
% Clear the plot axes
%==========================================================================
function clear_plot()
gui = getappdata(gcf, 'UsedByGUIData_m');
cla(gui.barPlotAxes,'reset');
gui.barPlotAxes.Box = 'on';
guidata(gcf,gui);
end


%==========================================================================
% Mouse movement
%==========================================================================
function mouse_motion(~,~)
if ~strcmp(get(gcf, 'Name'), 'MRM Post-estimation tools')
    return
end
gui = getappdata(gcf, 'UsedByGUIData_m');
if ~isempty(gui.dragging)
    switch gui.dragging
        case 1
            coronalPlot_mouseClick();
        case 2
            saggitalPlot_mouseClick();
        case 3
            axialPlot_mouseClick();
    end
elseif gui.plot3dTick.Value == 1 % Control the 3D rotation
    h = rotate3d(gui.barPlotAxes);
    if sum(in_slice_axes()) == 0
        if strcmp(h.Enable, 'off')
            h.Enable = 'on';
        end
    else
        h.Enable = 'off';
    end
end
end


%==========================================================================
% Mouse released
%==========================================================================
function mouse_release(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
gui.dragging = [];
guidata(gcf,gui);
if sum(in_slice_axes()) > 0
    add_error_bars();
    remove_table_highlight();
else
    table_resize();
end
end


%==========================================================================
% Window resizing
%==========================================================================
function win_resize(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
gui.resultsTable.ColumnWidth = {(gui.postEstToolsWin.Position(3)-25)/8};
if isfield(gui, 'oldWinSize')
    if gui.oldWinSize(2) - gui.postEstToolsWin.Position(4) == 0 || ...
       gui.oldWinSize(1) - gui.postEstToolsWin.Position(3) == 0
        return
    end
    vChange = gui.postEstToolsWin.Position(4)/gui.oldWinSize(2);
    gui.imagesGrid.Heights = gui.imagesGrid.Heights*vChange;
    gui.imagesGrid.Widths  = gui.imagesGrid.Widths*vChange;
    gui.oldWinSize = [gui.postEstToolsWin.Position(3) gui.postEstToolsWin.Position(4)];
end
gui.oldTabSize = [gui.resultsTable.Position(3) gui.resultsTable.Position(4)];
guidata(gcf,gui);                
end


%==========================================================================
% Scale the slices as the table changes height
%==========================================================================
function table_resize()
gui = getappdata(gcf, 'UsedByGUIData_m');
if gui.oldTabSize(2) ~= gui.resultsTable.Position(4)
    vChange = (gui.oldTabSize(2)-gui.resultsTable.Position(4))/2;
    gui.imagesGrid.Heights = gui.imagesGrid.Heights+vChange;
    gui.imagesGrid.Widths  = gui.imagesGrid.Widths+vChange;
    gui.oldTabSize = [gui.resultsTable.Position(3) gui.resultsTable.Position(4)];
    guidata(gcf,gui);
end
end


%==========================================================================
% CLAMP!
%==========================================================================
function out = clamp(x,low,high)
if x < low
    out = low;
    return;
elseif x > high
    out = high;
    return;
else
    out = x;
end       
end


%==========================================================================
% Clourmap options
%==========================================================================
function map = get_overlay_colourmap()
gui = getappdata(gcf, 'UsedByGUIData_m');

switch gui.colourMapList.Value
    case 1
        map = hot(64);
    case 2 
        map = cool(64);
    case 3
        map = jet(64);
    case 4
        map = spring(64); 
    case 5
        map = summer(64); 
    case 6
        map = autumn(64); 
    case 7
        map = winter(64); 
    case 8
        map = parula(64); 
    case 9
        map = hsv(64);
end

end


%==========================================================================
% Change overlay clourmap
%==========================================================================
function change_overlay_colourmap(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
colormap(gui.axialPlot,    [gray(64); get_overlay_colourmap()]);
colormap(gui.saggitalPlot, [gray(64); get_overlay_colourmap()]);
colormap(gui.coronalPlot,  [gray(64); get_overlay_colourmap()]);
colormap(gui.CBarPlot,     [gray(64); get_overlay_colourmap()]);
redraw_axial();
redraw_saggital();
redraw_coronal();
create_colourbar();
end


%==========================================================================
% Change overlay min and max contrast values
%==========================================================================
function change_overlay_contrast(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');

min = str2double(gui.overlayMinBox.String);
max = str2double(gui.overlayMaxBox.String);

if isnan(min)
   min = gui.images.overlayRange(1);
   gui.overlayMinBox.String = num2str(min);
end

if isnan(max)
   max = gui.images.overlayRange(2);
   gui.overlayMaxBox.String = num2str(max);
end

gui.images.overlayRange = [min max];

guidata(gcf,gui);

redraw_axial();
redraw_saggital();
redraw_coronal();
create_colourbar();
end


%==========================================================================
% Change background min and max contrast values
%==========================================================================
function change_background_contrast (~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');

min = str2double(gui.backgroundMinBox.String);
max = str2double(gui.backgroundMaxBox.String);

if isnan(min)
   min = gui.images.backgroundRange(1);
   gui.backgroundMinBox.String = num2str(min);
end

if isnan(max)
   max = gui.images.backgroundRange(2);
   gui.backgroundMaxBox.String = num2str(max);
end

gui.images.backgroundRange = [min max];

guidata(gcf,gui);

redraw_axial();
redraw_saggital();
redraw_coronal();
create_colourbar();
end


%==========================================================================
% Auto contrast
%==========================================================================
function auto_overlay_contrast(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
gui.images.overlayRange  = [min(gui.images.overlayVol(gui.images.overlayVol(:) ~= 0)) ...
                            quantile(gui.images.overlayVol(gui.images.overlayVol(:) ~= 0), 0.95)];
gui.overlayMinBox.String = num2str(gui.images.overlayRange(1));
gui.overlayMaxBox.String = num2str(gui.images.overlayRange(2));
guidata(gcf,gui);
redraw_axial();
redraw_saggital();
redraw_coronal();
create_colourbar();
end

function auto_background_contrast(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
gui.images.backgroundRange  = quantile(gui.images.backgroundVol(gui.images.backgroundVol(:) ~= 0), [0.05 0.95]);
gui.backgroundMinBox.String = num2str(gui.images.backgroundRange(1));
gui.backgroundMaxBox.String = num2str(gui.images.backgroundRange(2));
guidata(gcf,gui);
redraw_axial();
redraw_saggital();
redraw_coronal();
create_colourbar();
end


%==========================================================================
% Change overlay threshold
%==========================================================================
function change_overlay_thresh(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if isnan(str2double(gui.threshBox.String))
    gui.threshBox.String = num2str(gui.images.overlayThresh);
end
gui.images.overlayThresh = str2double(gui.threshBox.String);
guidata(gcf,gui);
redraw_axial();
redraw_saggital();
redraw_coronal();
end


%==========================================================================
% Switch off the crosshairs
%==========================================================================
function xhairs_off(~,~)
redraw_saggital();
redraw_axial();
redraw_coronal();
end


%==========================================================================
% Are we inside any of the slices?
%==========================================================================
function result = in_slice_axes()
gui = getappdata(gcf, 'UsedByGUIData_m');

result = [0 0 0];

C = get(gui.saggitalPlot,'currentpoint');
xlim = get(gui.saggitalPlot,'xlim');
ylim = get(gui.saggitalPlot,'ylim');
outX = ~any(diff([xlim(1) C(1,1) xlim(2)])<0);
outY = ~any(diff([ylim(1) C(1,2) ylim(2)])<0);
if outX && outY
    result(1) = 1;
    return
end

C = get(gui.coronalPlot,'currentpoint');
xlim = get(gui.coronalPlot,'xlim');
ylim = get(gui.coronalPlot,'ylim');
outX = ~any(diff([xlim(1) C(1,1) xlim(2)])<0);
outY = ~any(diff([ylim(1) C(1,2) ylim(2)])<0);
if outX && outY
    result(2) = 1;
    return
end

C = get(gui.axialPlot,'currentpoint');
xlim = get(gui.axialPlot,'xlim');
ylim = get(gui.axialPlot,'ylim');
outX = ~any(diff([xlim(1) C(1,1) xlim(2)])<0);
outY = ~any(diff([ylim(1) C(1,2) ylim(2)])<0);
if outX && outY
    result(3) = 1;
    return
end

end


%==========================================================================
% Change confidence level
%==========================================================================
function CI_level_callback(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
val = str2double(gui.CILevel.String);
if isnan(val) || val > 99 || val < 0
    gui.CILevel.String = '95';
    guidata(gcf,gui);
end
update_plot();
add_error_bars();
end


%==========================================================================
% SE and CI tick boxes
%==========================================================================
function CILevel_callback(obj,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if obj == gui.SECheck
   if gui.SECheck.Value == 1
       gui.CICheck.Value  = 0;
       gui.SECheck.Value  = 1;
       gui.CILevel.Enable = 'off';
   else
       gui.SECheck.Value = 0;
       gui.CICheck.Value = 1;
       gui.CILevel.Enable = 'on';
   end
elseif obj == gui.CICheck
   if gui.CICheck.Value == 1
       gui.CICheck.Value = 1;
       gui.SECheck.Value = 0;
       gui.CILevel.Enable = 'on';
   else
       gui.SECheck.Value = 1;
       gui.CICheck.Value = 0;
       gui.CILevel.Enable = 'off';
   end
end
guidata(gcf,gui);
update_plot();
add_error_bars();
end


%==========================================================================
% Double click the plot to open in external window
%==========================================================================
function plot_clicked(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if strcmp(gui.postEstToolsWin.SelectionType, 'open')
    fig = figure;
    copyobj(gui.barPlotAxes,fig); 
end
end


%==========================================================================
% Window close
%==========================================================================
function win_close(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
MRM = gui.MRM; %#ok<NASGU>
save([gui.MRM.Options.Out filesep 'MRM.mat'], 'MRM');
delete(gcf);
end


%==========================================================================
% Draw design matrix
%==========================================================================
function draw_design(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
MRM_drawDesign(gui.MRM);
end


%==========================================================================
% Draw covariance matrix
%==========================================================================
function draw_vcov(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');
if gui.overlay  == 1
    coords = round([gui.images.coordsMNI(1) gui.images.coordsMNI(2) gui.images.coordsMNI(3) 1] * gui.images.overlayiMat);
    coords(coords <= 0) = 1;
    V = gui.vcov(:,:,coords(1),coords(2),coords(3));
    MRM_drawVcov(gui.MRM, V, gui.images.coordsMNI);
else
    msgbox('Please select a result first');
end
end


%==========================================================================
% Draw the colourbar
%==========================================================================
function create_colourbar()
gui = getappdata(gcf, 'UsedByGUIData_m');

if gui.imagesTabPanel.Selection ~= 2 || gui.overlay  ~= 1;
    return
end

if ~isfield(gui, 'cbar')
    gui.cbar = colorbar(gui.CBarPlot, 'WestOutside');
    gui.cbar.Limits(1) = 0.5;
    gui.cbar.TickDirection = 'in';
    gui.cbar.YAxisLocation = 'right';
    gui.cbar.Position(1) = 0.5*gui.cbar.Position(1);
    gui.cbar.Position(3) = 1.5*gui.cbar.Position(3);
end

min   = gui.images.overlayRange(1);
max   = gui.images.overlayRange(2);
diff  = max - min;
step  = diff/4;
ticks = round(min:step:max,2);

gui.cbar.Ticks = [0.5 0.625 0.75 0.875 1];
gui.cbar.TickLabels = cellstr(num2str(ticks(:)));

guidata(gcf,gui);
end


%==========================================================================
% Manual coordinates entered
%==========================================================================
function coords_changed(~,~)
gui = getappdata(gcf, 'UsedByGUIData_m');

X = str2double(gui.xCoordBox.String);
Y = str2double(gui.yCoordBox.String);
Z = str2double(gui.zCoordBox.String);

if isnan(X)
   X = gui.images.coordsMNI(1);
   gui.xCoordBox.String = num2str(X);
end
if isnan(Y)
   Y = gui.images.coordsMNI(2);
   gui.yCoordBox.String = num2str(Y);
end
if isnan(Z)
   Z = gui.images.coordsMNI(3);
   gui.zCoordBox.String = num2str(Z);
end

coords = round([X Y Z 1] * gui.images.iMat);

for i = 1:3
    if coords(i) <= 0
        coords(i) = 1;
    elseif coords(i) > gui.images.backgroundHdr.dim(i)
        coords(i) = gui.images.backgroundHdr.dim(i);
    end
end

gui.images.coordsMat = [coords(1) coords(2) coords(3)];
gui.images.coordsMNI = [gui.images.backgroundHdr.mat * [gui.images.coordsMat 1]']';

gui.xCoordBox.String = num2str(gui.images.coordsMNI(1));
gui.yCoordBox.String = num2str(gui.images.coordsMNI(2));
gui.zCoordBox.String = num2str(gui.images.coordsMNI(3));

guidata(gcf,gui);

redraw_axial();
redraw_saggital();
redraw_coronal();
update_info_table();
update_plot();
add_error_bars();
fill_param_ests_table();

end


%==========================================================================
% Table row clicked
%==========================================================================
function table_click(obj,data)
gui = getappdata(gcf, 'UsedByGUIData_m');
if numel(class(data)) == 23 % Bit hacky, but it works
    if data.getKeyCode ~= 38 && data.getKeyCode ~= 40 % Up or down arrow
        return
    end
end
row = obj.SelectedRow + 1;
gui.xCoordBox.String = gui.resultsTable.Data{row,6};
gui.yCoordBox.String = gui.resultsTable.Data{row,7};
gui.zCoordBox.String = gui.resultsTable.Data{row,8};
guidata(gcf,gui);
coords_changed([],[]);
end


%==========================================================================
% Turns on highlighting whole rows rather than just cells
%==========================================================================
function table_row_highlight()
gui = getappdata(gcf, 'UsedByGUIData_m');
hJScroll = findjobj(gui.resultsTable);
hJTable = hJScroll.getViewport.getView;
hJTable.setNonContiguousCellSelection(false);
hJTable.setColumnSelectionAllowed(false);
hJTable.setRowSelectionAllowed(true);
hJTable = handle(hJTable, 'CallbackProperties');
set(hJTable, 'MousePressedCallback', @table_click);
set(hJTable, 'KeyPressedCallback',   @table_click);
guidata(gcf,gui);
end


%==========================================================================
% De-select the table row
%==========================================================================
function remove_table_highlight()
gui = getappdata(gcf, 'UsedByGUIData_m');
temp = get(gui.resultsTable,'Data');
gui.resultsTable.Data = {'dummy'};
gui.resultsTable.Data = temp;
guidata(gcf,gui);
end

