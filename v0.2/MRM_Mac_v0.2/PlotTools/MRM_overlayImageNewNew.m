function MRM_overlayImageNewNew(overlay, MNIcoords, thresh, position, maxima, crosshairs, infoPane, template)
%==========================================================================
% Interactive overlay viewer window
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This function creates the GUI for the interactive display of overlays to
% be used in conjunction with the MRM post-estimation tools
%
%=========================================================================%
% Copyright 2015 Martyn McFarquhar                                        %
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

global MRM
global MRMoptions
global sliceHandles
global OverlayMatCoords
global maximaMat

maximaMat = maxima;

if license('test', 'image_toolbox') ~= 1
     msgbox('The image processing toolbox is required to display overlays, sorry.');
     return
end

% Clear the window before plotting
%---------------------------------
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1 
    WholeFigHand = findobj('type','figure','name','MRMviewer');
    set(0, 'currentfigure', WholeFigHand);
    cla(sliceHandles(1));
    cla(sliceHandles(2));
    cla(sliceHandles(3));
    cla(sliceHandles(4));
end
      

subplot = @(m,n,p) subtightplot(m,n,p);

% Figure details
%---------------
if isempty(findobj('type','figure','name','MRMviewer')) == 1
    
    figure('Position',  position, ...
                           'Units',               'Pixels', ...
                           'NumberTitle',         'off',    ...
                           'Color',               'black',  ...
                           'MenuBar',             'None', ...
                           'DockControls',        'off', ...
                           'NumberTitle',         'off', ...
                           'ToolBar',             'none', ...
                           'WindowButtonDownFcn', @MouseClickCallback, ...
                           'CloseRequestFcn',     @CloseWindowCallback, ...
                           'Name',                'MRMviewer');

    % Drop-down menu
    toolMenu = uimenu('Label','Tools');
    
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
    
        uimenu(toolMenu,'Label','Go to nearest maximum','Callback',@goToMaximum);

        if strcmp(MRM.Options.Thresh.Level, 'Cluster') == 1
            uimenu(toolMenu,'Label','Go to largest cluster', 'Callback',@goToLargestClust);
        else
            uimenu(toolMenu,'Label','Go to global maximum', 'Callback',@goToGlobalMaximum);
        end
    end
    
    uimenu(toolMenu,'Label','Save image as jpeg','Callback',@saveImage);
end

sliceHandles = zeros(4,1);

templatePath = MRMoptions.Template;
tarhdr       = spm_data_hdr_read(templatePath);

% Find the matrix coordinates from MNI
TemplateMatCoords  = MRM_mni2mat(MNIcoords, templatePath);

% Slices
SagSlice = squeeze(template(TemplateMatCoords(1), :, :));
CorSlice = squeeze(template(:, TemplateMatCoords(2), :));
AxiSlice = squeeze(template(:,:,TemplateMatCoords(3)));

% Find the matrix coordinates from MNI
OverlayMatCoords    = [MNIcoords 1] * (inv(tarhdr.mat))';
OverlayMatCoords(4) = [];
OverlayMatCoords    = round(OverlayMatCoords);

for i = 1:3
    if sign(OverlayMatCoords(i)) == -1
        OverlayMatCoords(i) = 1;
    end
end

%==========================================================================
% Saggital
%==========================================================================

% Extract slice implied by coordinates
SagDisplay = squeeze(overlay(OverlayMatCoords(1),:,:));

% Threshold the overlay
SagDisplay(SagDisplay < thresh) = NaN;

% Create subplot
subplot(2,2,2);

% Reorient and display the template slice
imagesc(fliplr(rot90(SagSlice,3)));
colormap gray;
hold on; 

% Display the Overlay slice and set the alpha transparency to NaN
hImg = subimage(fliplr(rot90(SagDisplay,3)), hot);
set(hImg, 'AlphaData', 1 - isnan(fliplr(rot90(SagDisplay,3))));

% Reverse Y-axis so 0 is the bottom of the brain
set(gca,'YDir','normal');

% Save the handle for the saggital plot
sliceHandles(1) = gca;


%--------------------------------------------------------------------------
% Saggital crosshairs
%--------------------------------------------------------------------------
xLims = get(gca, 'xlim');
yLims = get(gca, 'ylim');

if crosshairs == 1
    
    % Sagital view - x-axis of subplot, Y-axis of the brain
    %              - y-axis of subplot, Z-axis of the brain

    % X
    line('XData', [0 xLims(2)],                                       ...
         'YData', [OverlayMatCoords(3) OverlayMatCoords(3)],          ...
         'LineStyle', '-',                                            ...
         'LineWidth', 1,                                              ...
         'Color', 'blue');

    % Y             
    line('XData', [OverlayMatCoords(2) OverlayMatCoords(2)],          ...
        'YData',  [0 yLims(2)],                                       ...
        'LineStyle', '-',                                             ...
        'LineWidth', 1,                                               ...
        'Color', 'blue');
end

axis off
hold off


%==========================================================================
% Coronal
%==========================================================================

% Extract slice implied by coordinates
CorDisplay = squeeze(overlay(:,OverlayMatCoords(2),:));

% Threshold the overlay
CorDisplay(CorDisplay < thresh) = NaN;

% Create subplot
subplot(2,2,1);

% Reorient and display the template slice
imagesc(fliplr(rot90(CorSlice,3)));
colormap gray

% Check for LR flips in the template image
if sign(tarhdr.mat(1,1)) == -1
    set(gca,'XDir','reverse');
end

hold on

% Display the Overlay slice and set the alpha transparency to NaN
hImg = subimage(fliplr(rot90(CorDisplay,3)),hot);
set(hImg, 'AlphaData', 1 - isnan(fliplr(rot90(CorDisplay,3))));

% Reverse Y-axis so 0 is the bottom of the brain
set(gca,'YDir','normal');

% Save the handle for the coronal plot
sliceHandles(2) = gca;

%--------------------------------------------------------------------------
% Coronal crosshairs
%--------------------------------------------------------------------------
xLims = get(gca, 'xlim');
yLims = get(gca, 'ylim');

if crosshairs == 1
    
    % Coronal view - x-axis of subplot, X-axis of the brain
    %              - y-axis of subplot, Z-axis of the brain
    
    % X 
    line('XData', [0 xLims(2)],                                       ...
         'YData', [OverlayMatCoords(3) OverlayMatCoords(3)],          ...
         'LineStyle', '-',                                            ...
         'LineWidth', 1,                                              ...
         'Color', 'blue');

    % Y           
    line('XData', [OverlayMatCoords(1) OverlayMatCoords(1)],           ...
        'YData',  [0 yLims(2)],                                        ...
        'LineStyle', '-',                                              ...
        'LineWidth', 1,                                                ...
        'Color', 'blue');
end

axis off 
hold off



%==========================================================================
% Axial image
%==========================================================================

AxiDisplay = overlay(:,:,OverlayMatCoords(3));
AxiDisplay(AxiDisplay < thresh) = NaN;

subplot(2,2,3);

imagesc(fliplr(rot90(AxiSlice,3)));
colormap gray

% Check for LR flips in the template image
if sign(tarhdr.mat(1,1)) == -1
    set(gca,'XDir','reverse');
end

hold on

hImg = subimage(fliplr(rot90(AxiDisplay,3)),hot);
set(hImg, 'AlphaData', 1 - isnan(fliplr(rot90(AxiDisplay,3))));

set(gca,'YDir','normal');

% Save the handle for the axial plot
sliceHandles(3) = gca;

%--------------------------------------------------------------------------
% Axial crosshairs
%--------------------------------------------------------------------------

xLims = get(gca, 'xlim');
yLims = get(gca, 'ylim');

if crosshairs == 1
    
    % Axial view - x-axis of subplot, X-axis of the brain
    %            - y-axis of subplot, Y-axis of the brain
    
    % X
    line('XData', [0 xLims(2)],                                       ...
         'YData', [OverlayMatCoords(2) OverlayMatCoords(2)],          ...
         'LineStyle', '-',                                            ...
         'LineWidth', 1,                                              ...
         'Color', 'blue');

    % Y             
    line('XData', [OverlayMatCoords(1) OverlayMatCoords(1)],           ...
        'YData',  [0 yLims(2)],                                        ... 
        'LineStyle', '-',                                              ...
        'LineWidth', 1,                                                ...
        'Color', 'blue');
end

axis off
hold off

%==========================================================================
% Lower pane info
%==========================================================================
subplot(2,2,4);
set(gca, 'color', 'black')
set(gca, 'units', 'normalized')

if infoPane == 1

    text(0.05, 0.6, ...
        ['Coordinates: ' num2str(MNIcoords(1)) ' ' num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3))], ...
        'HorizontalAlignment', 'left',       ...
        'VerticalAlignment',   'top',        ...
        'FontSize',             18,          ...
        'Clipping',            'off',        ...
        'Color',               'white');

    text(0.05, 0.5, ...
        ['F-value: ' num2str(overlay(OverlayMatCoords(1),    ...
                                     OverlayMatCoords(2),    ...
                                     OverlayMatCoords(3)))], ...
        'HorizontalAlignment', 'left',   ...
        'VerticalAlignment',   'top',    ...
        'FontSize',             18,      ...
        'Clipping',            'off',    ...
        'Color',               'white');

    text(0.05, 0.4, ...
        ['Label: ' MRM_getAtlasLabel(MRMoptions.Atlas, MNIcoords)], ...
        'HorizontalAlignment', 'left',   ...
        'VerticalAlignment',   'top',    ...
        'FontSize',             18,      ...
        'Clipping',            'off',    ...
        'Color',               'white');
end

sliceHandles(4) = gca;


%==========================================================================
% Mouse click
%==========================================================================
function MouseClickCallback(objectHandle, eventData)
    
   handles = guidata(findobj('type','figure','name','MRM Post-estimation Tools'));
   chs     = get(handles.crosshairsCheck, 'Value');
   infoPn  = get(handles.infoPaneBox, 'Value');
    
   temp    = handles.template;
   temphdr = spm_data_hdr_read(MRMoptions.Template);
    
   %=======================================================================
   % Saggital
   %=======================================================================
   if gca == sliceHandles(1) 
       
      coordinates = get(gca,'CurrentPoint');
      coordinates = coordinates(1,1:2);
      
      xLimits = get(gca, 'xlim');
      yLimits = get(gca, 'ylim');
      
      if coordinates(1) < xLimits(2) && coordinates(1) > xLimits(1) && ...
         coordinates(2) < yLimits(2) && coordinates(2) > yLimits(1)
          
          coords           = [OverlayMatCoords(1) coordinates(1) coordinates(2)];
          OverlayMatCoords = coords;
          
          MNICoords    = temphdr.mat * [coords 1]';
          MNICoords(4) = [];
          MNICoords    = round(MNICoords');
          
          if get(handles.overlayPopup, 'Value') == 2
              threshold = 0;
          else
              threshold  =  str2num(get(handles.thresholdText, 'String'));
          end
          
          if strcmp(MRM.Options.Thresh.Level, 'None') == 1
              threshold  =  str2num(get(handles.thresholdText, 'String'));
          end
          
          set(handles.XcoordText, 'String', num2str(MNICoords(1)));
          set(handles.YcoordText, 'String', num2str(MNICoords(2)));
          set(handles.ZcoordText, 'String', num2str(MNICoords(3)));
          
          % Remove all crosshairs and redraw for saggital
          if chs == 1
              
              for j = 1:3
                children = get(sliceHandles(j), 'children');
                delete(children(1));
                delete(children(2));
              end
         
              line('XData', [0 xLimits(2)],        ...
                   'YData', [coords(3) coords(3)], ...
                   'LineStyle', '-',               ...
                   'LineWidth', 1,                 ...
                   'Color', 'blue');

              line('XData', [coords(2) coords(2)],  ...
                   'YData',  [0 yLimits(2)],        ...
                   'LineStyle', '-',                ...
                   'LineWidth', 1,                  ...
                   'Color', 'blue');
          end
          
          % Redraw coronal
          slice = squeeze(temp(:, round(coords(2)), :));
          C     = squeeze(handles.reslicedOverlay(:,round(coords(2)),:));
          C(C < threshold) = NaN;
          redrawCoronal(slice, temphdr, C, coords, chs);
          
          % Redraw axial
          slice = squeeze(temp(:, :, round(coords(3))));
          A     = squeeze(handles.reslicedOverlay(:,:,round(coords(3))));
          A(A < threshold) = NaN;
          redrawAxial(slice, temphdr, A, coords, chs);
          
          % Redraw info pane
          if infoPn == 1
            redrawInfo([num2str(MNICoords(1)) ' ' num2str(MNICoords(2)) ' ' ...
                           num2str(MNICoords(3))], ...
                           num2str(handles.reslicedOverlay(round(coords(1)), ...
                                                           round(coords(2)), ...
                                                           round(coords(3)))));
          end
          
          % Redraw plot window (if open)
          redrawPlotWindow();
          
          % Tidy up
          clearvars A C
          
      end

   %=======================================================================
   % Coronal
   %=======================================================================   
   elseif gca == sliceHandles(2)
       
      coordinates = get(gca,'CurrentPoint');
      coordinates = coordinates(1,1:2);
      
      xLimits = get(gca, 'xlim');
      yLimits = get(gca, 'ylim');
      
      if coordinates(1) < xLimits(2) && coordinates(1) > xLimits(1) && ...
         coordinates(2) < yLimits(2) && coordinates(2) > yLimits(1)
          
          coords           = [coordinates(1) OverlayMatCoords(2) coordinates(2)];
          OverlayMatCoords = coords;
          
          MNICoords    = temphdr.mat * [coords 1]';
          MNICoords(4) = [];
          MNICoords    = round(MNICoords');
          
          if get(handles.overlayPopup, 'Value') == 2
              threshold = 0;
          else
              threshold = str2num(get(handles.thresholdText, 'String'));
          end
          
          if strcmp(MRM.Options.Thresh.Level, 'None') == 1
              threshold  =  str2num(get(handles.thresholdText, 'String'));
          end
          
          set(handles.XcoordText, 'String', num2str(MNICoords(1)));
          set(handles.YcoordText, 'String', num2str(MNICoords(2)));
          set(handles.ZcoordText, 'String', num2str(MNICoords(3)));

          % Remove all crosshairs and redraw for coronal
          if chs == 1
              
              for j = 1:3
                children = get(sliceHandles(j), 'children');
                delete(children(1));
                delete(children(2));
              end
         
              line('XData', [0 xLimits(2)],          ...
                  'YData', [coords(3) coords(3)],    ...
                  'LineStyle', '-',                  ...
                  'LineWidth', 1,                    ...
                  'Color', 'blue');
              
              line('XData', [coords(1) coords(1)],   ...
                  'YData',  [0 yLimits(2)],          ...
                  'LineStyle', '-',                  ...
                  'LineWidth', 1,                    ...
                  'Color', 'blue');
          end
           
          % Redraw axial
          slice = squeeze(temp(:, :, round(coords(3))));
          A     = squeeze(handles.reslicedOverlay(:,:,round(coords(3))));
          A(A < threshold) = NaN;
          redrawAxial(slice, temphdr, A, coords, chs);
          
          % Redraw saggital
          slice = squeeze(temp(round(coords(1)),:,:));
          S     = squeeze(handles.reslicedOverlay(round(coords(1)),:,:));
          S(S < threshold) = NaN;
          redrawSaggital(slice, [], S, coords, chs);
          
          % Redraw info pane
          if infoPn == 1
            redrawInfo([num2str(MNICoords(1)) ' ' num2str(MNICoords(2)) ' ' ...
                           num2str(MNICoords(3))], ...
                           num2str(handles.reslicedOverlay(round(coords(1)), ...
                                                           round(coords(2)), ...
                                                           round(coords(3)))));
          end
          
          % Redraw plot window if open
          redrawPlotWindow();
          
          % Tidy up
          clearvars A S
      end
   
   %-----------------------------------------------------------------------
   % Axial
   %-----------------------------------------------------------------------
   elseif gca == sliceHandles(3)
       
      coordinates = get(gca,'CurrentPoint');
      coordinates = coordinates(1,1:2);
      
      xLimits = get(gca, 'xlim');
      yLimits = get(gca, 'ylim');
      
      if coordinates(1) < xLimits(2) && coordinates(1) > xLimits(1) && ...
         coordinates(2) < yLimits(2) && coordinates(2) > yLimits(1)
          
          coords           = [coordinates(1) coordinates(2) OverlayMatCoords(3)];
          OverlayMatCoords = coords;
          
          MNICoords    = temphdr.mat * [coords 1]';
          MNICoords(4) = [];
          MNICoords    = round(MNICoords');
          
          if get(handles.overlayPopup, 'Value') == 2
              threshold = 0;
          else
              threshold = str2num(get(handles.thresholdText, 'String'));
          end
          
          if strcmp(MRM.Options.Thresh.Level, 'None') == 1
              threshold  =  str2num(get(handles.thresholdText, 'String'));
          end
          
          set(handles.XcoordText, 'String', num2str(MNICoords(1)));
          set(handles.YcoordText, 'String', num2str(MNICoords(2)));
          set(handles.ZcoordText, 'String', num2str(MNICoords(3)));
          
          % Remove all crosshairs and redraw for axial
          if chs == 1
              
              for j = 1:3
                children = get(sliceHandles(j), 'children');
                delete(children(1));
                delete(children(2));
              end
         
              line('XData', [0 xLimits(2)],        ...
                  'YData', [coords(2) coords(2)],  ...
                  'LineStyle', '-',                ...
                  'LineWidth', 1,                  ...
                  'Color', 'blue');
              
              line('XData', [coords(1) coords(1)], ...
                  'YData',  [0 yLimits(2)],        ...
                  'LineStyle', '-',                ...
                  'LineWidth', 1,                  ...
                  'Color', 'blue');
          end

          % Redraw saggital
          slice = squeeze(temp(round(coords(1)),:,:));
          S     = squeeze(handles.reslicedOverlay(round(coords(1)),:,:));
          S(S < threshold) = NaN;
          redrawSaggital(slice, [], S, coords, chs);
          
          % Redraw coronal
          slice = squeeze(temp(:, round(coords(2)), :));
          C     = squeeze(handles.reslicedOverlay(:,round(coords(2)),:));
          C(C < threshold) = NaN;
          redrawCoronal(slice, temphdr, C, coords, chs);
          
          % Redraw info pane
          if infoPn == 1
            redrawInfo([num2str(MNICoords(1)) ' ' num2str(MNICoords(2)) ' ' ...
                           num2str(MNICoords(3))], ...
                           num2str(handles.reslicedOverlay(round(coords(1)), ...
                                                           round(coords(2)), ...
                                                           round(coords(3)))));
          end
          
          % Redraw plot window if open
          redrawPlotWindow();
          
          % Tidy up
          clearvars S C
          
      end
   end
   
   clearvars temp temphdr
   
end



%==========================================================================
% Go to nearest maximum
%==========================================================================
function goToMaximum(objectHandle, eventData)

    if size(maximaMat,1) == 0 || size(maximaMat,2) == 0
        msgbox('There are no maxima in the image at the given thresholding level.')
        return
    end
    
    handles = guidata(findobj('type','figure','name','MRM Post-estimation Tools'));
    temp    = handles.template;
    temphdr = spm_data_hdr_read(MRMoptions.Template);
    
    coords    = [OverlayMatCoords(1) OverlayMatCoords(2) OverlayMatCoords(3)];
    
    MNICoords    = temphdr.mat * [coords 1]';
    MNICoords(4) = [];
    MNICoords    = round(MNICoords');
    
    % Euclidean distance
    D = sqrt(((maximaMat(:,2) - MNICoords(1)).^2) + ...
             ((maximaMat(:,3) - MNICoords(2)).^2) + ...
             ((maximaMat(:,4) - MNICoords(3)).^2));
         
    if isempty(D) == 1
        msgbox('No maxima in the image')
        return
    end
    
    chs      = get(handles.crosshairsCheck, 'Value');
    infoPn   = get(handles.infoPaneBox, 'Value');
    
    maximaMNICoords   = maximaMat(D == min(D),2:4);
     
    maximaMatCoords    = [maximaMNICoords 1] * (inv(temphdr.mat))';
    maximaMatCoords(4) = [];
    maximaMatCoords    = round(maximaMatCoords);
    OverlayMatCoords   = maximaMatCoords;
    
    if get(handles.overlayPopup, 'Value') == 2
        threshold = 0;
    else
        threshold = str2num(get(handles.thresholdText, 'String'));
    end
    
    set(handles.XcoordText, 'String', num2str(maximaMNICoords(1)));
    set(handles.YcoordText, 'String', num2str(maximaMNICoords(2)));
    set(handles.ZcoordText, 'String', num2str(maximaMNICoords(3)));
    
    % Redraw saggital
    slice = squeeze(temp(round(maximaMatCoords(1)),:,:));
    S     = squeeze(handles.reslicedOverlay(round(maximaMatCoords(1)),:,:));
    S(S < threshold) = NaN;
    redrawSaggital(slice, [], S, maximaMatCoords, chs);
    
    % Redraw coronal
    slice = squeeze(temp(:, round(maximaMatCoords(2)), :));
    C     = squeeze(handles.reslicedOverlay(:,round(maximaMatCoords(2)),:));
    C(C < threshold) = NaN;
    redrawCoronal(slice, temphdr, C, maximaMatCoords, chs);
    
    % Redraw axial
    slice = squeeze(temp(:, :, round(maximaMatCoords(3))));
    A     = squeeze(handles.reslicedOverlay(:,:,round(maximaMatCoords(3))));
    A(A < threshold) = NaN;
    redrawAxial(slice, temphdr, A, maximaMatCoords, chs);
    
    % Redraw info pane
    if infoPn == 1
        redrawInfo([num2str(maximaMNICoords(1)) ' ' num2str(maximaMNICoords(2)) ' ' ...
                    num2str(maximaMNICoords(3))], ...
                    num2str(handles.reslicedOverlay(round(maximaMatCoords(1)), ...
                                                    round(maximaMatCoords(2)), ...
                                                    round(maximaMatCoords(3)))));
    end
          
    
    % Update the table location
    if isempty(findobj('type','figure','name','Results Table')) ~= 1
        
        % See: http://uk.mathworks.com/matlabcentral/newsreader/view_thread/165066
        jscroll  = findjobj(findobj('type','uitable','tag','resultsTable'));
        h        = jscroll.getComponents;
        a        = h(1).getComponents;
        jtable   = a(1);
        row      = find(D == min(D));
        
        jtable.setRowSelectionAllowed(0);
        jtable.setColumnSelectionAllowed(0);
        jtable.changeSelection(row-1,0, false, false); % Java indexes from 0
        
    else
        % Redraw plot window if open
        redrawPlotWindow();
    end
    
    clearvars S C A
end


%==========================================================================
% Go to global maximum
%==========================================================================
function goToGlobalMaximum(objectHandle, eventData)

    if size(maximaMat,1) == 0 || size(maximaMat,2) == 0
        msgbox('There are no maxima in the image at the given thresholding level.')
        return
    end
    
    handles = guidata(findobj('type','figure','name','MRM Post-estimation Tools'));
    temp      = handles.template;
    vals      = maximaMat(:,1);
    MNICoords = maximaMat(vals == max(vals),2:4);
    
    if isempty(MNICoords) == 1
        msgbox('No maxima in the image')
        return
    end
    
    pos     = get(findobj('type','figure','name','MRMviewer'), 'Position');
    %handles = guidata(findobj('type','figure','name','MRM Post-estimation Tools'));
    chs     = get(handles.crosshairsCheck, 'Value');
    infoPn  = get(handles.infoPaneBox, 'Value');
    
    if get(handles.overlayPopup, 'Value') == 2
        threshold = 0;
    else
        threshold = str2num(get(handles.thresholdText, 'String'));
    end
    
    set(handles.XcoordText, 'String', num2str(MNICoords(1)));
    set(handles.YcoordText, 'String', num2str(MNICoords(2)));
    set(handles.ZcoordText, 'String', num2str(MNICoords(3)));
    
    MRM_overlayImageNewNew(handles.reslicedOverlay, MNICoords, threshold, ...
                           pos, maximaMat, chs, infoPn, temp)
    
    % Update the table location
    if isempty(findobj('type','figure','name','Results Table')) ~= 1
        
        % See: http://uk.mathworks.com/matlabcentral/newsreader/view_thread/165066
        
        jscroll  = findjobj(findobj('type','uitable','tag','resultsTable'));
        h        = jscroll.getComponents;
        a        = h(1).getComponents;
        jtable   = a(1);
        row      = find(vals == max(vals));
        
        jtable.setRowSelectionAllowed(0);
        jtable.setColumnSelectionAllowed(0);
        jtable.changeSelection(row-1,0, false, false); % Java indexes from 0
        
    else
       % Redraw plot window if open
       redrawPlotWindow();
    end
end


%==========================================================================
% Go to largest cluster
%==========================================================================
function goToLargestClust(objectHandle, eventData)

    if size(maximaMat,1) == 0 || size(maximaMat,2) == 0
        msgbox('There are no clusters in the image at the given thresholding level.')
        return
    end
    
    handles      = guidata(findobj('type','figure','name','MRM Post-estimation Tools'));
    temp         = handles.template;
    biggestClust = maximaMat(maximaMat == max(maximaMat(:,1)),:);
    MNICoords    = biggestClust(1,2:4);
    
    pos     = get(findobj('type','figure','name','MRMviewer'), 'Position');
    chs     = get(handles.crosshairsCheck, 'Value');
    infoPn  = get(handles.infoPaneBox, 'Value');
    
    if get(handles.overlayPopup, 'Value') == 2
        threshold = 0;
    else
        threshold = str2num(get(handles.thresholdText, 'String'));
    end
    
    set(handles.XcoordText, 'String', num2str(MNICoords(1)));
    set(handles.YcoordText, 'String', num2str(MNICoords(2)));
    set(handles.ZcoordText, 'String', num2str(MNICoords(3)));
    
    MRM_overlayImageNewNew(handles.reslicedOverlay, MNICoords, threshold, ...
                           pos, maximaMat, chs, infoPn, temp)
    
    % Update the table location
    if isempty(findobj('type','figure','name','Results Table')) ~= 1
        
        % See: http://uk.mathworks.com/matlabcentral/newsreader/view_thread/165066
        
        jscroll  = findjobj(findobj('type','uitable','tag','resultsTable'));
        h        = jscroll.getComponents;
        a        = h(1).getComponents;
        jtable   = a(1);
        row      = find(maximaMat == max(maximaMat(:,1)), 1);
        
        jtable.setRowSelectionAllowed(0);
        jtable.setColumnSelectionAllowed(0);
        jtable.changeSelection(row-1,0, false, false); % Java indexes from 0
        
    else
        % Redraw plot window if open
        redrawPlotWindow();
    end
end




%==========================================================================
% Redraw plot window
%==========================================================================
function redrawPlotWindow()
    
    if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1
        
        nDVs      = size(MRM.Design.Y.Y,2);
        coords    = [OverlayMatCoords(1) OverlayMatCoords(2) OverlayMatCoords(3)];
        MNICoords = MRM_mat2mni(coords, MRMoptions.Template);
        MatCoords = MRM_mni2mat(MNICoords, [MRM.Options.Out filesep 'MRM_PE_1_1.nii']);
        
        % Read in PEs
        %------------
        tempPE = zeros(size(MRM.Design.X.X,2),nDVs);
        
        for k = 1:size(MRM.Design.X.X,2)
            for j = 1:nDVs
                
                file = [MRM.Options.Out filesep 'MRM_PE_' num2str(k) '_' num2str(j) '.nii'];
                
                tempPE(k,j) = spm_data_read(spm_data_hdr_read(file), ...
                    'xyz', ...
                    [MatCoords(1) MatCoords(2) MatCoords(3)]');
            end
        end
        
        % Read in variance-covariance matrix
        %-----------------------------------
        tempVCOV = zeros(nDVs,nDVs);
        
        for k = 1:nDVs
            for j = 1:nDVs
                
                if k >= j
                    
                    file = [MRM.Options.Out filesep 'MRM_Covar_' num2str(k) '_' num2str(j) '.nii'];
                    
                    tempVCOV(k,j) = spm_data_read(spm_data_hdr_read(file), ...
                                                    'xyz', ...
                                                    [MatCoords(1) MatCoords(2) MatCoords(3)]');
                end
            end
        end
        
        VCOVlowerTri = tril(tempVCOV, -1);
        tempVCOV     = tempVCOV + VCOVlowerTri';
        
        hand    = findobj('type','figure','name','MRM Post-estimation Tools');
        handles = guidata(hand);
        
        selection = get(handles.conPlotSelectMenu,'Value');
        
        A = MRM.Contrasts.Plots.Con{selection - 1}.A;
        C = MRM.Contrasts.Plots.Con{selection - 1}.C;
        
        MRM_barPlot(A,                                           ...
            tempPE,                                      ...
            C,                                           ...
            tempVCOV,                                    ...
            get(handles.plotCIsButton, 'Value'),         ...
            MNICoords,                                   ...
            MRM.Contrasts.Plots.Con{selection - 1}.Name, ...
            0,                                           ...
            [],                                          ...
            str2double(get(handles.CILevel, 'String')) / 100);
    end
end


%==========================================================================
% Redraw coronal slice
%==========================================================================
function redrawCoronal(templateSlice, templateHDR, overlaySlice, coordinates, crosshairs)
   
    cla(sliceHandles(2));
    
    subplot(2,2,1);
    imagesc(fliplr(rot90(templateSlice,3)));
    colormap gray

    if sign(templateHDR.mat(1,1)) == -1
        set(gca,'XDir','reverse');
    end

    hold on
    Img = subimage(fliplr(rot90(overlaySlice,3)),hot);
    set(Img, 'AlphaData', 1 - isnan(fliplr(rot90(overlaySlice,3))));
    set(gca,'YDir','normal');
    
    xLimits = get(gca, 'xlim');
    yLimits = get(gca, 'ylim');

    if crosshairs == 1

        line('XData', [0 xLimits(2)],        ...
            'YData', [coordinates(3) coordinates(3)], ...
            'LineStyle', '-',               ...
            'LineWidth', 1,                 ...
            'Color', 'blue');

        line('XData', [coordinates(1) coordinates(1)], ...
            'YData',  [0 yLimits(2)],        ...
            'LineStyle', '-',                ...
            'LineWidth', 1,                  ...
            'Color', 'blue');
    end

    axis off
    hold off
end

%==========================================================================
% Redraw axial slice
%==========================================================================
function redrawAxial(templateSlice, templateHDR, overlaySlice, coordinates, crosshairs)
    
    cla(sliceHandles(3));
    
    subplot(2,2,3);
    imagesc(fliplr(rot90(templateSlice,3)));
    colormap gray

    if sign(templateHDR.mat(1,1)) == -1
        set(gca,'XDir','reverse');
    end
    
    hold on
    Img = subimage(fliplr(rot90(overlaySlice,3)),hot);
    set(Img, 'AlphaData', 1 - isnan(fliplr(rot90(overlaySlice,3))));
    set(gca,'YDir','normal');
    
    xLimits = get(gca, 'xlim');
    yLimits = get(gca, 'ylim'); 
    
    if crosshairs == 1

        line('XData', [0 xLimits(2)],                 ...
            'YData', [coordinates(2) coordinates(2)], ...
            'LineStyle', '-',                         ...
            'LineWidth', 1,                           ...
            'Color', 'blue');

        line('XData', [coordinates(1) coordinates(1)], ...
            'YData',  [0 yLimits(2)],                  ...
            'LineStyle', '-',                          ...
            'LineWidth', 1,                            ...
            'Color', 'blue');
    end
    
    axis off
    hold off
end


%==========================================================================
% Redraw saggital slice
%==========================================================================
function redrawSaggital(templateSlice, ~, overlaySlice, coordinates, crosshairs)
    
    cla(sliceHandles(1));
    
    subplot(2,2,2);
    imagesc(fliplr(rot90(templateSlice,3)));
    colormap gray
    
    hold on
    Img = subimage(fliplr(rot90(overlaySlice,3)),hot);
    set(Img, 'AlphaData', 1 - isnan(fliplr(rot90(overlaySlice,3))));
    set(gca,'YDir','normal');
    
    xLimits = get(gca, 'xlim');
    yLimits = get(gca, 'ylim'); 
    
    if crosshairs == 1

        line('XData', [0 xLimits(2)],                   ...                                                             ...
            'YData', [coordinates(3) coordinates(3)],   ...
            'LineStyle', '-',                           ...
            'LineWidth', 1,                             ...
            'Color', 'blue');

        line('XData', [coordinates(2) coordinates(2)],  ...
            'YData',  [0 yLimits(2)],                   ...
            'LineStyle', '-',                           ...
            'LineWidth', 1,                             ...
            'Color', 'blue');
    end
    
    axis off
    hold off
end



%==========================================================================
% Redraw info panel
%==========================================================================
function redrawInfo(coordsString, FvalString)

    cla(sliceHandles(4));
    
    subplot(2,2,4);
    set(gca, 'color', 'black')
    set(gca, 'units', 'normalized')
    
    hold on
    
    % Coordinates
    text(0.05, 0.6, ...
        ['Coordinates: ' coordsString], ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 18, ...
        'FontUnits', 'normalized', ...
        'Clipping', 'off', ...
        'Color', 'white');
    
    % F-value
    text(0.05, 0.5, ...
        ['F-value: ' FvalString], ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 18, ...
        'FontUnits', 'normalized', ...
        'Clipping', 'off', ...
        'Color', 'white');
    
    % Atlas label
    text(0.05, 0.4, ...
        ['Label: ' MRM_getAtlasLabel(MRMoptions.Atlas, str2num(coordsString))], ...
        'HorizontalAlignment', 'left',   ...
        'VerticalAlignment',   'top',    ...
        'FontSize',            18,       ...
        'Clipping',            'off',    ...
        'Color',               'white');
    
    hold off
end


%==========================================================================
% Save image
%==========================================================================
function saveImage(objectHandle, eventData)
    
    [FileName,PathName] = uiputfile('Image.jpg', 'Save Image as');
    
    hand = findobj('type','figure','name','MRMviewer');
    set(hand, 'InvertHardCopy', 'off');
    
    if ~isequal(FileName,0) || ~isequal(PathName,0)
        saveas(hand, [PathName FileName], 'jpeg')
    end

end



%==========================================================================
% Close window
%==========================================================================
function CloseWindowCallback(objectHandle, eventData)

    clearvars -global maximaMat OverlayMatCoords sliceHandles
    
    delete(findobj('type','figure','name','MRMviewer'));
    
end
 
end
