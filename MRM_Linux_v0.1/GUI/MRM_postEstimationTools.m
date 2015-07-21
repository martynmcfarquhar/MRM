%=========================================================================
% GUI for the post-estimation tools
%=========================================================================
% This script creates the GUI used to explore results after model and
% contrast estimation. This is a complicated GUI, reflected in the length 
% of the script. Be careful altering even small things here as a lot of the
% functionality is inter-dependent across multiple aspects of the GUI.
%
%=========================================================================
% Copyright 2015 Martyn McFarquhar                                        
%=========================================================================
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
%=========================================================================

function varargout = MRM_postEstimationTools(varargin)
%=========================================================================%
% Begin initialization code - DO NOT EDIT
%=========================================================================%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_postEstimationTools_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_postEstimationTools_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%=========================================================================%
% End initialization code - DO NOT EDIT
%=========================================================================%




%==========================================================================
% Opening function
%==========================================================================
function MRM_postEstimationTools_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM
global MRMoptions

% Choose default command line output for MRM_postEstimationTools
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Load options
MRM_getOptions();

%--------------------------------------------------------------------------
% Set the contrasts list
%--------------------------------------------------------------------------
if MRM.Contrasts.Number > 0
   
    conLabs = cell(MRM.Contrasts.Number + 1,1);
    
    conLabs{1} = 'Select contrast to explore ...';
    
    for i = 2:MRM.Contrasts.Number + 1   
        conLabs{i} = MRM.Contrasts.Con{i-1}.Name;  
    end
    
    set(handles.contrastExploreMenu, 'String', conLabs);
    
else
    
    set(handles.contrastExploreMenu, 'String', 'No contrasts defined');
    
end


%--------------------------------------------------------------------------
% Set the overlay list
%--------------------------------------------------------------------------
overlayLabs    = cell(3,1);
overlayLabs{1} = 'Overlay ...';
thresh         = MRM.Options.Thresh.PThresh;

switch MRM.Options.Thresh.Level  
    case 'None'
        overlayLabs{2} = 'Unthresholded F image';
        set(handles.thresholdText, 'Enable', 'On');
    
    case 'Uncorrected'  
        overlayLabs{2} = ['F image thresholded at uncorrected voxel p < ' num2str(thresh)];
        overlayLabs{3} = 'Unthresholded F image'; 
        
    case 'Voxel'
            switch MRM.Options.Thresh.Method 
                case 'FWE'
                    overlayLabs{2} = ['F image thresholded at FWE voxel p < ' num2str(thresh)];
                    overlayLabs{3} = 'Unthresholded F image';
                    
                case 'FDR'    
                    overlayLabs{2} = ['F image thresholded at FDR voxel q < ' num2str(thresh)];
                    overlayLabs{3} = 'Unthresholded F image'; 
            end
    case 'Cluster'  
        overlayLabs{2} = ['F image thresholded at FWE cluster p < ' num2str(thresh)];
        overlayLabs{3} = 'Cluster-forming thresholded F image';    
end

set(handles.overlayPopup, 'String', overlayLabs);
set(handles.overlayPopup, 'Enable', 'off');

%--------------------------------------------------------------------------
% Set the DV list
%--------------------------------------------------------------------------
DVlabs    = cell(size(MRM.Design.Y.Cell,2),1);
DVlabs{1} = 'All DVs';

if size(MRM.Design.Y.Cell,2) > 1
    ind = 1;
    for i = 2:size(MRM.Design.Y.Cell,2) + 1
        DVlabs{i} = MRM.Design.Y.Cell{ind}.Label;
        ind = ind + 1;
    end
else
    set(handles.DVselectMenu, 'Enable', 'Off');
end

set(handles.DVselectMenu, 'String', DVlabs);


%--------------------------------------------------------------------------
% Update text if contrasts have been previously specified
%--------------------------------------------------------------------------
if isfield(MRM.Contrasts, 'Plots')
    if MRM.Contrasts.Plots.Number > 0

        set(handles.numConsText, 'String', [num2str(MRM.Contrasts.Plots.Number) ...
            ' contrast(s) specified']);

        % Update the drop-down menu
        conNames = cell(MRM.Contrasts.Plots.Number + 1,1);

        conNames{1} = 'Select contrast for plotting ...';

        ind = 1;

        for i = 2:MRM.Contrasts.Plots.Number + 1
            conNames{i} = MRM.Contrasts.Plots.Con{ind}.Name;
            ind = ind + 1;
        end

        set(handles.conPlotSelectMenu, 'String', conNames);
    end
end


%--------------------------------------------------------------------------
% Add PlotTools and Utilities to MATLAB path
%--------------------------------------------------------------------------

path = fileparts(which('MRM_launcher.m'));
addpath(genpath([path filesep 'PlotTools']));
addpath(genpath([path filesep 'Utilities']));

%--------------------------------------------------------------------------
% Disable some of the menus if we are doing cluster-based inference
%--------------------------------------------------------------------------
if strcmp(MRM.Options.Thresh.Level, 'Cluster') == 1
   
    set(handles.DVselectMenu,     'Enable', 'off');
    set(handles.whichPlotMenu,    'Enable', 'off');
    set(handles.drawPlotsButton,  'Enable', 'off');
    set(handles.assumpTestsMenu,  'Enable', 'off');
    set(handles.assumpTestButton, 'Enable', 'off');
    
end




% --- Outputs from this function are returned to the command line.
function varargout = MRM_postEstimationTools_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



%==========================================================================
% X coordinate text box
%==========================================================================
function XcoordText_Callback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Check a number was entered, if it was round it to a whole number to
% prevent fractional coordinates
%--------------------------------------------------------------------------
if isnan(str2double(get(hObject,'String')))
    
    msgbox('Coordinates must be numeric');
    set(hObject, 'String', '0');
    
else
    
    num = str2double(get(hObject,'String'));
    num = round(num);
    set(hObject, 'String', num);
    
end

% Check whether a viewer window is open. If it is refresh with the new
% coordinates
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    crosshairs = get(handles.crosshairsCheck, 'Value');
    infoPane   = get(handles.infoPaneBox, 'Value');
    hand       = findobj('type','figure','name','MRMviewer');
    pos        = get(hand, 'Position');
    
    coords  = [str2num(get(handles.XcoordText, 'String')), ...
               str2num(get(handles.YcoordText, 'String')), ...
               str2num(get(handles.ZcoordText, 'String'))];
    
    if get(handles.overlayPopup, 'Value') == 2   
        thresh = 0;
    else     
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    handles = guidata(hObject);
    
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
        maximaMat  = MRM_getMaxima();
    else
        maximaMat = [];
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    try
        MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                               maximaMat, crosshairs, infoPane, handles.template)
    end  
end

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    plotConButton_Callback(handles.plotConButton, eventdata, handles);
    
end



%==========================================================================
% Y coordinate text box
%==========================================================================
function YcoordText_Callback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Check a number was entered, if it was round it to a whole number to
% prevent fractional coordinates
%--------------------------------------------------------------------------
if isnan(str2double(get(hObject,'String')))
    
    msgbox('Coordinates must be numeric');
    set(hObject, 'String', '0');
    
else
    
    num = str2double(get(hObject,'String'));
    num = round(num);
    set(hObject, 'String', num);
    
end

% Check whether a viewer window is open. If it is refresh with the new
% coordinates
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    crosshairs = get(handles.crosshairsCheck, 'Value');
    infoPane   = get(handles.infoPaneBox, 'Value');
    
    hand    = findobj('type','figure','name','MRMviewer');
    pos     = get(hand, 'Position');
    
    coords  = [str2num(get(handles.XcoordText, 'String')), ...
               str2num(get(handles.YcoordText, 'String')), ...
               str2num(get(handles.ZcoordText, 'String'))];
           
    if get(handles.overlayPopup, 'Value') == 2   
        thresh = 0;
    else     
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    handles    = guidata(hObject);
    
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
        maximaMat  = MRM_getMaxima();
    else
        maximaMat = [];
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    try
        MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                               maximaMat, crosshairs, infoPane, handles.template)
    end
    
    
end

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    plotConButton_Callback(handles.plotConButton, eventdata, handles);
    
end



%==========================================================================
% Z coordinate text box
%==========================================================================
function ZcoordText_Callback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Check a number was entered, if it was round it to a whole number to
% prevent fractional coordinates
%--------------------------------------------------------------------------
if isnan(str2double(get(hObject,'String')))
    
    msgbox('Coordinates must be numeric');
    set(hObject, 'String', '0');
    
else
    
    num = str2double(get(hObject,'String'));
    num = round(num);
    set(hObject, 'String', num);
    
end

% Check whether a viewer window is open. If it is refresh with the new
% coordinates
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    crosshairs = get(handles.crosshairsCheck, 'Value');
    infoPane   = get(handles.infoPaneBox, 'Value');
    hand       = findobj('type','figure','name','MRMviewer');
    pos        = get(hand, 'Position');
    
    coords     = [str2num(get(handles.XcoordText, 'String')), ...
                  str2num(get(handles.YcoordText, 'String')), ...
                  str2num(get(handles.ZcoordText, 'String'))];
           
    if get(handles.overlayPopup, 'Value') == 2   
        thresh = 0;
    else     
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    if strcmp(MRM.Options.Thresh.Level, 'None') == 1
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    handles = guidata(hObject);
    
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
        maximaMat  = MRM_getMaxima();
    else
        maximaMat = [];
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    try
        MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                               maximaMat, crosshairs, infoPane, handles.template)
    end
    
end

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    plotConButton_Callback(handles.plotConButton, eventdata, handles);
    
end



%==========================================================================
% 'Draw the selected plot' button
%==========================================================================
function drawPlotsButton_Callback(hObject, eventdata, handles)

global MRM

nSubs = size(MRM.Design.Y.Y,1);
nDVs  = size(MRM.Design.Y.Y,2);
cells = size(MRM.Design.X.Cell,2);

whichPlot = get(handles.whichPlotMenu, 'Value');

if whichPlot == 5
    
    if license('test', 'statistics_toolbox') ~= 1
        msgbox('The statistics toolbox is required to create boxplots.');
        return
    end
    
elseif whichPlot == 3
    
    if license('test', 'statistics_toolbox') ~= 1
        msgbox('The statistics toolbox is required to create the normal Q-Q plots.');
        return
    end
    
elseif whichPlot == 6
    
    if get(handles.DVselectMenu, 'Value') > 1
        msgbox('Cannot draw scatterplot pairs for only 1 DV.');
        return
        
    elseif (nDVs^2 + nDVs)/2 > 28
        msgbox(['Your covariance matrix has ' num2str((nDVs^2 + nDVs)/2)       ...
                ' unique elements, which is greater than the maximum number allowed ' ...
                'for drawing the scatterplot pairs.']);
        return
    end
end

%--------------------------------------------------------------------------
% Gather file paths
%--------------------------------------------------------------------------
files = cell(1,nDVs);

for i = 1:nDVs  
    for j = 1:cells
        files{i} = [files{i} MRM.Data.Y{i}.Cell{j}.Scans];
    end
end


%--------------------------------------------------------------------------
% Find the matrix coordinates
%--------------------------------------------------------------------------
MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
             str2double(get(handles.YcoordText, 'String'))  ...
             str2double(get(handles.ZcoordText, 'String'))];
  
try         
    [matCoords, ~] = MRM_mni2mat(MNIcoords, MRM.Data.Y{1}.Cell{1}.Scans{1});
catch
    msgbox('Cannot find raw data.')
    return
end


%--------------------------------------------------------------------------
% Check mask
%--------------------------------------------------------------------------
if isempty(MRM.Options.Mask) ~= 1

    maskPath      = [MRM.Options.Out filesep 'MRM_mask.nii'];
    
    maskImg = spm_data_read(spm_data_hdr_read(maskPath), ...
                            'xyz', ...
                            [matCoords(1) matCoords(2) matCoords(3)]');

    if maskImg == 0 || isnan(maskImg) == 1
        msgbox('The selected voxel is outside of the masked region');
        return;
    end

    clearvars maskImg
end


%--------------------------------------------------------------------------
% Get data from coordinate
%--------------------------------------------------------------------------
voxData = zeros(nSubs, nDVs);

for i = 1:nDVs
    for j = 1:nSubs
        try
            voxData(j,i) = spm_data_read(spm_data_hdr_read(files{i}{j}),            ...
                                         'xyz',                                     ...
                                         [matCoords(1) matCoords(2) matCoords(3)]');
        catch
            msgbox('Cannot find raw data.')
            return
        end
    end 
end

% Create plots
%-------------
if size(MRM.Design.Y.Cell,2) > 1
    which = get(handles.DVselectMenu, 'Value');
else
    which = 2;
end

if which == 1
    MRM_assumptionPlots(voxData, MNIcoords, whichPlot, 0);
else
    MRM_assumptionPlots(voxData(:,which - 1), MNIcoords, whichPlot, which - 1);
end
    

%-Tidy up
%----------------------
clearvars voxData


%==========================================================================
% 'Save model to text file' button
%==========================================================================
function saveModOutcomesButton_Callback(hObject, eventdata, handles)

global MRM

nSubs = size(MRM.Design.Y.Y,1);
nDVs  = size(MRM.Design.Y.Y,2);
cells = size(MRM.Design.X.Cell,2);
out   = MRM.Options.Out;


%--------------------------------------------------------------------------
% Gather file paths
%--------------------------------------------------------------------------
files = cell(1,nDVs);

for i = 1:nDVs  
    for j = 1:cells
        files{i} = [files{i} MRM.Data.Y{i}.Cell{j}.Scans];
    end
end


%--------------------------------------------------------------------------
% Find the matrix coordinates
%--------------------------------------------------------------------------
MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
             str2double(get(handles.YcoordText, 'String'))  ...
             str2double(get(handles.ZcoordText, 'String'))];

try
    [matCoords, ~] = MRM_mni2mat(MNIcoords, MRM.Data.Y{1}.Cell{1}.Scans{1});
catch
    msgbox('Cannot find raw data.');
    return
end


%--------------------------------------------------------------------------
% Check mask
%--------------------------------------------------------------------------
if isempty(MRM.Options.Mask) ~= 1

    maskPath      = [MRM.Options.Out filesep 'MRM_mask.nii'];
    
    maskImg = spm_data_read(spm_data_hdr_read(maskPath), ...
                            'xyz', ...
                            [matCoords(1) matCoords(2) matCoords(3)]');

    if maskImg == 0 || isnan(maskImg) == 1
        msgbox('The selected voxel is outside of the masked region');
        return;
    end

    clearvars maskImg
end


%--------------------------------------------------------------------------
% Output directory
%--------------------------------------------------------------------------
saveDir = uigetdir(pwd, 'Select output directory');


%--------------------------------------------------------------------------
% Get data from coordinate
%--------------------------------------------------------------------------
voxData = zeros(nSubs, nDVs);

for i = 1:nDVs
    for j = 1:nSubs
        try
            voxData(j,i) = spm_data_read(spm_data_hdr_read(files{i}{j}),            ...
                                         'xyz',                                     ...
                                         [matCoords(1) matCoords(2) matCoords(3)]');
        catch
            msgbox('Cannot find raw data.')
            return
        end
    end 
end

% Save to file
fileID = fopen([saveDir filesep 'RawData_' num2str(MNIcoords(1))             ...
                '_' num2str(MNIcoords(2)) '_' num2str(MNIcoords(3)) '.txt'], ...
                'w');

if fileID == -1
    return
end
            
if ispc == 1
    fprintf(fileID, [repmat('%.15f\t', 1, size(voxData,2)) '\r\n'], voxData');
else
    fprintf(fileID, [repmat('%.15f\t', 1, size(voxData,2)) '\n'], voxData');
end

fclose(fileID);


%--------------------------------------------------------------------------
% Read in PEs
%--------------------------------------------------------------------------
tempPE = zeros(size(MRM.Design.X.X,2), nDVs);

for i = 1:size(MRM.Design.X.X,2)
    for j = 1:nDVs

        tempPE(i,j) = spm_data_read(spm_data_hdr_read([                 ...
                                    out filesep 'MRM_PE_'               ...
                                    num2str(i) '_' num2str(j) '.nii']), ...
                                    'xyz',                              ...
                                    [matCoords(1) matCoords(2) matCoords(3)]');
    end
end

% Predicted values
predVals = MRM.Design.X.X * tempPE;

% Save to file
fileID = fopen([saveDir filesep 'PredictedValues_' num2str(MNIcoords(1))     ...
                '_' num2str(MNIcoords(2)) '_' num2str(MNIcoords(3)) '.txt'], ...
                'w');
if ispc == 1
    fprintf(fileID, [repmat('%.15f\t', 1, size(predVals,2)) '\r\n'], predVals');
else
    fprintf(fileID, [repmat('%.15f\t', 1, size(predVals,2)) '\n'], predVals');
end

fclose(fileID);


%--------------------------------------------------------------------------
% Residuals
%--------------------------------------------------------------------------
residVals = voxData - predVals;

% Save to file
fileID = fopen([saveDir filesep 'Residuals_' num2str(MNIcoords(1))           ...
                '_' num2str(MNIcoords(2)) '_' num2str(MNIcoords(3)) '.txt'], ...
                'w');
if ispc == 1
    fprintf(fileID, [repmat('%.15f\t', 1, size(residVals,2)) '\r\n'], residVals');
else
    fprintf(fileID, [repmat('%.15f\t', 1, size(residVals,2)) '\n'], residVals');
end

fclose(fileID);


%--------------------------------------------------------------------------
% Tidy up
%--------------------------------------------------------------------------
clearvars voxData predVals residVals data;




%==========================================================================
% 'View Design' button
%==========================================================================
function viewDesignButton_Callback(hObject, eventdata, handles)

MRM_drawDesign(0);




%==========================================================================
% 'Specify contrasts for plotting' button
%==========================================================================
function specifyPlotConsButton_Callback(hObject, eventdata, handles)

global MRM

if isfield(MRM.Contrasts, 'Plots') ~= 1
    MRM.Contrasts.Plots.Number = 0;
end

h = MRM_contrastsSelectorPlotting();
uiwait(h);

% Set the number of contrasts
set(handles.numConsText, 'String', [num2str(MRM.Contrasts.Plots.Number) ...
                                    ' contrast(s) specified']);

% Update the drop-down menu
set(handles.conPlotSelectMenu, 'Value', 1);

conNames = cell(MRM.Contrasts.Plots.Number + 1,1);

conNames{1} = 'Select contrast for plotting ...';

ind = 1;
    
for i = 2:MRM.Contrasts.Plots.Number + 1      
    conNames{i} = MRM.Contrasts.Plots.Con{ind}.Name;
    ind = ind + 1;       
end

set(handles.conPlotSelectMenu, 'String', conNames);



%==========================================================================
% 'Standard error' checkbox
%==========================================================================
function standardErrorCheck_Callback(hObject, eventdata, handles)
 
if get(hObject, 'Value') == 1
    set(handles.plotCIsButton, 'Value', 0); 
    set(handles.CILevel, 'Enable', 'Off');
else
    set(handles.plotCIsButton, 'Value', 1);
    set(handles.CILevel, 'Enable', 'On');
end

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    plotConButton_Callback(handles.plotConButton, eventdata, handles);
    
end
    


%==========================================================================
% '95% CIs' checkbox
%==========================================================================
function plotCIsButton_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.standardErrorCheck, 'Value', 0);
    set(handles.CILevel, 'Enable', 'On');
else
    set(handles.standardErrorCheck, 'Value', 1);
    set(handles.CILevel, 'Enable', 'Off');
end

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    plotConButton_Callback(handles.plotConButton, eventdata, handles);
    
end


%==========================================================================
% Contrast selection menu
%==========================================================================
function conPlotSelectMenu_Callback(hObject, eventdata, handles)

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    plotConButton_Callback(handles.plotConButton, eventdata, handles);
    
end



%==========================================================================
% 'Plot contrasts' button
%==========================================================================
function plotConButton_Callback(hObject, eventdata, handles)

global MRM

selection = get(handles.conPlotSelectMenu,'Value');

if selection == 1
    return
end
    
out   = MRM.Options.Out;
nDVs  = size(MRM.Design.Y.Y,2);

%--------------------------------------------------------------------------
% Find the matrix coordinates
%--------------------------------------------------------------------------
MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
             str2double(get(handles.YcoordText, 'String'))  ...
             str2double(get(handles.ZcoordText, 'String'))];

[matCoords, ~] = MRM_mni2mat(MNIcoords, [out filesep 'MRM_PE_1_1.nii']);


%--------------------------------------------------------------------------
% Read in PEs
%--------------------------------------------------------------------------
tempPE = zeros(size(MRM.Design.X.X,2),nDVs);

for i = 1:size(MRM.Design.X.X,2)
    for j = 1:nDVs
        
        file = [out filesep 'MRM_PE_' num2str(i) '_' num2str(j) '.nii'];
        
        tempPE(i,j) = spm_data_read(spm_data_hdr_read(file), ...
                                    'xyz', ...
                                    [matCoords(1) matCoords(2) matCoords(3)]');
    end
end

%--------------------------------------------------------------------------
% Read in variance-covariance matrix
%--------------------------------------------------------------------------
tempVCOV = zeros(nDVs,nDVs);

for i = 1:nDVs
    for j = 1:nDVs
        
        if i >= j
        
            file = [out filesep 'MRM_Covar_' num2str(i) '_' num2str(j) '.nii'];

            tempVCOV(i,j) = spm_data_read(spm_data_hdr_read(file), ...
                                          'xyz', ...
                                          [matCoords(1) matCoords(2) matCoords(3)]');
                                  
        end
    end
end

VCOVlowerTri = tril(tempVCOV, -1);
tempVCOV     = tempVCOV + VCOVlowerTri';



%--------------------------------------------------------------------------
% Get contrast
%--------------------------------------------------------------------------
A = MRM.Contrasts.Plots.Con{selection - 1}.A;
C = MRM.Contrasts.Plots.Con{selection - 1}.C;

%--------------------------------------------------------------------------
% Create plot
%--------------------------------------------------------------------------
MRM_barPlot(A, tempPE, C, tempVCOV, get(handles.plotCIsButton, 'Value'), ...
            MNIcoords, MRM.Contrasts.Plots.Con{selection - 1}.Name,      ...
            0, [], str2double(get(handles.CILevel, 'String')) / 100);
    


%==========================================================================
% 'Help' button
%==========================================================================
function helpButton_Callback(hObject, eventdata, handles)

% Find the scripts directory
scriptDir = fileparts(mfilename('fullpath'));

% Find out if we are in Windows. If we are use the 'winopen' command,
% otherwise use a 'system' call to open the PDF
os = computer;

if strfind(os, 'WIN') ~= 0
    winopen([scriptDir filesep 'Help' filesep 'Manual.pdf'])
elseif strfind(os, 'MAC') ~= 0
    system(['open ' scriptDir filesep 'Help' filesep 'Manual.pdf']);
else
    system(['xdg-open ' scriptDir filesep 'Help' filesep 'Manual.pdf']);
end


%==========================================================================
% 'Exit' button
%==========================================================================
function exitButton_Callback(hObject, eventdata, handles)

figure1_CloseRequestFcn(handles.figure1, eventdata, handles)



%==========================================================================
% 'Save plot values' button
%==========================================================================
function savePlotValsButton_Callback(hObject, eventdata, handles)

global MRM

selection = get(handles.conPlotSelectMenu,'Value');

if selection == 1
    return
end

handles = guidata(hObject);
saveDir = uigetdir(pwd, 'Select output directory');   
out     = MRM.Options.Out;
nDVs    = size(MRM.Design.Y.Y,2);


%--------------------------------------------------------------------------
% Find the matrix coordinates
%--------------------------------------------------------------------------
MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
             str2double(get(handles.YcoordText, 'String'))  ...
             str2double(get(handles.ZcoordText, 'String'))];

[matCoords, ~] = MRM_mni2mat(MNIcoords, [out filesep 'MRM_PE_1_1.nii']);


%--------------------------------------------------------------------------
% Read in PEs
%--------------------------------------------------------------------------
tempPE = zeros(size(MRM.Design.X.X,2),nDVs);

for i = 1:size(MRM.Design.X.X,2)
    for j = 1:nDVs
        
        file = [out filesep 'MRM_PE_' num2str(i) '_' num2str(j) '.nii'];
        
        tempPE(i,j) = spm_data_read(spm_data_hdr_read(file), ...
            'xyz', ...
            [matCoords(1) matCoords(2) matCoords(3)]');
    end
end

%--------------------------------------------------------------------------
% Read in variance-covariance matrix
%--------------------------------------------------------------------------
tempVCOV = zeros(nDVs,nDVs);

for i = 1:nDVs
    for j = 1:nDVs
        
        file = [out filesep 'MRM_Covar_' num2str(i) '_' num2str(j) '.nii'];
        
        tempVCOV(i,j) = spm_data_read(spm_data_hdr_read(file), ...
            'xyz', ...
            [matCoords(1) matCoords(2) matCoords(3)]');
    end
end


%--------------------------------------------------------------------------
% Get contrast
%--------------------------------------------------------------------------
A = MRM.Contrasts.Plots.Con{selection - 1}.A;
C = MRM.Contrasts.Plots.Con{selection - 1}.C;

%--------------------------------------------------------------------------
% Create plot
%--------------------------------------------------------------------------
MRM_barPlot(A, tempPE, C, tempVCOV, get(handles.plotCIsButton, 'Value'), ...
            MNIcoords, MRM.Contrasts.Plots.Con{selection - 1}.Name,      ...
            1, saveDir, str2double(get(handles.CILevel, 'String')) / 100);

clearvars tempPE slice tempVCOV
 


%==========================================================================
% 'CIs' level textbox
%==========================================================================
function CILevel_Callback(hObject, eventdata, handles)

val = str2double(get(hObject, 'String'));

if isnan(val)
    set(hObject, 'String', '95'); 
elseif val >= 100
    set(hObject, 'String', '99');
end

% If a plot window already exists then just update it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1
    plotConButton_Callback(handles.plotConButton, eventdata, handles);   
end




%==========================================================================
% 'Display' button
%==========================================================================
function displayButton_Callback(hObject, eventdata, handles)

global MRM
global MRMoptions

if get(handles.overlayPopup, 'Value') == 1
    return
end

if strcmp(MRM.Options.Thresh.Level, 'None') == 1 && get(handles.overlayPopup, 'Value') == 3
    return
end

handles = guidata(hObject);
con     = get(handles.contrastExploreMenu, 'Value') - 1;
name    = regexprep(MRM.Contrasts.Con{con}.Name, ' ', '_');
start   = [MRM.Options.Out filesep];
f       = filesep;

% Get the path to the overlay
if get(handles.overlayPopup, 'Value') == 1
    return
    
elseif get(handles.overlayPopup, 'Value') == 2
    
    switch MRM.Options.Thresh.Level
        case 'None'
            switch MRM.Options.Thresh.Pvals 
                case 'Permutation'
                    handles.currentOverlay = [start 'PermutationContrasts' f name f ...
                                                'MRM_' name '_refF.nii'];
                case 'Approximate'
                    handles.currentOverlay = [start 'Contrasts' f name f ...
                                                'MRM_' name '_F.nii'];
            end
        case 'Uncorrected'
            switch MRM.Options.Thresh.Pvals 
                case 'Permutation'
                    handles.currentOverlay = [start 'PermutationContrasts' f name f ...
                                              'MRM_' name '_F_Thresholded_UCperm.nii']; 
                case 'Approximate'
                    handles.currentOverlay = [start 'Contrasts' f name f ...
                                              'MRM_' name '_F_Thresholded_UC.nii'];
            end
        case 'Voxel' 
            switch MRM.Options.Thresh.Pvals
                case 'Permutation'
                    switch MRM.Options.Thresh.Method
                        case 'FWE'
                            handles.currentOverlay = [start 'PermutationContrasts' f name f ...
                                                      'MRM_' name '_F_FWER_thresholded.nii'];
                        case 'FDR'
                            handles.currentOverlay = [start 'PermutationContrasts' f name f ...
                                                      'MRM_' name '_F_FDRperm_thresholded.nii'];
                    end
                case 'Approximate'
                    handles.currentOverlay = [start 'Contrasts' f name f ...
                                              'MRM_' name '_F_Thresholded_FDR.nii'];
            end
        case 'Cluster'
            handles.currentOverlay = [start 'PermutationContrasts' f name f ...
                                      'MRM_' name '_F_FWER_cluster_thresholded.nii'];
    end
    
    thresh = 0;
    
else
    
    switch MRM.Options.Thresh.Pvals
        case 'Permutation'
            handles.currentOverlay = [start 'PermutationContrasts' f name f ...
                'MRM_' name '_refF.nii'];
        case 'Approximate'
            handles.currentOverlay = [start 'Contrasts' f name f ...
                'MRM_' name '_F.nii'];
    end
    
    thresh  =  str2num(get(handles.thresholdText, 'String'));
    
end

overlayPath = fileparts(handles.currentOverlay);

if exist(overlayPath,'dir') == 0 
   msgbox('There were no results for the selected contrast');
   return
end

% Get list of maxima coordinates
if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
    maximaMat = MRM_getMaxima();   
else
    maximaMat = [];
    thresh    =  str2num(get(handles.thresholdText, 'String'));
end

coords  = [str2num(get(handles.XcoordText,    'String')), ...
           str2num(get(handles.YcoordText,    'String')), ...
           str2num(get(handles.ZcoordText,    'String'))];


%--------------------------------------------------------------------------
% Reslice the selected overlay to the template and save the resliced image
% to the handles structure
%--------------------------------------------------------------------------
% Thanks to Chris Rorden for providing help with the following code
%--------------------------------------------------------------------------
templatePath    = MRMoptions.Template;
tarhdr          = spm_data_hdr_read(templatePath);
inhdr           = spm_data_hdr_read(handles.currentOverlay);
inimg           = spm_read_vols(inhdr);
imgdim          = tarhdr.dim(1:3);

outhdr            = inhdr;
outhdr.dim(1:3)   = imgdim(1:3);
outhdr.mat        = tarhdr.mat;

overlay = NaN(outhdr.dim(1:3));

for i = 1:imgdim(3)
    
    M = inv(spm_matrix([0 0 -i]) / (outhdr.mat) * inhdr.mat);
    
    % Nearest neighbour interpolation
    overlay(:,:,i) = spm_slice_vol(inimg, M, imgdim(1:2), 0);
end

% MATLAB pads the edges with zeros, so we replace these with NaN
overlay(overlay == 0) = NaN;

% Put resliced overlay and template into handles so we can recall them quickly
handles.reslicedOverlay = overlay;
handles.template        = spm_data_read(tarhdr);
      
win        = get(0, 'ScreenSize');
crosshairs = get(handles.crosshairsCheck, 'Value');
infoPane   = get(handles.infoPaneBox, 'Value');
             
MRM_overlayImageNewNew(overlay, coords, thresh, [win(3)/4, win(4)/4,    ...
                       win(3)/2, win(4)/1.6], maximaMat, crosshairs,    ...
                       infoPane, handles.template);
% Update handles                                     
guidata(hObject, handles);

% Tidy-up
clearvars overlay M outhdr tarhdr
                            

%==========================================================================
% Threshold text box
%==========================================================================
function thresholdText_Callback(hObject, eventdata, handles)

global MRM

% Check whether a viewer window is open. If it is refresh with the new
% coordinates
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    hand      = findobj('type','figure','name','MRMviewer');
    pos       = get(hand, 'Position');
    coords    = [str2num(get(handles.XcoordText, 'String')), ...
                 str2num(get(handles.YcoordText, 'String')), ...
                 str2num(get(handles.ZcoordText, 'String'))];      
    thresh    =  str2num(get(handles.thresholdText, 'String'));
    handles   = guidata(hObject);
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
        maximaMat  = MRM_getMaxima();
    else
        maximaMat = [];
    end
    crosshairs = get(handles.crosshairsCheck, 'Value');
    infoPane   = get(handles.infoPaneBox, 'Value');
    
    MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                           maximaMat, crosshairs, infoPane, handles.template)
    
end


%==========================================================================
% Results table button
%==========================================================================
function resultsTableButton_Callback(hObject, eventdata, handles)

MRM_displayResultsTable(get(handles.contrastExploreMenu, 'Value') - 1);



%==========================================================================
% Close button
%==========================================================================
function figure1_CloseRequestFcn(hObject, eventdata, handles)

% If results table is open then close it
if isempty(findobj('type','figure','name','Results Table')) ~= 1
    
    h = findobj('type','figure','name','Results Table');
    delete(h);
    
end

% If viewer window is open then close it
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    h = findobj('type','figure','name','MRMviewer');
    delete(h);
    
    handles = guidata(hObject);
    handles.reslicedOverlay = [];
    guidata(hObject, handles);
    
end

% If plot window is open then close it
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1
    
    h = findobj('type','figure','name','MRM Bar Plot');
    delete(h);
    
end

% If dLDA window is open then close it
if isempty(findobj('type','figure','name','MRM dLDA')) ~= 1
    
    h = findobj('type','figure','name','MRM dLDA');
    delete(h);
    
end


delete(hObject);



%==========================================================================
% 'Select contrast to explore' menu
%==========================================================================
function contrastExploreMenu_Callback(hObject, eventdata, handles)

global MRM

if get(hObject, 'Value') == 1
    
    set(handles.resultsTableButton, 'Enable', 'off');
    set(handles.overlayPopup,       'Enable', 'off');
    set(handles.displayButton,      'Enable', 'off');
    set(handles.thresholdText,      'Enable', 'off');
    
else
    
    set(handles.resultsTableButton, 'Enable', 'on');
    set(handles.overlayPopup,       'Enable', 'on');
    set(handles.displayButton,      'Enable', 'on');
    
    if get(handles.overlayPopup, 'Value') == 3
        set(handles.thresholdText, 'Enable', 'on');
    else
        set(handles.thresholdText, 'Enable', 'off');
    end
    
    if strcmp(MRM.Options.Thresh.Level, 'None') == 1
        set(handles.thresholdText, 'Enable', 'on');
    end
    
    % Refresh the results table if it's open
    
    if isempty(findobj('type','figure','name','Results Table')) ~= 1
        
        h = findobj('type','figure','name','Results Table');
        set(h, 'Units', 'pixels');
        pos = get(h, 'Position');
        delete(h);
        MRM_displayResultsTable(get(handles.contrastExploreMenu, 'Value') - 1, pos);
        
    end
    
    % Refresh the overlay image if it's open
    if isempty(findobj('type','figure','name','MRMviewer')) ~= 1  
       displayButton_Callback(handles.displayButton, eventdata, handles);  
    end
    
end
    


%==========================================================================
% Overlay results menu
%==========================================================================
function overlayPopup_Callback(hObject, eventdata, handles)

global MRM

if get(hObject, 'Value') == 1  
    set(handles.thresholdText, 'Enable', 'off');    
elseif get(handles.overlayPopup, 'Value') == 3   
    set(handles.thresholdText, 'Enable', 'on');
elseif strcmp(MRM.Options.Thresh.Level, 'None') == 1
    set(handles.thresholdText, 'Enable', 'on');
else
    set(handles.thresholdText, 'Enable', 'off');
end

% If an image is open then change it
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1  
   displayButton_Callback(handles.displayButton, eventdata, handles);  
end



%==========================================================================
% 'Run test' button for the assumptions tests
%==========================================================================
function assumpTestButton_Callback(hObject, eventdata, handles)

global MRM

nSubs = size(MRM.Design.Y.Y,1);
nDVs  = size(MRM.Design.Y.Y,2);
cells = size(MRM.Design.X.Cell,2);

%-Load the data
%--------------
files = cell(1,nDVs);

for i = 1:nDVs  
    for j = 1:cells
        files{i} = [files{i} MRM.Data.Y{i}.Cell{j}.Scans];
    end
end

% Find the matrix coordinates
MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
             str2double(get(handles.YcoordText, 'String'))  ...
             str2double(get(handles.ZcoordText, 'String'))];

try         
    [matCoords, ~] = MRM_mni2mat(MNIcoords, MRM.Data.Y{1}.Cell{1}.Scans{1});
catch
    msgbox('Some of the raw data appears to have been moved. Cannot compute test.')
    return
end


% Check mask
if isempty(MRM.Options.Mask) ~= 1

    maskPath      = [MRM.Options.Out filesep 'MRM_mask.nii'];

    maskImg       = spm_data_read(spm_data_hdr_read(maskPath), ...
                                  'xyz', ...
                                  [matCoords(1) matCoords(2) matCoords(3)]');

    if maskImg == 0 || isnan(maskImg) == 1
        msgbox('The selected voxel is outside of the masked region');
        return;
    end

    clearvars maskImg
end


% Get data from coordinate 
voxData = zeros(nSubs, nDVs);

for i = 1:nDVs
    for j = 1:nSubs

        try
            voxData(j,i) = spm_data_read(spm_data_hdr_read(files{i}{j}), ...
                                        'xyz', ...
                                       [matCoords(1) matCoords(2) matCoords(3)]');
        catch
            msgbox('Cannot find the raw data.')
            return
        end
    end 
end
    
%--------------------------------------------------------------------------
% Box's M Test
%--------------------------------------------------------------------------    
if get(handles.assumpTestsMenu, 'Value') == 1
    
    % Run the test
    %----------------------------------------------------------
    X = MRM.Design.X.X(:,1:(size(MRM.Design.X.X,2) - MRM.Covariates.Number));
    
    [M, F, df1, df2, pval] = MRM_boxMtest(voxData, X);
    
    if isnan(M) ~= 1
        
        fprintf('%s\n', ...
        '----------------------------------------------------------------', ...
        'Box''s M Test for equality of covariance matrices'               , ...
        '----------------------------------------------------------------', ...
        '',                                                                 ...
        ['Results for voxel: ' num2str(MNIcoords)],                          ...
        '');
    
        displaytable([M F df1 df2 pval], {'Box''s M', 'F', 'df1', 'df2', 'p-value'})
        
        fprintf('\n');
        
        fprintf('%s\n', ...
                'Tests the null hypothesis that the population covariance ', ...
                'matrices are equal across the cells of the design.',        ...
                '',                                                          ...
                'NOTE: This test can be notoriously sensitive to departures',...
                'from normality. It is also usual to evaluate the test at',  ...
                'an alpha of 0.001.');

        fprintf('\n');
        
        msgbox('Results returned in the MATLAB command window');
    
    else
        msgbox('There was a problem computing the test. No output produced.')
    end

    % Tidy up
    clearvars voxData

%--------------------------------------------------------------------------
% Levene's Test
%--------------------------------------------------------------------------
else
    
    %----------------------------------------------------------------------
    % Run the test
    %----------------------------------------------------------------------
    [F, df1, df2, pval] = MRM_leveneTest(voxData, MRM.Design.X.X);
    
    fprintf('%s\n', ...
            '----------------------------------------------------------------', ...
            'Levene''s test for equality of variance',                          ...
            '----------------------------------------------------------------', ...
            '',                                                                 ...
            ['Results for voxel: ' num2str(MNIcoords)],                          ...
            '');
    
    rowLabs = cell(1,size(F,1));
    
    for i = 1:size(F,1)
        rowLabs{1,i} = MRM.Design.Y.Cell{i}.Label;
    end
        
    displaytable([F df1 df2 pval], {'F', 'df1', 'df2', 'p-value'}, [], [], rowLabs)
    
    fprintf('\n');
    
    fprintf('%s\n', ...                                              
            'Tests the null hypothesis that the population variances are equal', ...
            'across the cells of the design for each DV.',                       ...
            '',                                                                  ...
            'NOTE: This test uses the mean as the center of each group, and so', ...
            'can be sensitive to outliers. It is also usual to evaluate the',    ...
            'test at an alpha of 0.001.');
    
    fprintf('\n');
    
    msgbox('Results returned in the MATLAB command window');
    
    % Tidy up
    clearvars voxData
end


%==========================================================================
% View covariance structure
%==========================================================================
function viewCovButton_Callback(hObject, eventdata, handles)

global MRM

nDVs  = size(MRM.Design.Y.Y,2);

%-Find the matrix coordinates
%----------------------------
MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
             str2double(get(handles.YcoordText, 'String'))  ...
             str2double(get(handles.ZcoordText, 'String'))];

[matCoords, ~] = MRM_mni2mat(MNIcoords, [MRM.Options.Out filesep 'MRM_Covar_1_1.nii']);


%-Check mask
%-----------
if isempty(MRM.Options.Mask) ~= 1
    
    maskImg = spm_data_read(spm_data_hdr_read([MRM.Options.Out filesep 'MRM_mask.nii']), ...
                            'Slice', matCoords(3));

    if maskImg(matCoords(1), matCoords(2)) == 0 || isnan(maskImg(matCoords(1), matCoords(2)))
        msgbox('The selected voxel is outside of the masked region');
        return;
    end

    clearvars maskImg
end

% Read in variance-covariance matrix
%-----------------------------------
tempVCOV = zeros(nDVs,nDVs);

for i = 1:nDVs
    for j = 1:nDVs

        slice       = spm_data_read(spm_data_hdr_read([MRM.Options.Out filesep 'MRM_Covar_' ...
                                    num2str(i) '_' num2str(j) '.nii']),                     ...
                                    'Slice', matCoords(3));

        tempVCOV(i,j) = slice(matCoords(1), matCoords(2));
    end
end

clearvars slice

%--------------------------------------------------------------------------
% Whole model
%--------------------------------------------------------------------------
minVal = min(min(tempVCOV));
maxVal = max(max(tempVCOV));

nSubs   = size(MRM.Design.X.X,1);

win = get(0, 'ScreenSize');

    
% Univariate form of the CV matrix
visVCOV = kron(eye(nSubs), tempVCOV);
visVCOV(visVCOV == 0) = NaN;
pos = [win(3)/4, win(4)/4, win(3)/2, win(4)/1.2];

figure('Position',      pos,                             ...
       'Units',         'Pixels',                        ...
       'Color',         [1 1 1],                         ...
       'Name',          'Variance-covariance structure', ...
       'MenuBar',       'None',                          ...
       'DockControls',  'off',                           ...
       'NumberTitle',   'off',                           ...
       'ToolBar',       'none',                          ...
       'Resize',        'on');    

imagesc(visVCOV);
caxis([(minVal - 0.1) maxVal]);
set(gca,'xtick', []);
set(gca,'xticklabel', []);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
cmap = gray;
cmap(1,:) = [0.85 0.85 0.85];
colormap(cmap);

covar = tril(tempVCOV,-1);
covar(covar == 0) = NaN;
minCovar = min(min(covar));
maxCovar = max(max(covar));

if isnan(minCovar) == 1 
    minCovar = 0;
    maxCovar = 0;
end

minVar   = min(diag(tempVCOV));
maxVar   = max(diag(tempVCOV));

string = {'Displayed structure relates to re-expressing $\mathbf{Y}$ as a single vector.', ...
    '',                                                             ...
    ['Covariance = ' num2str(minCovar) ' to ' num2str(maxCovar)],   ...
    ['Variance   = ' num2str(minVar) ' to ' num2str(maxVar)],       ...
    '',                                                             ...
    'Brighter is larger',                                           ...
    ''};

xlabel(string,'interpreter','latex', 'FontSize', 12);

title(['Modelled variance-covariance structure for voxel ' num2str(MNIcoords(1)) ' ' ...
    num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3))],                             ...
    'interpreter','latex',                                                        ...
    'FontSize', 15);
    
%--------------------------------------------------------------------------
% Single subject
%--------------------------------------------------------------------------

visVCOV = tempVCOV;
pos     = [win(3)/4, win(4)/4, win(3)/2, win(4)/1.5];

figure('Position',      pos,                             ...
    'Units',         'Pixels',                        ...
    'Color',         [1 1 1],                         ...
    'Name',          'Variance-covariance structure', ...
    'MenuBar',       'None',                          ...
    'DockControls',  'off',                           ...
    'NumberTitle',   'off',                           ...
    'ToolBar',       'none',                          ...
    'Resize',        'on');

imagesc(visVCOV);
caxis([(minVal - 0.1) maxVal]);
set(gca,'xtick', []);
set(gca,'xticklabel', []);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
cmap = gray;
cmap(1,:) = [0.85 0.85 0.85];
colormap(cmap);

covar = tril(tempVCOV,-1);
covar(covar == 0) = NaN;

labels = cell(1, size(MRM.Design.Y.Cell, 2));

for i = 1:size(MRM.Design.Y.Cell, 2)
    labels{i} = MRM.Design.Y.Cell{i}.Label;
end

for i = 1:nDVs
    for j = 1:nDVs
        if j >= i
            
            text(i,j,sprintf('%.2f',tempVCOV(i,j)), 'Color', 'red', 'FontSize', 12);
        end
    end
end

set(gca,'xtick',      1:size(MRM.Design.Y.Cell, 2));
set(gca,'xticklabel', labels);
set(gca,'FontSize',   12);

rotateXLabels(gca(), 45);

set(gca,'ytick',      1:size(MRM.Design.Y.Cell, 2));
set(gca,'yticklabel', labels);
set(gca,'FontSize',   12);

xlabel(['Modelled variance-covariance structure for voxel ' ...
    num2str(MNIcoords(1)) ' ' num2str(MNIcoords(2)) ' ' ...
    num2str(MNIcoords(3))], ...
    'interpreter','latex', 'FontSize', 15);



  
%==========================================================================
% Crosshairs checkbox
%==========================================================================
function crosshairsCheck_Callback(hObject, eventdata, handles)

global MRM

% Check whether a viewer window is open. If it is refresh
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    crosshairs = get(hObject, 'Value');
    infoPane   = get(handles.infoPaneBox, 'Value');
    hand       = findobj('type','figure','name','MRMviewer');
    pos        = get(hand, 'Position');
    coords     = [str2num(get(handles.XcoordText, 'String')), ...
                  str2num(get(handles.YcoordText, 'String')), ...
                  str2num(get(handles.ZcoordText, 'String'))];
           
    if get(handles.overlayPopup, 'Value') == 2   
        thresh = 0;
    else     
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
        maximaMat  = MRM_getMaxima();
    else
        maximaMat = [];
        thresh    =  str2num(get(handles.thresholdText, 'String'));
    end
    
    handles = guidata(hObject);
    
    MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                           maximaMat, crosshairs, infoPane, handles.template)
    
end


%==========================================================================
% Info pane checkbox
%==========================================================================
function infoPaneBox_Callback(hObject, eventdata, handles)

global MRM

% Check whether a viewer window is open. If it is refresh
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    crosshairs = get(handles.crosshairsCheck, 'Value');
    infoPane   = get(hObject, 'Value');
    hand       = findobj('type','figure','name','MRMviewer');
    pos        = get(hand, 'Position');
    coords     = [str2num(get(handles.XcoordText, 'String')), ...
                  str2num(get(handles.YcoordText, 'String')), ...
                  str2num(get(handles.ZcoordText, 'String'))];
    
    if get(handles.overlayPopup, 'Value') == 2   
        thresh = 0;
    else     
        thresh  =  str2num(get(handles.thresholdText, 'String'));
    end
    
    if strcmp(MRM.Options.Thresh.Level, 'None') ~= 1
        maximaMat  = MRM_getMaxima();
    else
        maximaMat = [];
        thresh    =  str2num(get(handles.thresholdText, 'String'));
    end
    
    handles = guidata(hObject);
    
    MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                           maximaMat, crosshairs, infoPane, handles.template)
    
end



%==========================================================================
% Up and down arrow presses inside the Z coordinate box
%==========================================================================
function ZcoordText_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'uparrow') == 1
   
    coord = str2num(get(handles.ZcoordText, 'String'));
    
    if coord == 103
        return
    end
    
    coord = coord + 1;
    set(handles.ZcoordText, 'String', num2str(coord));
    
    ZcoordText_Callback(handles.ZcoordText, eventdata, handles);
    
elseif strcmp(eventdata.Key, 'downarrow') == 1
    
    coord = str2num(get(handles.ZcoordText, 'String'));
    
    if coord == -71
        return
    end
    
    coord = coord - 1;
    set(handles.ZcoordText, 'String', num2str(coord));
    
    ZcoordText_Callback(handles.ZcoordText, eventdata, handles);
    
end


%==========================================================================
% Up and down arrow presses inside the Y coordinate box
%==========================================================================
function YcoordText_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'uparrow') == 1
   
    coord = str2num(get(handles.YcoordText, 'String'));
    
    if coord == 91
        return
    end
    
    coord = coord + 1;
    set(handles.YcoordText, 'String', num2str(coord));
    
    YcoordText_Callback(handles.YcoordText, eventdata, handles);
    
elseif strcmp(eventdata.Key, 'downarrow') == 1
    
    coord = str2num(get(handles.YcoordText, 'String'));
    
    if coord == -123
        return
    end
    
    coord = coord - 1;
    set(handles.YcoordText, 'String', num2str(coord));
    
    YcoordText_Callback(handles.YcoordText, eventdata, handles);
    
end


%==========================================================================
% Up and down arrow presses inside the X coordinate box
%==========================================================================
function XcoordText_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'uparrow') == 1
   
    coord = str2num(get(handles.XcoordText, 'String'));
    
    if coord == 90
        return
    end
    
    coord = coord + 1;
    set(handles.XcoordText, 'String', num2str(coord));
    
    XcoordText_Callback(handles.XcoordText, eventdata, handles);
    
elseif strcmp(eventdata.Key, 'downarrow') == 1
    
    coord = str2num(get(handles.XcoordText, 'String'));
    
    if coord == -90
        return
    end
    
    coord = coord - 1;
    set(handles.XcoordText, 'String', num2str(coord));
    
    XcoordText_Callback(handles.XcoordText, eventdata, handles);
    
end



%==========================================================================
% Linear Discriminant button
%==========================================================================
function LDAbutton_Callback(hObject, eventdata, handles)

MRM_LDAsetup();







%==========================================================================
% CreateFcn
%==========================================================================

% --- Executes during object creation, after setting all properties.
function CILevel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function DVselectMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function overlayPopup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function assumpTestsMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function contrastExploreMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function whichPlotMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function thresholdText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function conPlotSelectMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ZcoordText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function XcoordText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function YcoordText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
