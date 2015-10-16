%=========================================================================
% GUI for the display of an interactive results table
%=========================================================================
% This script creates the GUI used to display the results read-in from a
% results text file
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

function varargout = MRM_displayResultsTable(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_displayResultsTable_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_displayResultsTable_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


%==========================================================================
% Opening function
%==========================================================================
function MRM_displayResultsTable_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM
global MRMoptions

MRM_getOptions();

% Choose default command line output for MRM_displayResultsTable
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% If window position is provided
if size(varargin,2) == 2
    set(hObject, 'Unit', 'pixels');
    set(hObject, 'Position', varargin{2});   
end

% Selected contrast
con = varargin{1};

name = regexprep(MRM.Contrasts.Con{con}.Name, ' ', '_');

if strcmp(MRM.Options.Thresh.Pvals, 'Permutation') == 1
    filepath = [MRM.Options.Out filesep 'PermutationContrasts' filesep name filesep ...
                'MRM_' name '_Results.txt'];
elseif strcmp(MRM.Options.Thresh.Pvals, 'Approximate') == 1
    filepath = [MRM.Options.Out filesep 'Contrasts' filesep name filesep            ...
                'MRM_' name '_Results.txt'];
end

tablePath = fileparts(filepath);

if exist(tablePath,'dir') == 0 
   msgbox('There were no results for the selected contrast');
   delete(handles.resultsTableWindow);
   return
end

if strcmp(MRM.Options.Thresh.Level, 'None') == 1
   msgbox('No results table is available when no thresholding has been used');
   delete(handles.resultsTableWindow);
   return
end

colWidth = 0;

if strcmp(MRM.Options.Thresh.Level, 'Voxel') == 1 || strcmp(MRM.Options.Thresh.Level, 'Uncorrected') == 1

    % Create table data
    fileID = fopen(filepath);
    input = textscan(fileID,'%s %d %.3f %.3f %.3f %d %d %d', ...
                     'Delimiter','|', ...
                     'HeaderLines',1);

    inputNew = cell(size(input{1},1), 3);
    
    % See whether the cluster column is whole numbers. If it is then we
    % have cluster extent, if not it's cluster mass
    temp = input{3}(:);
    tot  = sum(mod(temp,1) == 0);
    
    names = get(handles.resultsTable, 'ColumnName');

    if tot ~= size(temp,1)
        names{3} = 'Mass';
    else
        names{3} = 'Extent';
    end

    for i = 1:size(input{1},1)

        % Atlas label
        inputNew{i,1} = input{1}{i};
        
        if (size(inputNew{i,1},2) * 6) > colWidth % Guess of 6 pixels per character
            colWidth = size(inputNew{i,1},2) * 6;
        end

        % Cluster number
        inputNew{i,2} = sprintf('%d',input{2}(i));

        % Extent/Mass
        if strcmp(names{3}, 'Mass') == 1
            inputNew{i,3} = sprintf('%.3f',input{3}(i));
        else
            inputNew{i,3} = sprintf('%d',input{3}(i));
        end

        % F-value
        inputNew{i,4} = sprintf('%.3f',input{4}(i));

        % P-value
        if strcmp(sprintf('%.3f',input{5}(i)), '0.000') == 1;  
            inputNew{i,5} = '< 0.001';
        else
            inputNew{i,5} = ['   ' sprintf('%.3f',input{5}(i))];
        end

        % Coordinates
        for j = 6:8   
            inputNew{i,j} = input{j}(i);   
        end

    end

    columnformat = {'char', 'char', 'char', 'char',   ...
                    'char', 'char', 'char', 'char'};
    columnWidth  = {colWidth, 'auto', 'auto', 'auto', ...
                    'auto', 'auto', 'auto', 'auto'}; 

elseif strcmp(MRM.Options.Thresh.Level, 'Cluster') == 1
    
    % Create table data
    fileID = fopen(filepath);
    input = textscan(fileID,'%s %d %.3f %.3f %d %d %d', ...
                     'Delimiter','|',                 ...
                     'HeaderLines',1);

    inputNew = cell(size(input{1},1), 3);
    
    % See whether the cluster column is whole numbers. If it is then we
    % have cluster extent, if not it's cluster mass
    temp = input{3}(:);
    tot = sum(mod(temp,1) == 0);
    
    names = get(handles.resultsTable, 'ColumnName');

    if tot ~= size(temp,1)
        names{3} = 'Mass';
    else
        names{3} = 'Extent';
    end

    for i = 1:size(input{1},1)

        % Atlas label
        inputNew{i,1} = input{1}{i};
        
        if (size(inputNew{i,1},2) * 4) > colWidth
            colWidth = size(inputNew{i,1},2) * 4;
        end

        % Cluster number
        inputNew{i,2} = sprintf('%d',input{2}(i));

        % Extent/Mass
        if strcmp(names{3}, 'Mass') == 1
            inputNew{i,3} = sprintf('%.3f',input{3}(i));
        else
            inputNew{i,3} = sprintf('%d',input{3}(i));
        end

        % P-value
        if strcmp(sprintf('%.3f',input{4}(i)), '0.000') == 1;  
            inputNew{i,4} = '< 0.001';
        else
            inputNew{i,4} = ['   ' sprintf('%.3f',input{4}(i))];
        end

        % Coordinates
        for j = 5:7   
            inputNew{i,j} = input{j}(i);   
        end

    end

    columnformat = {'char', 'char', 'char', 'char', ...
                    'char', 'char', 'char'};
    columnWidth  = {colWidth, 'auto', 'auto', 'auto', ...
                    'auto', 'auto', 'auto'}; 
end

set(handles.resultsTable, 'ColumnWidth',  columnWidth);            
set(handles.resultsTable, 'ColumnFormat', columnformat);            
set(handles.resultsTable, 'Data',         inputNew);

% Change the p-value labelling

switch MRM.Options.Thresh.Level
        case 'None'
              names{5} = 'Uncorrected p-value (voxel)';
        case 'Uncorrected'
              names{5} = 'Uncorrected p-value (voxel)';
        case 'Voxel' 
            switch MRM.Options.Thresh.Method
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

set(handles.resultsTable, 'ColumnName', names);

fclose(fileID);
 
set(handles.resultsTable, 'Units', 'pixels');
set(hObject, 'Units', 'pixels');

tablePos   = get(handles.resultsTable, 'Position');
winPos     = get(hObject, 'Position');
winPos(3)  = tablePos(3) + 47;
set(hObject, 'Position', winPos);

set(handles.resultsTable, 'Units', 'normalized');
set(hObject, 'Units', 'normalized');




% --- Outputs from this function are returned to the command line.
function varargout = MRM_displayResultsTable_OutputFcn(hObject, eventdata, handles) 

try
    varargout{1} = handles.output;
end


%==========================================================================
% Cell selected
%==========================================================================
function resultsTable_CellSelectionCallback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Get X, Y, and Z coordinates
%--------------------------------------------------------------------------
row       = eventdata.Indices(1);
tableData = get(handles.resultsTable, 'data');

if strcmp(MRM.Options.Thresh.Level, 'Cluster') == 1
    X = tableData{row,5};
    Y = tableData{row,6};
    Z = tableData{row,7};
else
    X = tableData{row,6};
    Y = tableData{row,7};
    Z = tableData{row,8};
end    

h          = findobj('type','figure','name','MRM Post-estimation Tools');
handles    = guidata(h);
crosshairs = get(handles.crosshairsCheck, 'Value');
infoPn     = get(handles.infoPaneBox, 'Value');

set(handles.XcoordText, 'String', X);
set(handles.YcoordText, 'String', Y); 
set(handles.ZcoordText, 'String', Z);

%--------------------------------------------------------------------------
% If a viewer window is open then update it with the new coordinates
%--------------------------------------------------------------------------
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    
    maximaMat = MRM_getMaxima();
    
    
    hand    = findobj('type','figure','name','MRMviewer');
    pos     = get(hand, 'Position');
    
    coords  = [str2num(get(handles.XcoordText, 'String')), ...
               str2num(get(handles.YcoordText, 'String')), ...
               str2num(get(handles.ZcoordText, 'String'))];
    
    if get(handles.overlayPopup, 'Value') == 2   
        thresh = 0;
    else     
        thresh = str2num(get(handles.thresholdText, 'String'));
    end
    
    MRM_overlayImageNewNew(handles.reslicedOverlay, coords, thresh, pos, ...
                           maximaMat, crosshairs, infoPn, handles.template)
    
end

%--------------------------------------------------------------------------
% If a plot window is open then update it with the new coordinates
%--------------------------------------------------------------------------
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1

    selection = get(handles.conPlotSelectMenu,'Value');
    
    out   = MRM.Options.Out;
    nDVs  = size(MRM.Design.Y.Y,2);

    %----------------------------------------------------------------------
    % Find the matrix coordinates
    %----------------------------------------------------------------------
    MNIcoords = [str2double(get(handles.XcoordText, 'String')) ...
                 str2double(get(handles.YcoordText, 'String'))  ...
                 str2double(get(handles.ZcoordText, 'String'))];

    [matCoords, ~] = MRM_mni2mat(MNIcoords, [out filesep 'MRM_PE_1_1.nii']);
    
    %----------------------------------------------------------------------
    % Read in PEs
    %----------------------------------------------------------------------
    tempPE = zeros(size(MRM.Design.X.X,2),nDVs);
    
    for i = 1:size(MRM.Design.X.X,2)
        for j = 1:nDVs
            
            file = [out filesep 'MRM_PE_' num2str(i) '_' num2str(j) '.nii'];
            
            tempPE(i,j) = spm_data_read(spm_data_hdr_read(file), ...
                                        'xyz', ...
                                        [matCoords(1) matCoords(2) matCoords(3)]');
        end
    end
    
    %----------------------------------------------------------------------
    % Read in variance-covariance matrix
    %----------------------------------------------------------------------
    tempVCOV = zeros(nDVs,nDVs);
    
    for i = 1:nDVs
        for j = 1:nDVs
            
            file = [out filesep 'MRM_Covar_' num2str(i) '_' num2str(j) '.nii'];
            
            tempVCOV(i,j) = spm_data_read(spm_data_hdr_read(file), ...
                                          'xyz', ...
                                          [matCoords(1) matCoords(2) matCoords(3)]');
        end
    end
    
    %----------------------------------------------------------------------
    % Get contrasts
    %----------------------------------------------------------------------
    A = MRM.Contrasts.Plots.Con{selection - 1}.A;
    C = MRM.Contrasts.Plots.Con{selection - 1}.C;

    %----------------------------------------------------------------------
    % Draw plot
    %----------------------------------------------------------------------
    MRM_barPlot(A, tempPE, C, tempVCOV, get(handles.plotCIsButton, 'Value'), ...
                MNIcoords, MRM.Contrasts.Plots.Con{selection - 1}.Name,      ...
                0, [], str2double(get(handles.CILevel, 'String')) / 100);

    clearvars tempVCOV tempPE 
end



%==========================================================================
% Close function
%==========================================================================
function resultsTableWindow_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);
