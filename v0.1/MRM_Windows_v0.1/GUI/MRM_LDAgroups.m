%=========================================================================
% GUI to define new groups for a descriptive linear discriminant analysis
%=========================================================================
% This script creates the GUI used to specify a different grouping
% structure for calculating an dLDA e.g. for main effects
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

function varargout = MRM_LDAgroups(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_LDAgroups_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_LDAgroups_OutputFcn, ...
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
function MRM_LDAgroups_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_LDAgroups
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

nSubs = size(MRM.Design.X.X,1);

if isfield(MRM.Design.X, 'XldaVec') == 1
    data = num2cell(MRM.Design.X.XldaVec);
else
    data(1:nSubs,1) = {[]};
end

set(handles.groupsTable, 'Data', data);

set(hObject, 'Name', 'MRM LDA Groups');


% --- Outputs from this function are returned to the command line.
function varargout = MRM_LDAgroups_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%==========================================================================
% 'Paste' button
%==========================================================================
function pasteButton_Callback(hObject, eventdata, handles)

global MRM

values = importdata('-pastespecial');

% Check the data type
if isa(values, 'double') ~= 1 
    msgbox('The data on the clipboard does not appear to be numeric');
    return   
end

% Check the size of the data
if size(values,1) ~= size(MRM.Design.X.X,1)
    
    msgbox(['You need to provide ' num2str(size(MRM.Design.X.X,1))    ... 
            ' values for this covariate. The data on the clipboard'   ...
            ' contains only ' num2str(size(values,1)) ' values.']);  
    return
end

% If all the above is fine then add to the table
data = cell(size(MRM.Design.X.X,1),1);

for i = 1:size(MRM.Design.X.X,1)
    data(i) = {values(i)}; 
end

set(handles.groupsTable, 'Data', data);


%==========================================================================
% 'View files' button
%==========================================================================
function viewFilesButton_Callback(hObject, eventdata, handles)

global MRM

files = cell(size(MRM.Design.X.X,1),1);

idx = 1;

for i = 1:size(MRM.Data.Y{1}.Cell,2)
    for j = 1:size(MRM.Data.Y{1}.Cell{i}.Scans,2)
        
        files{idx} = MRM.Data.Y{1}.Cell{i}.Scans{j};
        idx = idx + 1;
        
    end
end

MRM_viewFiles(files);


%==========================================================================
% 'Done' button
%==========================================================================
function doneButton_Callback(hObject, eventdata, handles)

global MRM

groupVec = get(handles.groupsTable, 'Data');
groupVec = cell2mat(groupVec);

if size(groupVec,1) ~= size(MRM.Design.X.X,1)
    msgbox('Not enough values entered');
    return
end

newGroups = unique(groupVec);

newX = zeros(size(MRM.Design.X.X,1), size(newGroups,1));

for i = 1:size(newGroups,1)
    vec                  = newX(:,i);
    vec(groupVec == i,:) = 1;
    newX(:,i)            = vec;
end

MRM.Design.X.Xlda    = newX;
MRM.Design.X.XldaVec = groupVec;

delete(handles.figure1);


%==========================================================================
% Executes when user attempts to close
%==========================================================================
function figure1_CloseRequestFcn(hObject, eventdata, handles)

quit = questdlg('Are you sure you want to close this window? All changes will be lost.', ...
                'Quit?', 'Yes', 'No', 'No');

if strcmp(quit, 'Yes') == 1
    delete(handles.figure1)
end


