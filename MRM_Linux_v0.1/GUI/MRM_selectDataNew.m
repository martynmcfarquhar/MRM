%=========================================================================
% GUI for the selection of data files
%=========================================================================
% This script creates the GUI used to select the data on a per-DV per-cell
% basis, adding the filepaths to the correct place in the global MRM
% structure.
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

function varargout = MRM_selectDataNew(varargin)
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_selectDataNew_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_selectDataNew_OutputFcn, ...
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
% 'Select data' button
%==========================================================================
function MRM_selectDataNew_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM
global labels

% Choose default command line output for MRM_selectDataNew
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Window title
set(hObject, 'Name', 'Select files');


% Create the labels and mapping from the labels to the cells of Y and X and
% populate the listbox with the label names
k = 1;

for i = 1:size(MRM.Design.Y.Cell, 2)
    
    if MRM.Factors.Between.Number ~= 0
    
        for j = 1:size(MRM.Design.X.Cell, 2)
       
            labels.name{k} = [MRM.Design.Y.Cell{i}.Label ' ' MRM.Design.X.Cell{j}.Label];
            labels.Y{k} = i; % The cell of Y
            labels.X{k} = j; % The cell of X
        
            k = k + 1;
        end 
    else
       
        labels.name{k} = MRM.Design.Y.Cell{i}.Label;
        labels.Y{k} = i; % The cell of Y
        labels.X{k} = 1; % The cell of X
        
        k = k + 1; 
    end
end

set(handles.cellsList, 'String', labels.name');

cellsList_Callback(handles.cellsList, eventdata, handles);



% --- Outputs from this function are returned to the command line.
function varargout = MRM_selectDataNew_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



% --- Executes on selection change in cellsList.
function cellsList_Callback(hObject, eventdata, handles)

global MRM
global labels

cellNum = get(hObject,'Value');

Yval = labels.Y{cellNum};
Xval = labels.X{cellNum};

% Update number of files text as the selection changes
if Yval <= size(MRM.Data.Y, 2)
    if Xval <= size(MRM.Data.Y{Yval}.Cell, 2)
        
        numFiles = size(MRM.Data.Y{Yval}.Cell{Xval}.Scans,2);
        set(handles.numFilesText, 'String', [num2str(numFiles) ' file(s) selected']);
        set(handles.viewFilesButton, 'Enable', 'on');
        
    else
       
        set(handles.numFilesText, 'String', '0 file(s) selected');
        set(handles.viewFilesButton, 'Enable', 'off');
    end
else
    set(handles.numFilesText, 'String', '0 file(s) selected');
    set(handles.viewFilesButton, 'Enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function cellsList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function selectFilesButton_Callback(hObject, eventdata, handles)

global MRM
global labels

files = uigetfile_n_dir();

% Check that all files were .nii
for i = 1:size(files,2)
    [~, ~, ext] = fileparts(files{i});

    if strcmp(ext, '.nii') == 0
        msgbox('At least one file does not appear to be NIfTI.', ...
               'File error');

        return 
    end
end

% Add to MRM struct
cellNum = get(handles.cellsList,'Value');

Yval = labels.Y{cellNum};
Xval = labels.X{cellNum};

MRM.Data.Y{Yval}.Cell{Xval}.Scans = files;
MRM.Data.Y{Yval}.Cell{Xval}.Label = [MRM.Design.Y.Cell{Yval}.Label ' ' ...
                                     MRM.Design.X.Cell{Xval}.Label];
MRM.Data.Y{Yval}.Cell{Xval}.Levels = MRM.Design.X.Cell{Xval}.Levels;

% Update the text feedback on how many files were selected
set(handles.numFilesText, 'String', [num2str(size(files,2)) ' file(s) selected']);

% If files were selected let the user press the 'View' button
if size(files,2) > 0
    set(handles.viewFilesButton, 'Enable', 'on');
end



% --- Executes on button press in viewFilesButton.
function viewFilesButton_Callback(hObject, eventdata, handles)

global MRM
global labels

cellNum = get(handles.cellsList,'Value');

Yval = labels.Y{cellNum};
Xval = labels.X{cellNum};

MRM_viewFiles(MRM.Data.Y{Yval}.Cell{Xval}.Scans);

% --- Executes on button press in doneButton.
function doneButton_Callback(hObject, eventdata, handles)

% Check the data and files
specOK = MRM_checkSpecification('All', 'both');

if specOK == 1

    MRM_buildDesignMatrix();
    delete(handles.figure1);

elseif specOK == 2
    
    msgbox('There is a mismatch in the number of files for each group across the within-subject variables');
    
else

    msgbox('Some data is missing');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
%
cancelButton_Callback(handles.cancelButton, eventdata, handles);


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)

delete(handles.figure1);
