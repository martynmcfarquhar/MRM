%=========================================================================
% GUI for displaying file paths
%=========================================================================
% This script creates the GUI used to display the filepaths for a selected
% DV and cell of the design.
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

function varargout = MRM_viewFiles(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_viewFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_viewFiles_OutputFcn, ...
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


% --- Executes just before MRM_viewFiles is made visible.
function MRM_viewFiles_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for MRM_viewFiles
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(hObject, 'Name', 'File list');

set(handles.fileList, 'String', varargin{1});


% --- Outputs from this function are returned to the command line.
function varargout = MRM_viewFiles_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fileList.
function fileList_Callback(hObject, eventdata, handles)
% 


% --- Executes during object creation, after setting all properties.
function fileList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)

figure1_CloseRequestFcn(handles.figure1, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);
