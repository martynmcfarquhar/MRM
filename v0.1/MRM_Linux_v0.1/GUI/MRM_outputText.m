%=========================================================================
% GUI for output text box
%=========================================================================
% This script creates a window where text can be printed -- this isn't used
% anymore as output is now printed to the MATLAB console, however, I've
% kept it in case there's need in the future
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

function varargout = MRM_outputText(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_outputText_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_outputText_OutputFcn, ...
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
function MRM_outputText_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set the display to the input string
set(handles.outputText, 'String', varargin{2});



% --- Outputs from this function are returned to the command line.
function varargout = MRM_outputText_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



% --- Executes during object creation, after setting all properties.
function outputText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%==========================================================================
% 'Done' button
%==========================================================================
function doneButton_Callback(hObject, eventdata, handles)

delete(handles.figure1);


%==========================================================================
% Save button
%==========================================================================
function saveButton_Callback(hObject, eventdata, handles)

[FileName,PathName] = uiputfile('*.txt');

% Check the user didn't cancel
if FileName == 0
    return
else
    FileID = fopen([PathName filesep FileName], 'w'); 
    str = get(handles.outputText, 'String');
    for i = 1:size(str,1)
        if ispc == 1
            fprintf(FileID, [str(i,:) '\r\n']);
        else
            fprintf(FileID, [str(i,:) '\n']);
        end
    end
    fclose(FileID);
end
