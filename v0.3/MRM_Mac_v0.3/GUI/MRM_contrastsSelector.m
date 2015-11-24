%=========================================================================
% GUI for the specification of the multivariate contrast selector
%=========================================================================
% This script creates the GUI used to display the currently available
% contrasts for a general linear hypothesis test
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

function varargout = MRM_contrastsSelector(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_contrastsSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_contrastsSelector_OutputFcn, ...
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



%=========================================================================%
% Opening function
%=========================================================================%
function MRM_contrastsSelector_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_contrastsSelector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% If there are no contrasts diable the 'View contrast' button and the
% 'Delete contrasts' button
if MRM.Contrasts.Number ~= 0
%    
%     set(handles.viewContrastButton, 'Enable', 'off');
%     set(handles.deleteConButton, 'Enable', 'off');
    
%else
    
    %---------------------------------------------------------------------%
    % Get the contrast names together in a cell and then set the values of
    % the listbox text to that cell so that it displays all the contrast
    % names
    %---------------------------------------------------------------------%
    
    conNames = cell(MRM.Contrasts.Number,1);
    
    for i = 1:MRM.Contrasts.Number
       
        conNames{i} = MRM.Contrasts.Con{i}.Name;
        
    end
    
    set(handles.contrastList, 'String', conNames);
    
end





function varargout = MRM_contrastsSelector_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



%=========================================================================%
% Contrast list
%=========================================================================%
function contrastList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%=========================================================================%
% 'View' button
%=========================================================================%
function viewContrastButton_Callback(hObject, eventdata, handles)

if size(get(handles.contrastList, 'String'),1) ~= 0

    % Get the number of the contrast selected in the list
    conNum = get(handles.contrastList,'Value');

    % Pass the number of the selected contrast to the MRM_contrasts() function
    MRM_contrastsNew('View', conNum);
    
else
    
    msgbox('No contrasts available to view');
    
end




%=========================================================================%
% 'New' button
%=========================================================================%
function newContrastButton_Callback(hObject, eventdata, handles)

MRM_contrastsNew('New');



%=========================================================================%
% 'Delete' button
%=========================================================================%
function deleteConButton_Callback(hObject, eventdata, handles)

global MRM

if size(get(handles.contrastList, 'String'),1) ~= 0

    answer = questdlg('Are you sure you want to delete this contrast?', ... 
                      'Delete contrast?', 'Yes', 'No', 'No');

    if strcmp(answer, 'Yes') == 1
        % Get the number of the contrast selected in the list
        conNum = get(handles.contrastList,'Value');

        % Delete that contrast from the MRM structure
        MRM.Contrasts.Con(conNum) = [];

        % Reduce the number of contrasts in the MRM structure by 1
        MRM.Contrasts.Number = MRM.Contrasts.Number - 1;

        % Refresh the list
        conNames = cell(MRM.Contrasts.Number,1);

        for i = 1:MRM.Contrasts.Number
            conNames{i} = MRM.Contrasts.Con{i}.Name;
        end

        set(handles.contrastList, 'String', conNames);
        set(handles.contrastList, 'Value', 1);

    end
else
    
    msgbox('No contrasts available to delete');
    
end


% --- Executes on button press in generateConsButton.
function generateConsButton_Callback(hObject, eventdata, handles)

global MRM

MRM_generateContrasts();

% Refresh the list
conNames = cell(MRM.Contrasts.Number,1);

for i = 1:MRM.Contrasts.Number
    conNames{i} = MRM.Contrasts.Con{i}.Name;
end

set(handles.contrastList, 'String', conNames);
set(handles.contrastList, 'Value', 1);
