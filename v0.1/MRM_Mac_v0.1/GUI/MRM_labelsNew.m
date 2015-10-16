%=========================================================================
% GUI for the specification of factor level names
%=========================================================================
% This script creates the GUI used to specify the labels for the levels of
% a specified factor in the design
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

function varargout = MRM_labelsNew(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_labelsNew_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_labelsNew_OutputFcn, ...
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


%==========================================================================
% Opening function
%==========================================================================
function MRM_labelsNew_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_labelsNew
handles.output = hObject;

% Factor number, levels, and type 
handles.facNum   = varargin{1};
handles.levelNum = varargin{2};
handles.Which    = varargin{3};

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Setup table
%--------------------------------------------------------------------------

% See if the factor has already been specified
if strcmp(handles.Which, 'within') == 1
    
   if isempty(MRM.Factors.Within.Factors{handles.facNum}.Levels) == 1

       data(1:handles.levelNum, 1) = {''};
        
   else
       data(1:handles.levelNum, 1) = {''};
       
       for i = 1:handles.levelNum    
           data{i,1} = MRM.Factors.Within.Factors{handles.facNum}.Levels{i};    
       end
       
   end
   
   set(handles.labelsTable, 'Data', data);
   
   switch MRM.Model
       case 'Repeated'
           set(handles.labelsTable, 'ColumnName', {'Labels'});
           set(handles.labelsText,  'String',     ['Labels for levels of ' varargin{4}]);
           set(hObject,             'Name',       'Level labels');
       case 'MANOVA'
           set(handles.labelsTable, 'ColumnName', {'Names'});
           set(handles.labelsText,  'String',     'DV names');
           set(hObject,             'Name',       'DVs');
   end

   
elseif strcmp(handles.Which, 'between') == 1
    
    if isempty(MRM.Factors.Between.Factors{handles.facNum}.Levels) == 1

        % Number of rows
        data(1:handles.levelNum, 1) = {''};
        set(handles.labelsTable, 'Data', data);
        set(handles.labelsTable, 'ColumnName', {'Labels'});
        
   else
       data(1:handles.levelNum, 1) = {''};
       
       for i = 1:handles.levelNum   
           data{i,1} = MRM.Factors.Between.Factors{handles.facNum}.Levels{i}; 
       end
       
       set(handles.labelsTable, 'Data', data);
       set(handles.labelsTable, 'ColumnName', {'Labels'});
    end
    
    set(handles.labelsText, 'String', ['Labels for levels of ' varargin{4}]);
    set(hObject,            'Name',   'Level labels');
    
end





% --- Outputs from this function are returned to the command line.
function varargout = MRM_labelsNew_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in doneButton.
function doneButton_Callback(hObject, eventdata, handles)

global MRM

data = get(handles.labelsTable, 'Data');

%--------------------------------------------------------------------------
% Check a label has been entered for each cell
%--------------------------------------------------------------------------
for i = 1:handles.levelNum
    
    if strcmp(data{i,1}, '') == 1
        
        if strcmp(handles.Which, 'within') == 1
        
            switch MRM.Model
                case 'Repeated'
                    msgbox('You need to enter a label for each level');
                    return
                case 'MANOVA'
                    msgbox('You need to enter a name for each DV');
                    return
            end
            
        else
            msgbox('You need to enter a label for each level');
           return
        end
    end
end

%--------------------------------------------------------------------------
% Add labels to MRM structure
%--------------------------------------------------------------------------
if strcmp(handles.Which, 'within') == 1
    
    for i = 1:handles.levelNum
        MRM.Factors.Within.Factors{handles.facNum}.Levels{i} =  regexprep(data{i,1},'[^\w'']','');
    end
    
    if strcmp(MRM.Model, 'MANOVA') == 1

        % Make the design if this is a MANOVA (because we missed out MRM_factorsNew())
        MRM_createDesign();

    end
    
elseif strcmp(handles.Which, 'between') == 1
    
    for i = 1:handles.levelNum
        MRM.Factors.Between.Factors{handles.facNum}.Levels{i} = regexprep(data{i,1},'[^\w'']','');
    end
    
end


% Close window
delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% 


% --- Executes when entered data in editable cell(s) in labelsTable.
function labelsTable_CellEditCallback(hObject, eventdata, handles)
% 
    
