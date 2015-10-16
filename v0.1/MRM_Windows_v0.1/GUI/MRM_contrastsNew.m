%=========================================================================
% GUI for the specification of a multivariate contrast
%=========================================================================
% This Script creates the GUI used to specify the contrast weights for a
% general linear hypothesis of the form ABC' = 0
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

function varargout = MRM_contrastsNew(varargin)
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_contrastsNew_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_contrastsNew_OutputFcn, ...
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
function MRM_contrastsNew_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_contrastsNew
handles.output = hObject;

handles.action = varargin{1};

if strcmp(varargin{1}, 'View')
    handles.conNum = varargin{2};
end

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------%
% If the input argument is 'View' then the user clicked the 'View contrast'
% button and we need to fill in the information in the GUI and set all
% cells to editable
%-------------------------------------------------------------------------%
if strcmp(varargin{1}, 'View')
   
    % Use the number of the contrast to fill in the details
    conNum = varargin{2};
    
    % Set the name and existing contrast values
    set(handles.conNameEdit, 'String', MRM.Contrasts.Con{conNum}.Name);
    
    cells = {MRM.Contrasts.Con{conNum}.A};
    set(handles.betweenSubjectsConsTable, 'Data', cells{:});
    
    if MRM.Factors.Between.Number == 0
       numCells = 1;
    else
       numCells = 1;
       for i = 1:MRM.Factors.Between.Number
           numCells = numCells * MRM.Factors.Between.Factors{i}.LevelsNum;
       end 
    end
    
    % Add space for the covariates
    numCells = numCells + MRM.Covariates.Number;
    
    cellWidth = repmat({40}, 1, numCells);
    set(handles.betweenSubjectsConsTable, 'ColumnWidth', cellWidth);
    edit = ones(1, numCells);
    edit = logical(edit);
    set(handles.betweenSubjectsConsTable, 'ColumnEditable', edit);
    
    cells = {MRM.Contrasts.Con{conNum}.C};
    set(handles.withinSubjectsConsTable, 'Data', cells{:});
    
    numCells = 1;

    for i = 1:MRM.Factors.Within.Number  
        numCells = numCells * MRM.Factors.Within.Factors{i}.LevelsNum;  
    end
    
    cellWidth = repmat({40}, 1, numCells);
    set(handles.withinSubjectsConsTable, 'ColumnWidth', cellWidth);
    edit = ones(1,numCells);
    edit = logical(edit);
    set(handles.withinSubjectsConsTable, 'ColumnEditable', edit);
    
end


%-------------------------------------------------------------------------%
% Set the size and editable columns of the contrast tables
%-------------------------------------------------------------------------%
if strcmp(varargin{1}, 'New') == 1

    % Between-subject
    if MRM.Factors.Between.Number == 0
       cells = 1;
    else
       cells = 1;
       for i = 1:MRM.Factors.Between.Number
           cells = cells * MRM.Factors.Between.Factors{i}.LevelsNum;
       end 
    end
    
    % Add space for the covariates
    cells = cells + MRM.Covariates.Number;

    % Create the correct number of table cells
    data = cell(cells,cells);
    set(handles.betweenSubjectsConsTable, 'Data', data);
    
    % Set their width to 40px
    cellWidth = repmat({40}, 1, cells);
    set(handles.betweenSubjectsConsTable, 'ColumnWidth', cellWidth);
    
    % Set all cells to be editable
    edit = ones(1,cells);
    edit = logical(edit);
    set(handles.betweenSubjectsConsTable, 'ColumnEditable', edit);

    % Within-subject
    cells = 1;

    for i = 1:MRM.Factors.Within.Number  
        cells = cells * MRM.Factors.Within.Factors{i}.LevelsNum;  
    end

    % Create the correct number of table cells
    data = cell(cells,cells);
    set(handles.withinSubjectsConsTable, 'Data', data);
    
    % Set their width to 40px
    cellWidth = repmat({40}, 1, cells);
    set(handles.withinSubjectsConsTable, 'ColumnWidth', cellWidth);
    
    % Set all cells to be editable
    edit = ones(1,cells);
    edit = logical(edit);
    set(handles.withinSubjectsConsTable, 'ColumnEditable', edit);
    
    if strcmp(MRM.Model, 'MANOVA') == 1    
        identityMatButtonWS_Callback(handles.identityMatButtonWS, eventdata, handles);
    end
    
end


%-------------------------------------------------------------------------%
% Window name
%-------------------------------------------------------------------------%
set(hObject, 'Name', 'Contrast specifier');




% --- Outputs from this function are returned to the command line.
function varargout = MRM_contrastsNew_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



%==========================================================================
% Add contrast button
%==========================================================================
function addConsButton_Callback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Check a name has been entered
%--------------------------------------------------------------------------
if isempty(get(handles.conNameEdit, 'String')) == 1
    msgbox('You need to specify a name for this contrast.')
    return
end

% Check it is unique
for i = 1:MRM.Contrasts.Number
     if strcmp(MRM.Contrasts.Con{i}.Name, get(handles.conNameEdit, 'String')) == 1
        
         if strcmp(handles.action, 'View')
             
             overwrite = questdlg('Are you sure you want to overwrite the current contrast?', ...
                                  'Overwite?', 'Yes', 'No', 'No');
             
             if strcmp(overwrite, 'No')
                 return
             end
                               
         else
         
            msgbox('Please select a unique name for this contrast.')
            return
         end
        
     end
end

%--------------------------------------------------------------------------
% Check the between-subjects contrasts
%--------------------------------------------------------------------------
BSdata = cell2mat(get(handles.betweenSubjectsConsTable, 'Data'));

% Make sure a value has been entered for each column
if size(BSdata, 2) ~= size(MRM.Design.X.X,2)
    msgbox('Incomplete specification for the ''Between-subjects contrasts''');
    return
elseif sum(isnan(BSdata(:))) ~= 0
    msgbox('NaN''s detected for the ''Between-subjects contrasts''');
    return
end

% Make sure there are no redundant rows - check for all zeros first and
% remove them
BSdata(all(BSdata == 0, 2), :) = [];

if rank(BSdata) ~= size(BSdata,1)
    msgbox([num2str(size(BSdata,1) - rank(BSdata)) ' redundant row(s) detected ' ...
            'in the ''Between-subjects contrasts''']);
    return
end

%--------------------------------------------------------------------------
% Check the within-subjects contrasts
%--------------------------------------------------------------------------
data = cell2mat(get(handles.withinSubjectsConsTable, 'Data'));

% Count the within subjects cells (i.e. columns of Y)
cells = 1;
    
for i = 1:MRM.Factors.Within.Number  
    cells = cells * MRM.Factors.Within.Factors{i}.LevelsNum;  
end

% Make sure a value has been entered for each column
if size(data, 2) ~= cells
    msgbox('Incomplete specification for the ''Within-subjects contrasts''');
    return
elseif sum(isnan(data(:))) ~= 0
    msgbox('NaN''s detected for the ''Within-subjects contrasts''');
    return
end

% Make sure there are no redundant rows - check for all zeros first and
% remove them
data(all(data == 0, 2), :) = [];

if rank(data) ~= size(data,1)
    msgbox([num2str(size(data,1) - rank(data)) ' redundant row(s) detected ' ...
            'in the ''Within-subjects contrasts''']);
    return
end

%--------------------------------------------------------------------------
% If we've made it this far add the name and contrasts to the MRM structure
%--------------------------------------------------------------------------

if exist('overwrite') ~= 0
    
    MRM.Contrasts.Con{handles.conNum}.A  = num2cell(BSdata);
    MRM.Contrasts.Con{handles.conNum}.C  = num2cell(data);
    
else

    % Get the current number of saved contrasts
    num = MRM.Contrasts.Number;

    % Add 1 to it
    newNum = num + 1;

    % Save the new number to the MRM structure
    MRM.Contrasts.Number = newNum;

    % Set the name and values of the new contrast
    MRM.Contrasts.Con{newNum}.Name = get(handles.conNameEdit, 'String');
    MRM.Contrasts.Con{newNum}.A    = num2cell(BSdata);
    MRM.Contrasts.Con{newNum}.C    = num2cell(data);
end

% Reopen the selector
MRM_contrastsSelector();

% Close the window
delete(handles.figure1);



% --- Executes during object creation, after setting all properties.
function conNameEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);

MRM_contrastsSelector();


%--------------------------------------------------------------------------
% Between-subject identity matrix
%--------------------------------------------------------------------------
function identityMatButtonBS_Callback(hObject, eventdata, handles)

global MRM

% Number of cells (i.e. rows/columns of matrix)
if MRM.Factors.Between.Number == 0
   cells = 1;
else
   cells = 1;
   for i = 1:MRM.Factors.Between.Number
       cells = cells * MRM.Factors.Between.Factors{i}.LevelsNum;
   end 
end

% Add space for the covariates
cells = cells + MRM.Covariates.Number;

% Pop an identity matrix into the table
set(handles.betweenSubjectsConsTable, 'Data', num2cell(eye(cells)));


%--------------------------------------------------------------------------
% Within-subject identity matrix
%--------------------------------------------------------------------------
function identityMatButtonWS_Callback(hObject, eventdata, handles)

global MRM

% Number of cells (i.e. rows/columns of matrix)
cells = 1;

for i = 1:MRM.Factors.Within.Number  
    cells = cells * MRM.Factors.Within.Factors{i}.LevelsNum;  
end

idMat = eye(cells);

% Pop an identity matrix into the table
set(handles.withinSubjectsConsTable, 'Data', num2cell(idMat));
