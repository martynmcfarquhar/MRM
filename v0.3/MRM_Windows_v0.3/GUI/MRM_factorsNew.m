%=========================================================================
% GUI for the specification of factor names and levels
%=========================================================================
% This script creates the GUI used to specify the names and levels of
% factors in the design
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

function varargout = MRM_factorsNew(varargin)
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_factorsNew_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_factorsNew_OutputFcn, ...
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
% Initisation code for execution before the window is made visible
%=========================================================================%
function MRM_factorsNew_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_factorsNew
handles.output = hObject;

% Which factors are we specifying
handles.Which = varargin{2};

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------%
% Window title
%-------------------------------------------------------------------------%
if strcmp(handles.Which, 'between') == 1
    
    set(hObject, 'Name', 'Between-subjects factor specification');
    
elseif strcmp(handles.Which, 'within') == 1
    
    set(hObject, 'Name', 'Within-subjects factor specification'); 
    
end


%-------------------------------------------------------------------------%
% Table setup                                                             %
%-------------------------------------------------------------------------%

if strcmp(handles.Which, 'within') == 1
    
    if isempty(MRM.Factors.Within.Factors) == 1

        % Number of rows
        data(1:varargin{1},1) = {''};
        data(1:varargin{1},2) = {[]};
        set(handles.factorTable, 'Data', data);
        
    else
        
        data(1:varargin{1},1) = {''};
        data(1:varargin{1},2) = {[]};
        
        for i = 1:size(MRM.Factors.Within.Factors, 2)
            
            data{i,1} = MRM.Factors.Within.Factors{i}.Name;
            data{i,2} = MRM.Factors.Within.Factors{i}.LevelsNum;
            
        end
        
        set(handles.factorTable, 'Data', data);
        
    end
    
elseif strcmp(handles.Which, 'between') == 1
    
    if isempty(MRM.Factors.Between.Factors) == 1

        % Number of rows
        data(1:varargin{1},1) = {''};
        data(1:varargin{1},2) = {[]};
        set(handles.factorTable, 'Data', data);
        
    else
        
        data(1:varargin{1},1) = {''};
        data(1:varargin{1},2) = {[]};
        
        for i = 1:size(MRM.Factors.Between.Factors, 2)
            
            data{i,1} = MRM.Factors.Between.Factors{i}.Name;
            data{i,2} = MRM.Factors.Between.Factors{i}.LevelsNum;
            
        end
        
        set(handles.factorTable, 'Data', data);
        
    end
end

% Data formats for columns
set(handles.factorTable, 'ColumnFormat', {'char', 'numeric'});






% --- Outputs from this function are returned to the command line.
function varargout = MRM_factorsNew_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;






%==========================================================================
% 'Set labels...' button
%==========================================================================
function setLabelsButton_Callback(hObject, eventdata, handles)

global MRM

if isfield(handles, 'currentRow') == 1

    row  = handles.currentRow;
    data = get(handles.factorTable, 'Data');

    %----------------------------------------------------------------------
    % Check for a name and levels specified
    %----------------------------------------------------------------------
    if isempty(data{row, 1}) == 1

        msgbox(['Please enter a name for factor ' num2str(row)]);
        return

    elseif isempty(data{row, 2}) == 1

        msgbox(['Please specify the number of levels for factor ' num2str(row)]);
        return
    
    end
    
    %----------------------------------------------------------------------
    % Check that levels != 1
    %----------------------------------------------------------------------
    if data{row, 2} < 2
       
       msgbox('You need more than 1 level for a factor');
       return
    
    end
    
    %----------------------------------------------------------------------
    % Add some initial info to the MRM structure
    %----------------------------------------------------------------------
    if strcmp(handles.Which, 'within') == 1
        
        if isfield(MRM.Factors.Within, 'Factors') == 1 
            
            if iscell(MRM.Factors.Within.Factors) == 1
                
                if size(MRM.Factors.Within.Factors,2) < row  
                    
                    MRM.Factors.Within.Factors{row}.Name      = data{row,1};
                    MRM.Factors.Within.Factors{row}.Levels    = [];
                    MRM.Factors.Within.Factors{row}.LevelsNum = data{row,2};
                    
                elseif isfield(MRM.Factors.Within.Factors{row}, 'Levels') == 1
                    
                    if data{row,2} ~= MRM.Factors.Within.Factors{row}.LevelsNum

                        MRM.Factors.Within.Factors{row}.Name      = data{row,1};
                        MRM.Factors.Within.Factors{row}.Levels    = [];
                        MRM.Factors.Within.Factors{row}.LevelsNum = data{row,2}; 
                    end
                else 
                    MRM.Factors.Within.Factors{row}.Name      = data{row,1};
                    MRM.Factors.Within.Factors{row}.Levels    = [];
                    MRM.Factors.Within.Factors{row}.LevelsNum = data{row,2}; 
                end
            else 
                MRM.Factors.Within.Factors{row}.Name      = data{row,1};
                MRM.Factors.Within.Factors{row}.Levels    = [];
                MRM.Factors.Within.Factors{row}.LevelsNum = data{row,2}; 
            end
        else
            MRM.Factors.Within.Factors{row}.Name      = data{row,1};
            MRM.Factors.Within.Factors{row}.Levels    = [];
            MRM.Factors.Within.Factors{row}.LevelsNum = data{row,2}; 
        end
        
    elseif strcmp(handles.Which, 'between') == 1
        
        if isfield(MRM.Factors.Between, 'Factors') == 1 
            
            if iscell(MRM.Factors.Between.Factors) == 1
            
                if size(MRM.Factors.Between.Factors,2) < row  
                    
                    MRM.Factors.Between.Factors{row}.Name      = data{row,1};
                    MRM.Factors.Between.Factors{row}.Levels    = [];
                    MRM.Factors.Between.Factors{row}.LevelsNum = data{row,2};
                
                elseif isfield(MRM.Factors.Between.Factors{row}, 'Levels') == 1

                    if data{row,2} ~= MRM.Factors.Between.Factors{row}.LevelsNum

                        MRM.Factors.Between.Factors{row}.Name      = data{row,1};
                        MRM.Factors.Between.Factors{row}.Levels    = [];
                        MRM.Factors.Between.Factors{row}.LevelsNum = data{row,2}; 
                    end
                else 
                    MRM.Factors.Between.Factors{row}.Name      = data{row,1};
                    MRM.Factors.Between.Factors{row}.Levels    = [];
                    MRM.Factors.Between.Factors{row}.LevelsNum = data{row,2}; 
                end
            else
                MRM.Factors.Between.Factors{row}.Name      = data{row,1};
                MRM.Factors.Between.Factors{row}.Levels    = [];
                MRM.Factors.Between.Factors{row}.LevelsNum = data{row,2};                 
            end
        else
            MRM.Factors.Between.Factors{row}.Name      = data{row,1};
            MRM.Factors.Between.Factors{row}.Levels    = [];
            MRM.Factors.Between.Factors{row}.LevelsNum = data{row,2}; 
        end

    end  
    
    %----------------------------------------------------------------------
    % Launch the labels GUI
    %----------------------------------------------------------------------
    MRM_labelsNew(row, data{row,2}, handles.Which, data{row,1});
    
else
    
    msgbox('Please select a row by clicking in one of the cells');
    
end




%=========================================================================%
% Done button
%=========================================================================%
function doneButton_Callback(hObject, eventdata, handles)

% Loop through the table column setting the names, just in case some one
% has changed the name of the factor but not set any new levels.

if strcmp(handles.Which, 'within') == 1
    specOK = MRM_checkSpecification('Factors', 'within');
else
    specOK = MRM_checkSpecification('Factors', 'between');
end

    
if specOK == 1

    MRM_createDesign();
    delete(gcf);
    
else
    
    msgbox('Some factor information is missing');
    
end
        
    


% --- Executes when selected cell(s) is changed in factorTable.
function factorTable_CellSelectionCallback(hObject, eventdata, handles)

if numel(eventdata.Indices) > 0

    % Update the handles to contain the last selected row
    handles.currentRow = eventdata.Indices(1);
    guidata(hObject, handles);

end


%==========================================================================
% Cells edited in the factor table
%==========================================================================
function factorTable_CellEditCallback(hObject, eventdata, handles)

global MRM
global labels

row  = handles.currentRow;

%--------------------------------------------------------------------------
% Name edited
%--------------------------------------------------------------------------
if eventdata.Indices(2) == 1
    
    if strcmp(handles.Which, 'within') == 1
        
        MRM.Factors.Within.Factors{row}.Name = eventdata.EditData;
        
    elseif strcmp(handles.Which, 'between') == 1
        
        MRM.Factors.Between.Factors{row}.Name = eventdata.EditData;
    
    end  

%--------------------------------------------------------------------------
% Levels edited  
%--------------------------------------------------------------------------
elseif eventdata.Indices(2) == 2
    
    %-Within
    if strcmp(handles.Which, 'within') == 1
        if isfield(MRM.Factors.Within, 'Factors') == 1
            if isfield(MRM.Factors.Within.Factors, 'Levels') == 1
                if str2double(eventdata.EditData) ~= size(MRM.Factors.Within.Factors{row}.Levels, 2)
                    
                    MRM.Factors.Within.Factors{row}.Levels    = [];
                    MRM.Factors.Within.Factors{row}.LevelsNum = str2double(eventdata.EditData);
                    MRM.Design                 = [];
                    MRM.Data.Y                 = [];
                    MRM.Contrasts              = [];
                    MRM.Contrasts.Number       = 0;
                    
                    labels.name = [];
                    labels.Y    = [];
                    labels.X    = [];
                    
                    hand = findobj('type','figure','name','Repeated measures for neuroimaging');
                    
                    if isempty(hand)
                        hand = findobj('type','figure','name','MANOVA for neuroimaging');
                    end
                
                    h =  guidata(hand);
                    
                    set(h.fileNumText,     'String', '0 file(s) selected');
                    set(h.conNumText,      'String', '0 contrast(s) specified');
                    set(h.testStatMenu,    'Enable', 'on');
                    set(h.exactTestButton, 'Enable', 'on');
                    
                end
            end
        end

          if isempty(eventdata.NewData) ~= 1
              
              MRM.Factors.Within.Factors{row}.Levels    = [];
              MRM.Factors.Within.Factors{row}.LevelsNum = str2double(eventdata.EditData);
              MRM.Design                 = [];
              MRM.Data.Y                 = [];
              MRM.Contrasts              = [];
              MRM.Contrasts.Number       = 0;
              
              labels.name = [];
              labels.Y    = [];
              labels.X    = [];
              
              hand = findobj('type','figure','name','Repeated measures for neuroimaging');
              
              if isempty(hand)
                  hand = findobj('type','figure','name','MANOVA for neuroimaging');
              end
              
              h =  guidata(hand);
              
              set(h.fileNumText,     'String', '0 file(s) selected');
              set(h.conNumText,      'String', '0 contrast(s) specified');
              set(h.testStatMenu,    'Enable', 'on');
              set(h.exactTestButton, 'Enable', 'on');
              
          end
              
     %-Between
     elseif strcmp(handles.Which, 'between') == 1
         if isfield(MRM.Factors.Between, 'Factors') == 1 
            if isfield(MRM.Factors.Between.Factors, 'Levels') == 1
                if str2double(eventdata.EditData) ~= size(MRM.Factors.Between.Factors{row}.Levels, 2)
                    
                    MRM.Factors.Between.Factors{row}.Levels    = [];
                    MRM.Factors.Between.Factors{row}.LevelsNum = str2double(eventdata.EditData);
                    MRM.Design                 = [];
                    MRM.Data.Y                 = [];
                    MRM.Contrasts              = [];
                    MRM.Contrasts.Number       = 0;
                    
                    labels.name = [];
                    labels.Y    = [];
                    labels.X    = [];
                    
                    hand = findobj('type','figure','name','Repeated measures for neuroimaging');
                    
                    if isempty(hand)
                        hand = findobj('type','figure','name','MANOVA for neuroimaging');
                    end
                    
                    h =  guidata(hand);
                    
                    set(h.fileNumText,     'String', '0 file(s) selected');
                    set(h.conNumText,      'String', '0 contrast(s) specified');
%                     set(h.testStatMenu,    'Enable', 'off');
%                     set(h.exactTestButton, 'Enable', 'off');
              
                    
                end
            end
        end

          if isempty(eventdata.NewData) ~= 1
              
              MRM.Factors.Between.Factors{row}.Levels    = [];
              MRM.Factors.Between.Factors{row}.LevelsNum = str2double(eventdata.EditData);
              MRM.Design                 = [];
              MRM.Data.Y                 = [];
              MRM.Contrasts              = [];
              MRM.Contrasts.Number       = 0;
              
              labels.name = [];
              labels.Y    = [];
              labels.X    = [];
              
              hand = findobj('type','figure','name','Repeated measures for neuroimaging');
              
              if isempty(hand)
                  hand = findobj('type','figure','name','MANOVA for neuroimaging');
              end
              
              h =  guidata(hand);
              
              set(h.fileNumText,     'String', '0 file(s) selected');
              set(h.conNumText,      'String', '0 contrast(s) specified');
%              set(h.testStatMenu,    'Enable', 'off');
%              set(h.exactTestButton, 'Enable', 'off');
              
          
          end
    end
    
end
