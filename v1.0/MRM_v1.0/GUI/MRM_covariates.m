%=========================================================================
% GUI for the covariates selector
%=========================================================================
% This script creates the GUI used to display the currently available 
% covariates during the model specification
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

function varargout = MRM_covariates(varargin)
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_covariates_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_covariates_OutputFcn, ...
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
function MRM_covariates_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_covariates
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------%
% Window title
%-------------------------------------------------------------------------%
set(hObject, 'Name', 'Covariate specification');

%-------------------------------------------------------------------------%
% Table setup                                                             %
%-------------------------------------------------------------------------%
    
if isfield(MRM.Covariates, 'CV') == 1

    % Populate with existing details
    data = cell(MRM.Covariates.Number,1);    
    
    for i = 1:MRM.Covariates.Number   
        data(i,1) = {MRM.Covariates.CV{i}.Name};  
    end
    
    set(handles.factorTable, 'Data', data);    
else

    data(1:MRM.Covariates.Number,1) = {''};
    set(handles.factorTable, 'Data', data);
end

% Data formats for columns
set(handles.factorTable, 'ColumnFormat', {'char'});






% --- Outputs from this function are returned to the command line.
function varargout = MRM_covariates_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;






%==========================================================================
% 'Set values...' button
%==========================================================================
function setCVvaluesButton_Callback(hObject, eventdata, handles)

if isfield(handles, 'currentRow') == 1

    row  = handles.currentRow;
    data = get(handles.factorTable, 'Data');
   
    % Check for a name
    if isempty(data{row, 1}) == 1

        msgbox(['Please enter a name for covariate ' num2str(row)]);
        return
    
    end
    
    % Launch the value specification GUI
    MRM_covariateValue(row, data{row,1});
    
else
    
    msgbox('Please select a covariate by clicking in one of the rows.');
    
end




%=========================================================================%
% Done button
%=========================================================================%
function doneButton_Callback(hObject, eventdata, handles)

global MRM

% Check values have been added for each specified covariate
for i = 1:MRM.Covariates.Number
    
    if isfield(MRM.Covariates, 'CV') == 1
        
        if isfield(MRM.Covariates.CV{i}, 'Name')
            
            if isempty(MRM.Covariates.CV{i}.Name) == 1
                
                msgbox('You need to name all the covariates');
                return
                
            elseif isfield(MRM.Covariates.CV{i}, 'Value') == 1
                
                if isempty(MRM.Covariates.CV{i}.Value) == 1
                   
                    msgbox('The value has not been set for at least one covariate');
                    return
                    
                end 
            else
                msgbox('The value has not been set for at least one covariate');
                return
            end
        else
             msgbox('You need to name all the covariates');
             return
        end
    else   
        msgbox('Some information is missing');
        return
    end 
end

MRM_buildDesignMatrix();

delete(handles.figure1);


        
    

%==========================================================================
% Different cell selected in the table
%==========================================================================
function factorTable_CellSelectionCallback(hObject, eventdata, handles)

if numel(eventdata.Indices) > 0

    % Update the handles to contain the last selected row
    handles.currentRow = eventdata.Indices(1);
    guidata(hObject, handles);

end


%==========================================================================
% Name edited in table
%==========================================================================
function factorTable_CellEditCallback(hObject, eventdata, handles)

global MRM

% Name edited
MRM.Covariates.CV{eventdata.Indices(1)}.Name = eventdata.EditData;


%==========================================================================
% 'Cancel' button
%==========================================================================
function cancelButton_Callback(hObject, eventdata, handles)

global MRM

answer = questdlg(['If you cancel all changes will be lost and nothing '  ... 
                   'will be added to your design. Are you sure you want ' ...
                   'to cancel?'], 'Really cancel?', 'Yes', 'No', 'No');
               
if strcmp(answer, 'Yes') == 1 
    
    MRM.Covariates           = [];
    MRM.Covariates.Number    = 0;
    
    delete(handles.figure1);
    
    switch MRM.Model
        case 'Repeated'
            h = findobj('type','figure','name','Repeated measures for neuroimaging');
        case 'MANOVA'
            h = findobj('type','figure','name','MANOVA for neuroimaging');
    end
    handles = guidata(h);
    set(handles.numCVsEdit, 'String', '0');
    set(handles.specifyCVsButton, 'Enable', 'off');
    
end


%==========================================================================
% Close request function
%==========================================================================
function figure1_CloseRequestFcn(hObject, eventdata, handles)

global MRM

everythingIsFine = 1;

% Check values have been added for each specified covariate
for i = 1:MRM.Covariates.Number
    
    if isfield(MRM.Covariates, 'CV') == 1
        
        if isfield(MRM.Covariates.CV{i}, 'Name')
            
            if isempty(MRM.Covariates.CV{i}.Name) == 1
                
                everythingIsFine = 0;
                
            elseif isfield(MRM.Covariates.CV{i}, 'Value') == 1
                
                if isempty(MRM.Covariates.CV{i}.Value) == 1
                   
                    everythingIsFine = 0;
                    
                end 
            else
                everythingIsFine = 0;
            end
        else
             everythingIsFine = 0;
        end
    else   
        everythingIsFine = 0;
    end 
end

if everythingIsFine == 1
    
    MRM_buildDesignMatrix();
    delete(handles.figure1);
    
else
   
    answer = questdlg(['Some information is missing, if you close the '   ...
                       'window now all the changes will be lost. Are '    ...
                       'you sure you want to do this?'], 'Close window?', ...
                       'Yes', 'No', 'No');
                   
    if strcmp(answer, 'Yes') == 1
       
        MRM.Covariates           = [];
        MRM.Covariates.Number    = 0;
        
        delete(handles.figure1);
        
        switch MRM.Model
            case 'Repeated'
                h = findobj('type','figure','name','Repeated measures for neuroimaging');
            case 'MANOVA'
                h = findobj('type','figure','name','MANOVA for neuroimaging');
        end
        handles = guidata(h);
        set(handles.numCVsEdit, 'String', '0');
        set(handles.specifyCVsButton, 'Enable', 'off');
        
        MRM_buildDesignMatrix();
        
    end               
    
end

%==========================================================================
% Delete CV button
%==========================================================================
function deleteCVbutton_Callback(hObject, eventdata, handles)

global MRM

if isfield(handles, 'currentRow') == 1
    
    choice = questdlg('Are you sure you want to delete this covariate?', ...
                      'Delete CV', 'Yes', 'No', 'No');
                  
    if strcmp(choice, 'Yes') == 1

        row  = handles.currentRow;
        MRM.Covariates.CV(row) = [];
        MRM.Covariates.Number = MRM.Covariates.Number - 1;
        
        if MRM.Covariates.Number == 0
            delete(handles.figure1);
            switch MRM.Model
                case 'Repeated'
                    h = findobj('type','figure','name','Repeated measures for neuroimaging');
                case 'MANOVA'
                    h = findobj('type','figure','name','MANOVA for neuroimaging');
            end
            handles = guidata(h);
            set(handles.numCVsEdit, 'String', '0');
            set(handles.specifyCVsButton, 'Enable', 'off');
            MRM_buildDesignMatrix();
        else
            switch MRM.Model
                case 'Repeated'
                    h = findobj('type','figure','name','Repeated measures for neuroimaging');
                case 'MANOVA'
                    h = findobj('type','figure','name','MANOVA for neuroimaging');
            end
            handles = guidata(h);
            set(handles.numCVsEdit, 'String', num2str(MRM.Covariates.Number));
            MRM_covariates();
        end
        
    end
    
else
    
    msgbox('Please select a covariate by clicking in one of the rows.');
    
end
