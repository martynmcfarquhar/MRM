%=========================================================================
% GUI for the specification of the covariate values 
%=========================================================================
% This script creates the GUI used to enter the values for a specified
% covariate
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

function varargout = MRM_covariateValue(varargin)
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_covariateValue_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_covariateValue_OutputFcn, ...
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
function MRM_covariateValue_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

handles.output = hObject;

% Row number
handles.which = varargin{1};

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------%
% Window title
%-------------------------------------------------------------------------%
set(hObject, 'Name', ['Enter values for ' varargin{2}]);

%-------------------------------------------------------------------------%
% Table setup                                                             %
%-------------------------------------------------------------------------%
    
nSubs = size(MRM.Design.X.X,1);

if isfield(MRM.Covariates, 'CV') == 1
    
    if isfield(MRM.Covariates.CV{varargin{1}}, 'Value')
        
        if isfield(MRM.Covariates.CV{varargin{1}}, 'Centering')
        
            if MRM.Covariates.CV{varargin{1}}.Centering == 1
                if isfield(MRM.Covariates.CV{varargin{1}}, 'ValueDM')
                    for i = 1:nSubs
                        data(i,1) = {MRM.Covariates.CV{varargin{1}}.ValueDM(i)};
                    end
                else
                    MRM.Covariates.CV{varargin{1}}.ValueDM = MRM.Covariates.CV{varargin{1}}.Value - mean(MRM.Covariates.CV{varargin{1}}.Value);

                    for i = 1:nSubs
                        data(i,1) = {MRM.Covariates.CV{varargin{1}}.ValueDM(i)};
                    end
                end

            else
                for i = 1:nSubs
                    data(i,1) = {MRM.Covariates.CV{varargin{1}}.Value(i)};
                end
            end
            
            set(handles.demeanCheck, 'Value', MRM.Covariates.CV{varargin{1}}.Centering);
            
        else
            for i = 1:nSubs
                    data(i,1) = {MRM.Covariates.CV{varargin{1}}.Value(i)};
            end
            
            set(handles.demeanCheck, 'Value', 1);
            
        end
            
        set(handles.CVvalueTable, 'Data', data);
        
        
    else 
        data(1:nSubs,1) = {[]};
        set(handles.CVvalueTable, 'Data', data);  
    end    
    
else

    data(1:nSubs,1) = {[]};
    set(handles.CVvalueTable, 'Data', data);
    set(handles.demeanCheck, 'Value', 1);

end

% Data formats for columns
set(handles.CVvalueTable, 'ColumnFormat', {'numeric'});
set(handles.CVvalueTable, 'ColumnName', {varargin{2}});


% --- Outputs from this function are returned to the command line.
function varargout = MRM_covariateValue_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



%==========================================================================
% 'Done' button
%==========================================================================
function doneButton_Callback(hObject, eventdata, handles, varargin)

global MRM

data = cell2mat(get(handles.CVvalueTable, 'Data'));

if sum(isnan(data)) ~= 0 || isempty(data) == 1 || size(data,1) ~= size(MRM.Design.X.X,1)    
    msgbox('There are missing values in the covariate');
    return 
end

MRM.Covariates.CV{handles.which}.Centering = get(handles.demeanCheck, 'Value');
MRM.Covariates.CV{handles.which}.Value     = data;
MRM.Covariates.CV{handles.which}.ValueDM   = data - mean(data);

delete(handles.figure1);



%==========================================================================
% 'Paste from clipboard' Button
%==========================================================================
function pasteClipboardButton_Callback(hObject, eventdata, handles)

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

set(handles.CVvalueTable, 'Data', data);

data = cell2mat(data);

MRM.Covariates.CV{handles.which}.Centering = get(handles.demeanCheck, 'Value');
MRM.Covariates.CV{handles.which}.Value     = data;
MRM.Covariates.CV{handles.which}.ValueDM   = data - mean(data);

% if isfield(MRM.Covariates, 'CV') == 1
%     
%     if isfield(MRM.Covariates.CV{handles.which}, 'Centering')
% 
%         if MRM.Covariates.CV{handles.which}.Centering == 1
%     
%             MRM.Covariates.CV{handles.which}.ValueDM = data - mean(data);
%    
%             if mean(data) > eps*1e5
%                     MRM.Covariates.CV{handles.which}.Value   = data;
%             end
%         else
%              MRM.Covariates.CV{handles.which}.Value = data;
%         end
%     else
%          MRM.Covariates.CV{handles.which}.Value = data;
%     end
%         
% else
%    MRM.Covariates.CV{handles.which}.Value = data;
% end


%==========================================================================
% Window close button has been pressed
%==========================================================================
function figure1_CloseRequestFcn(hObject, eventdata, handles)

answer = questdlg('If you close the window all changes will be lost. Do you want to continue?', ...
                  'Close window?','Yes', 'No', 'No');
              
if strcmp(answer, 'Yes') == 1 
    delete(handles.figure1);
end


%==========================================================================
% Mean centering button
%==========================================================================
function demeanCheck_Callback(hObject, eventdata, handles)

global MRM

MRM.Covariates.CV{handles.which}.Centering = get(hObject, 'Value');

nSubs = size(MRM.Design.X.X,1);

if isfield(MRM.Covariates, 'CV') == 1
    
    if isfield(MRM.Covariates.CV{handles.which}, 'Value')
        
        % Populate with existing details
        
        if MRM.Covariates.CV{handles.which}.Centering == 1
            if isfield(MRM.Covariates.CV{handles.which}, 'ValueDM')
                for i = 1:nSubs
                    data(i,1) = {MRM.Covariates.CV{handles.which}.ValueDM(i)};
                end
            else
                MRM.Covariates.CV{handles.which}.ValueDM = MRM.Covariates.CV{handles.which}.Value - mean(MRM.Covariates.CV{handles.which}.Value);
                for i = 1:nSubs
                    data(i,1) = {MRM.Covariates.CV{handles.which}.ValueDM(i)};
                end
            end
        else
            for i = 1:nSubs
                data(i,1) = {MRM.Covariates.CV{handles.which}.Value(i)};
            end
        end
              
        set(handles.CVvalueTable, 'Data', data);
    end
end
