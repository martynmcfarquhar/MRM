%=========================================================================
% Load the MRM launcher                 
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% This function defines the launcher GUI and its various functions for
% launching the difference components of MRM.
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

function varargout = MRM_launcher(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_launcher_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_launcher_OutputFcn, ...
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
function MRM_launcher_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for MRM_launcher
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Add folders to path
%--------------------------------------------------------------------------
path = fileparts(which('MRM_launcher'));
addpath(genpath([path filesep 'GUI']));
addpath(genpath([path filesep 'Utilities']));
addpath(genpath([path filesep 'PlotTools']));



% --- Outputs from this function are returned to the command line.
function varargout = MRM_launcher_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



%=========================================================================%
% 'Repeated measures' button 
%=========================================================================%
function modelButtonRM_Callback(hObject, eventdata, handles)

global MRM

if isempty(findobj('type','figure','name','MRM Post-estimation Tools')) ~= 1
    
    answer = questdlg(['Opening the model specification will close the ' ... 
                       'post-estimation tools. Do you want to continue?'], ...
                       'Continue?', 'Yes', 'No', 'No');
                   
    if strcmp(answer, 'Yes')     
        MRM_closeAllWindows();
        MRM = [];
    else
        return
    end
end

if isempty(findobj('type','figure','name','MANOVA for fMRI')) ~= 1
    
    answer = questdlg(['Opening a different model specification will lose the ' ... 
                       'information currently entered. Do you want to continue?'], ...
                       'Continue?', 'Yes', 'No', 'No');
                   
    if strcmp(answer, 'Yes')     
        MRM_closeAllWindows();
        MRM = [];
    else
        return
    end
end


MRM.Model = 'Repeated';

% Launch the model specification window
MRM_main();


%=========================================================================%
% 'MANOVA' button 
%=========================================================================%
function modelButtonMANOVA_Callback(hObject, eventdata, handles)

global MRM

if isempty(findobj('type','figure','name','MRM Post-estimation Tools')) ~= 1
    
    answer = questdlg(['Opening the model specification will close the ' ... 
                       'post-estimation tools. Do you want to continue?'], ...
                       'Continue?', 'Yes', 'No', 'No');
                   
    if strcmp(answer, 'Yes')     
        MRM_closeAllWindows();
        MRM = [];
    else
        return
    end
end

if isempty(findobj('type','figure','name','Repeated measures for fMRI')) ~= 1
    
    answer = questdlg(['Opening a different model specification will lose the ' ... 
                       'information currently entered. Do you want to continue?'], ...
                       'Continue?', 'Yes', 'No', 'No');
                   
    if strcmp(answer, 'Yes')     
        MRM_closeAllWindows();
        MRM = [];
    else
        return
    end
end

MRM.Model = 'MANOVA';

% Launch the model specification window
MRM_main();


%=========================================================================%
% 'Post-estimation tools' button 
%=========================================================================%
function postEstButton_Callback(hObject, eventdata, handles)

global MRM

if isempty(findobj('type','figure','name','Repeated measures for fMRI')) ~= 1 || ...
        isempty(findobj('type','figure','name','MANOVA for fMRI')) ~= 1
    
    answer = questdlg(['Opening the post-estimation tools will lose the information ' ...
        'currently entered in the model specification. Do you want to continue?'], ...
        'Continue?', 'Yes', 'No', 'No');
    
    if strcmp(answer, 'Yes')
        MRM_closeAllWindows();
    else
        return
    end
end


% Get the filename
[FileName,PathName] = uigetfile('.mat', 'Select MRM.mat file');

% Check for cancel
if isequal(FileName,0) 
    return
end

% Try and load the specified file
try
    temp = load([PathName FileName]);
catch error
    msgbox(['There was an error loading the file. ' error.message], ...
            'File error')
    return
end

if isfield(temp, 'MRM')
   
    if temp.MRM.Estimated == 0
        msgbox('This model has not yet been estimated', 'File error')  
        return
    elseif temp.MRM.Contrasts.Number == 0
        msgbox('This model has no contrasts specified', 'File error')  
        return       
    end
    
    MRM = temp.MRM;
    
    clearvars temp
    
else
    msgbox('The loaded file does not appear to be a valid MRM object', ...
           'Error loading MRM')   
    clearvars temp
    return
end

% If results table is open then close it
% if isempty(findobj('type','figure','name','Results Table')) ~= 1
%     
%     h = findobj('type','figure','name','Results Table');
%     delete(h);
%     
% end
% 
% % If viewer window is open then close it
% if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
%     
%     h = findobj('type','figure','name','MRMviewer');
%     delete(h);
%     
%     handles = guidata(hObject);
%     handles.reslicedOverlay = [];
%     guidata(hObject, handles);
%     
% end
% 
% % If plot window is open then close it
% if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1
%     
%     h = findobj('type','figure','name','MRM Bar Plot');
%     delete(h);
%     
% end
% 
% % If post-estimation window is open then close it
% if isempty(findobj('type','figure','name','MRM Post-estimation Tools')) ~= 1
%     
%     h = findobj('type','figure','name','MRM Post-estimation Tools');
%     delete(h);
% end

MRM_closeAllWindows();

% If the output directory does not match the selected file then change it
if strcmp([MRM.Options.Out filesep], PathName) ~= 1  
    MRM.Options.Out = PathName(1,1:size(PathName,2)-1);
    save([PathName 'MRM.mat'], 'MRM');
end


% Launch the post-estimation tools
MRM_postEstimationTools();



%==========================================================================
% 'Help' button
%==========================================================================
function helpButton_Callback(hObject, eventdata, handles)


scriptDir = fileparts(mfilename('fullpath'));

os = computer;

if strfind(os, 'WIN') ~= 0
    winopen([scriptDir filesep 'GUI' filesep 'Help' filesep 'Manual.pdf'])
elseif strfind(os, 'MAC') ~= 0
    system(['open ' scriptDir filesep 'GUI' filesep 'Help' filesep 'Manual.pdf']);
else
    system(['xdg-open ' scriptDir filesep 'GUI' filesep 'Help' filesep 'Manual.pdf']);
end



%==========================================================================
% 'Exit' button
%==========================================================================
function exitButton_Callback(hObject, eventdata, handles)

figure1_CloseRequestFcn(handles.figure1, eventdata, handles);



%==========================================================================
% Close request function
%==========================================================================
function figure1_CloseRequestFcn(hObject, eventdata, handles)

% Say goodbye, it's only polite
fprintf('Bye for now!\n');

% Delete the figure first to prevent warning about recursive close calls
delete(hObject); 
close all;
clear;


%==========================================================================
% 'Options' button
%==========================================================================
function optionsButton_Callback(hObject, eventdata, handles)

% Launch the options window
MRM_optionsGUI();
