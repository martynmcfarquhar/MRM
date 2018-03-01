%=========================================================================
% GUI for the global options
%=========================================================================
% This script creates the GUI used to change the global options as saved in
% the MRMoptions.mat file. Defaults are set using the MRM_defaultOptions()
% function
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

function varargout = MRM_optionsGUI(varargin)
%==========================================================================
% Initialisation 
%==========================================================================
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_optionsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_optionsGUI_OutputFcn, ...
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



%==========================================================================
% Opening function
%==========================================================================
function MRM_optionsGUI_OpeningFcn(hObject, eventdata, handles, varargin)

global MRMoptions

% Choose default command line output for MRM_optionsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Find the MRM directory
if exist('MRM_estimate.m', 'file') == 2
        MRMdir = fileparts(which('MRM_estimate'));
        MRMdir = [MRMdir filesep]; 
    else 
        msgbox('Can''t find the MRM directory. Resetting options aborted.');
        return
end   

% If there is no MRMoptions.mat file then make a new one
if exist([MRMdir 'MRMoptions.mat'], 'file') == 0
    MRM_defaultOptions();
end

% Load the options
load([MRMdir 'MRMoptions.mat']);

% Set the Value in the GUI
set(handles.UncorrPthreshDefaultVal, 'String', num2str(MRMoptions.UncorrThresh));
set(handles.FDRqThreshDefaultVal,    'String', num2str(MRMoptions.FDRqThresh));
set(handles.FWEpThreshDefaultVal,    'String', num2str(MRMoptions.FWEpThresh));
%set(handles.bootQbutton,             'Value',  MRMoptions.BootPiZero);
%set(handles.bootNumBox,              'String', num2str(MRMoptions.BootPiResamps));
%set(handles.setPiToOneButton,        'Value',  MRMoptions.PiZeroOne);
set(handles.checkPermsBox,           'Value',  MRMoptions.checkPerms);
%if MRMoptions.pFDR == 0
%    set(handles.typePFDRbutton,      'Value', 0);
%    set(handles.typeFDRbutton,       'Value', 1);
%else
%    set(handles.typePFDRbutton,      'Value', 1);
%    set(handles.typeFDRbutton,       'Value', 0);
%end
set(handles.saveSSCPcheck,           'Value', MRMoptions.SaveSSCP);
set(handles.saveResidsCheck,         'Value', MRMoptions.SaveResid);
set(handles.saveMultivarStatsCheck,  'Value', MRMoptions.SaveMultiStat);
%if get(handles.bootQbutton, 'Value') == 0
%    set(handles.bootNumBox, 'Enable', 'off');
%end
set(handles.defaultPermsBox,         'String', num2str(MRMoptions.nPerms));
set(handles.defaultClusterThreshBox, 'String', num2str(MRMoptions.ClusterP));
set(handles.templateFilepathEdit,    'String', MRMoptions.Template);

if MRMoptions.ClustStat == 1
    set(handles.clusterMassButton,      'Value', 0);
    set(handles.clusterSizeButton,      'Value', 1);
else
    set(handles.clusterMassButton,      'Value', 1);
    set(handles.clusterSizeButton,      'Value', 0);
end

if MRMoptions.Flips == 1
    set(handles.signFlipON,      'Value', 1);
    set(handles.signFlipOFF,     'Value', 0);
else
    set(handles.signFlipON,      'Value', 0);
    set(handles.signFlipOFF,     'Value', 1);
end

switch MRMoptions.Atlas
    case 'TD'
        set(handles.labellingMapMenu, 'Value', 1);
    case 'NM'
        set(handles.labellingMapMenu, 'Value', 2);
    case 'AAL'
        set(handles.labellingMapMenu, 'Value', 3);
end
        


% --- Outputs from this function are returned to the command line.
function varargout = MRM_optionsGUI_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;




%==========================================================================
% 'Bootstrap pi_0' radio button
%==========================================================================
% function bootQbutton_Callback(hObject, eventdata, handles)
% 
% if get(hObject, 'Value') == 1
%     set(handles.setPiToOneButton, 'Value', 0);
%     set(handles.bootNumBox, 'Enable', 'on');
% else
%     set(handles.setPiToOneButton, 'Value', 1);
%     set(handles.bootNumBox, 'Enable', 'off');
% end




%==========================================================================
% 'Set pi_0 = 1' radio button
%==========================================================================
% function setPiToOneButton_Callback(hObject, eventdata, handles)
% 
% if get(hObject, 'Value') == 1
%     set(handles.bootQbutton, 'Value', 0);
%     set(handles.bootNumBox, 'Enable', 'off');
% else
%     set(handles.bootQbutton, 'Value', 1);
%     set(handles.bootNumBox, 'Enable', 'on');
% end




%==========================================================================
% 'Number of resamples' box
%==========================================================================
% function bootNumBox_Callback(hObject, eventdata, handles)
% % 
% 
% % --- Executes during object creation, after setting all properties.
% function bootNumBox_CreateFcn(hObject, eventdata, handles)
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end




%==========================================================================
% Default uncorrected p-value threshold box
%==========================================================================
function UncorrPthreshDefaultVal_Callback(hObject, eventdata, handles)

global MRMoptions

% Check a number has been entered
%--------------------------------
if isnan(str2double(get(hObject, 'String'))) == 1
   msgbox('Please enter a numeric value for the p-value threshold');
   set(handles.UncorrPthreshDefaultVal, 'String', num2str(MRMoptions.UncorrThresh));
   return
end

% Check the value <= 1
%---------------------
if str2double(get(hObject, 'String')) > 1
   msgbox('The p-value threshold cannot be greater than 1')
   set(handles.UncorrPthreshDefaultVal, 'String', num2str(MRMoptions.UncorrThresh));
   return
end 

% --- Executes during object creation, after setting all properties.
function UncorrPthreshDefaultVal_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%==========================================================================
% Default FDR q-value threshold box
%==========================================================================
function FDRqThreshDefaultVal_Callback(hObject, eventdata, handles)

% Check a number has been entered
%--------------------------------
if isnan(str2double(get(hObject, 'String'))) == 1
   msgbox('Please enter a numeric value for the q-value threshold')
   set(handles.FDRqThreshDefaultVal, 'String', num2str(MRMoptions.FDRqThresh));
   return
end

% Check the value <= 1
%---------------------
if str2double(get(hObject, 'String')) > 1
   msgbox('The q-value threshold cannot be greater than 1')
   set(handles.FDRqThreshDefaultVal, 'String', num2str(MRMoptions.FDRqThresh));
   return
end 

% --- Executes during object creation, after setting all properties.
function FDRqThreshDefaultVal_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%==========================================================================
% Default FWE p-value threshold box
%==========================================================================
function FWEpThreshDefaultVal_Callback(hObject, eventdata, handles)

% Check a number has been entered
%--------------------------------
if isnan(str2double(get(hObject, 'String'))) == 1
   msgbox('Please enter a numeric value for the p-value threshold')
   set(handles.FWEpThreshDefaultVal, 'String', num2str(MRMoptions.FWEpThresh));
   return
end

% Check the value <= 1
%---------------------
if str2double(get(hObject, 'String')) > 1
   msgbox('The p-value threshold cannot be greater than 1')
   set(handles.FWEpThreshDefaultVal, 'String', num2str(MRMoptions.FWEpThresh));
   return
end 

% --- Executes during object creation, after setting all properties.
function FWEpThreshDefaultVal_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%==========================================================================
% 'Save' button
%==========================================================================
function saveButton_Callback(hObject, eventdata, handles)

% Pop all the options into the MRMoptions struct and save it back into the
% MRM folder
global MRMoptions

if exist('MRM_estimate.m', 'file') == 2
    MRMdir = fileparts(which('MRM_estimate'));
    MRMdir = [MRMdir filesep];
else
    msgbox('Can''t find the MRM directory. Resetting options aborted.');
    return
end

MRMoptions.UncorrThresh  = str2num(get(handles.UncorrPthreshDefaultVal, 'String'));
MRMoptions.FDRqThresh    = str2num(get(handles.FDRqThreshDefaultVal,    'String'));
MRMoptions.FWEpThresh    = str2num(get(handles.FWEpThreshDefaultVal,    'String'));
MRMoptions.BootPiResamps = 'None';
MRMoptions.BootPiZero    = 0;
MRMoptions.PiZeroOne     = 0;
MRMoptions.pFDR          = get(handles.typePFDRbutton,                  'Value');
MRMoptions.SaveSSCP      = get(handles.saveSSCPcheck,                   'Value');
MRMoptions.SaveResid     = get(handles.saveResidsCheck,                 'Value');
MRMoptions.SaveMultiStat = get(handles.saveMultivarStatsCheck,          'Value');
MRMoptions.nPerms        = str2num(get(handles.defaultPermsBox,         'String'));
MRMoptions.checkPerms    = get(handles.checkPermsBox,                   'Value');
MRMoptions.ClusterP      = str2num(get(handles.defaultClusterThreshBox, 'String'));
MRMoptions.Template      = get(handles.templateFilepathEdit,            'String');
MRMoptions.ClustStat     = get(handles.clusterSizeButton,               'Value');
MRMoptions.Flips         = get(handles.signFlipON,                      'Value');

switch get(handles.labellingMapMenu, 'Value')
    case 1
        MRMoptions.Atlas = 'TD';
    case 2
        MRMoptions.Atlas = 'NM';
    case 3
        MRMoptions.Atlas = 'AAL';
end


save([MRMdir 'MRMoptions.mat'], 'MRMoptions');
delete(handles.figure1);




%==========================================================================
% 'Reset to defaults' button
%==========================================================================
function resetOptionsButton_Callback(hObject, eventdata, handles)

global MRMoptions

% You sure?
answer = questdlg(['All currently set options will be lost. Are you sure you ' ...
                   'wish t continue?'], 'Reset?', 'Yes', 'No', 'No');

if strcmp(answer, 'Yes') == 1
    
    % Find the MRM directory
    if exist('MRM_estimate.m', 'file') == 2
        MRMdir = fileparts(which('MRM_estimate'));
        MRMdir = [MRMdir filesep]; 
    else 
        msgbox('Can''t find the MRM directory. Resetting options aborted.');
        return
    end               

    % Reset the options
    MRM_defaultOptions();

    % Load the new options
    load([MRMdir 'MRMoptions.mat']);
    
    % Set the Value in the GUI
    set(handles.UncorrPthreshDefaultVal, 'String', num2str(MRMoptions.UncorrThresh));
    set(handles.FDRqThreshDefaultVal,    'String', num2str(MRMoptions.FDRqThresh));
    set(handles.FWEpThreshDefaultVal,    'String', num2str(MRMoptions.FWEpThresh));
    %set(handles.bootQbutton,             'Value',  MRMoptions.BootPiZero);
    %set(handles.bootNumBox,              'String', num2str(MRMoptions.BootPiResamps));
    %if get(handles.bootQbutton, 'Value') == 0
    %    set(handles.bootNumBox, 'Enable', 'off');
    %end
    %set(handles.setPiToOneButton,        'Value',  MRMoptions.PiZeroOne);
    
    set(handles.checkPermsBox,           'Value',  MRMoptions.checkPerms);
    
    if MRMoptions.pFDR == 0
        set(handles.typePFDRbutton,      'Value', 0);
        set(handles.typeFDRbutton,       'Value', 1);
    else
        set(handles.typePFDRbutton,      'Value', 1);
        set(handles.typeFDRbutton,       'Value', 0);
    end
    
    set(handles.saveSSCPcheck,           'Value', MRMoptions.SaveSSCP);
    set(handles.saveResidsCheck,         'Value', MRMoptions.SaveResid);
    set(handles.saveMultivarStatsCheck,  'Value', MRMoptions.SaveMultiStat);
    set(handles.defaultPermsBox,         'String', num2str(MRMoptions.nPerms));
    set(handles.defaultClusterThreshBox, 'String', num2str(MRMoptions.ClusterP));
    set(handles.templateFilepathEdit,    'String', MRMoptions.Template);
    
    if MRMoptions.ClustStat == 1
        set(handles.clusterMassButton,      'Value', 0);
        set(handles.clusterSizeButton,      'Value', 1);
    else
        set(handles.clusterMassButton,      'Value', 1);
        set(handles.clusterSizeButton,      'Value', 0);
    end
    
    if MRMoptions.Flips == 1
        set(handles.signFlipON,      'Value', 1);
        set(handles.signFlipOFF,     'Value', 0);
    else
        set(handles.signFlipON,      'Value', 0);
        set(handles.signFlipOFF,     'Value', 1);
    end

    switch MRMoptions.Atlas
        case 'TD'
            set(handles.labellingMapMenu, 'Value', 1);
        case 'NP'
            set(handles.labellingMapMenu, 'Value', 2);
        case 'AAL'
            set(handles.labellingMapMenu, 'Value', 3);
    end
    
end




%==========================================================================
% 'FDR' radio button
%==========================================================================
function typeFDRbutton_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.typePFDRbutton, 'Value', 0);
else
    set(handles.typePFDRbutton, 'Value', 1);
end




%==========================================================================
% 'pFDR' radio button
%==========================================================================
% function typePFDRbutton_Callback(hObject, eventdata, handles)
% 
% if get(hObject, 'Value') == 1
%     set(handles.typeFDRbutton, 'Value', 0);
% else
%     set(handles.typeFDRbutton, 'Value', 1);
% end




%==========================================================================
% 'Save SSCP' checkbox
%==========================================================================
function saveSSCPcheck_Callback(hObject, eventdata, handles)
% 

%==========================================================================
% 'Save residuals' checkbox
%==========================================================================
function saveResidsCheck_Callback(hObject, eventdata, handles)
% 

%==========================================================================
% 'Save multivariate stat' checkbox
%==========================================================================
function saveMultivarStatsCheck_Callback(hObject, eventdata, handles)
%


%==========================================================================
% 'Default no. permutations' box
%==========================================================================
function defaultPermsBox_Callback(hObject, eventdata, handles)
% 

% --- Executes during object creation, after setting all properties.
function defaultPermsBox_CreateFcn(hObject, eventdata, handles)
% 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%==========================================================================
% 'Default cluster threshold' box
%==========================================================================
function defaultClusterThreshBox_Callback(hObject, eventdata, handles)

global MRMoptions

% Check a number has been entered
%--------------------------------
if isnan(str2double(get(hObject, 'String'))) == 1
   msgbox('Please enter a numeric value for the cluster threshold')
   set(handles.hObject, 'String', num2str(MRMoptions.ClusterP));
   return
end

% Check the value <= 1
%---------------------
if str2double(get(hObject, 'String')) > 1
   msgbox('The cluster forming p-value threshold cannot be greater than 1')
   set(handles.hObject, 'String', num2str(MRMoptions.ClusterP));
   return
end 

% --- Executes during object creation, after setting all properties.
function defaultClusterThreshBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function templateFilepathEdit_Callback(hObject, eventdata, handles)

global MRMoptions

filePath = get(hObject, 'String');

[path, ~, ext] = fileparts(filePath);

if isempty(path) == 1
    msgbox('The entered text does not appear to be a file path');
    set(hObject, 'String', MRMoptions.Template);
    return
end

if strcmp(ext, '.gz') == 1
    msgbox('Please unzip the file first');
    set(hObject, 'String', MRMoptions.Template);
    return
elseif strcmp(ext, '.nii') ~= 1
    msgbox('The entered file does not appear to be NifTi');
    set(hObject, 'String', MRMoptions.Template);
    return
end

try
    spm_data_hdr_read(get(hObject, 'String'));
catch error
    msgbox(['There was a problem reading the file: ' error.message]);
    set(hObject, 'String', MRMoptions.Template);
end


% --- Executes during object creation, after setting all properties.
function templateFilepathEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%==========================================================================
% Select template file button
%==========================================================================
function templateSelectFileButton_Callback(hObject, eventdata, handles)

MRMdir = fileparts(which('MRM_estimate'));
MRMdir = [MRMdir filesep]; 

filePath = uigetfile_n_dir([MRMdir 'Utilities' filesep 'Template']);

if isempty(filePath)
    return
end

[path, name, ext] = fileparts(filePath{1});

if strcmp(ext, '.gz') == 1
    msgbox('Please unzip the file first');
    return
elseif strcmp(ext, '.nii') ~= 1
    msgbox('The entered file does not appear to be NifTi');
    return
end

set(handles.templateFilepathEdit, 'String', [path filesep name ext]);


% --- Executes on selection change in labellingMapMenu.
function labellingMapMenu_Callback(hObject, eventdata, handles)
% hObject    handle to labellingMapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns labellingMapMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from labellingMapMenu


% --- Executes during object creation, after setting all properties.
function labellingMapMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labellingMapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkPermsBox.
function checkPermsBox_Callback(hObject, eventdata, handles)
% hObject    handle to checkPermsBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkPermsBox


% --- Executes on button press in clusterSizeButton.
function clusterSizeButton_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.clusterMassButton, 'Value', 0);
else
    set(handles.clusterMassButton, 'Value', 1);
end


% --- Executes on button press in clusterMassButton.
function clusterMassButton_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.clusterSizeButton, 'Value', 0);
else
    set(handles.clusterSizeButton, 'Value', 1);
end


% --- Executes on button press in signFlipON.
function signFlipON_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.signFlipOFF, 'Value', 0);
else
    set(handles.signFlipOFF, 'Value', 1);
end


% --- Executes on button press in signFlipOFF.
function signFlipOFF_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.signFlipON, 'Value', 0);
else
    set(handles.signFlipON, 'Value', 1);
end
