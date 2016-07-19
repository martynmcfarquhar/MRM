%=========================================================================
% GUI for the model specification
%=========================================================================
% This script creates the core GUI used to specify a multivariate linear
% model, contrasts, and muliple comparison corrections. The global MRM 
% structure is also defined here and so can be used as reference for other
% functions.
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

function varargout = MRM_main(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_main_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_main_OutputFcn, ...
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
function MRM_main_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for MRM_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Use the current working directory as a save point by default
set(handles.outputDirText, 'String', pwd);


%==========================================================================
% Load options
%==========================================================================
global MRMoptions

% Find the MRM directory
if exist('MRM_estimate.m', 'file') == 2
        MRMdir = fileparts(which('MRM_estimate'));
        MRMdir = [MRMdir filesep]; 
    else 
        msgbox('Can''t find the MRM directory for loading options.');
        return
end   

% If there is no MRMoptions.mat file then make a new one
if exist([MRMdir 'MRMoptions.mat'], 'file') == 0
    MRM_defaultOptions();
end

% Load the options
load([MRMdir 'MRMoptions.mat']);

% Set number of perms
set(handles.permuNumber, 'String', num2str(MRMoptions.nPerms));

% Set default p-value threshold
set(handles.PThreshold, 'String', num2str(MRMoptions.FWEpThresh));

% Set default cluster threshold
set(handles.clusterFormingP, 'String', num2str(MRMoptions.ClusterP));



%==========================================================================%
% Switch model
%==========================================================================%
global MRM

switch MRM.Model  
    case 'Repeated'   
        set(handles.WSFactorsPanel,        'Title',  'Within-subjects factors');
        set(handles.WSSpecifyLevelsButton, 'String', 'Specify levels');
        set(hObject,                       'Name',   'Repeated measures for neuroimaging');
        set(handles.uipanel8,              'Title',  'Factor specification');
        
    case 'MANOVA'    
        set(handles.WSFactorsPanel,        'Title',  'Dependent variables');
        set(handles.WSSpecifyLevelsButton, 'String', 'Set names');
        set(hObject,                       'Name',   'MANOVA for neuroimaging');
        set(handles.uipanel8,              'Title',  'Model specification');       
end

%=========================================================================%
% MRM structure definition and defaults
%=========================================================================%
  MRM.Factors                               = [];     
  MRM.Factors.Within.Number                 = 0;    % Number of within-subject factors
  MRM.Factors.Between.Number                = 0;    % Number of between-subject factors
% MRM.Factors.Within.Factors{i}.Name                - Factor i name
% MRM.Factors.Within.Factors{i}.LevelsNum           - Number levels of Factor i
% MRM.Factors.Within.Factors{i}.Levels{j}           - Factor i level j name

if strcmp(MRM.Model, 'MANOVA') == 1
    MRM.Factors.Within.Number                  = 1;
    MRM.Factors.Within.Factors{1}.Name         = 'DVs';
    MRM.Factors.Within.Factors{1}.LevelsNum    = 1;
    MRM.Factors.Within.Factors{1}.Levels{1}    = '';
end    

  MRM.Covariates                            = [];
  MRM.Covariates.Number                     = 0;
% MRM.Covariates.CV{i}.Name                         % Covariate name
% MRM.Covariates.CV{i}.Value                        % Non mean centered value
% MRM.Covariates.CV{i}.ValueDM                      % Mean centered value
% MRM.Covariates.CV{i}.Centering                    % Centering flag for CV i

  MRM.Contrasts.Number                      = 0;    % Number of contrasts
% MRM.Contrasts.Con{i}.Name                         % Name of contrast i
% MRM.Contrasts.Con{i}.A                            % Between-subject contrast
% MRM.Contrasts.Con{i}.C                            % Within-subject contrast

  MRM.Data                                  = [];   % File paths to the scans
  MRM.Data.Y = [];
% MRM.Data.Y{i}.Label                               - Name of the within-subject cell that
%                                                     Y{i} represents 
%                                                     e.g. 'Time 1 Session 1'
% MRM.Data.Y{i}.Cell{j}.Label                       - First column of Y, cell 1 of X label
% MRM.Data.Y{i}.Cell{j}.Levels                      - The levels of the between-subjects
%                                                     variables that the cell represents 
%                                                     e.g. [1 2] is the first level of
%                                                     Factor 1 and the second level of 
%                                                     Factor 2
% MRM.Data.Y{i}.Cell{j}.Scans                       - Cell strings of file paths to the
%                                                     scans for this cell of the design 
%                                                     for Y{i}

  MRM.Design                                = [];   % Design information                                  
  MRM.Design.X.X                            = [];   % Design matrix
% MRM.Design.X.Cell{i}.Label                        - Label for cell i of the design
% MRM.Design.X.Cell{i}.Levels                       - Levels of the between-subjects
%                                                     factors that code the cell e.g. 
%                                                     [1 1 3] is first level of Factor 1
%                                                     first levels of Factor 2 and third
%                                                     levels of Factor 3

  MRM.Design.Y.Y                            = [];   % Factorial structure of Y
% MRM.Design.Y.Cell{i}.Label                        - Cell i of Y (e.g. the first column
%                                                     of Y) label
% MRM.Design.Y.Cell{i}.Levels                       - Levels of the within-subjects
%                                                     factors that code the cell

MRM.Options.Stat                            = [];                        % Multivariate test statistic
MRM.Options.Stat.Name                       = 'PT';                      % Which test statistic
MRM.Options.Thresh.Level                    = 'Voxel';
MRM.Options.Thresh.Pvals                    = 'Permutation';
MRM.Options.Thresh.Method                   = 'FWE';
MRM.Options.Thresh.nPerms                   = MRMoptions.nPerms;
MRM.Options.Thresh.PThresh                  = MRMoptions.FWEpThresh;     % P-value threshold for correction
MRM.Options.Thresh.ClustThresh              = MRMoptions.ClusterP;
MRM.Options.Out                             = pwd;                       % Output directory
MRM.Options.AllScans                        = 1;                         % Load all scans into memory
MRM.Options.Mask                            = '';                        % Path to a mask file
MRM.Estimated                               = 0;                         % Flag for whether the model has been estimated
        




%=========================================================================%
% Command line outputs
%=========================================================================%

function varargout = MRM_main_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;





%=========================================================================%
% GUI elements (buttons, menus etc.)
%=========================================================================%
% "Number of within-subjects factors" editable text                       %
%-------------------------------------------------------------------------%
function WSFactorsNum_Callback(hObject, eventdata, handles)

global MRM
global labels

if isempty(get(hObject, 'string')) == 1
    set(hObject, 'string', '0')
end

% If 0 is entered clear out the within structure and re-make the design
if str2double(get(hObject, 'string')) == 0
    
    MRM.Factors.Within.Number  = 0;
    MRM.Factors.Within.Factors = [];
    MRM.Design                 = [];
    MRM.Data.Y                 = [];
    MRM.Contrasts              = [];
    MRM.Contrasts.Number       = 0;
    
    if strcmp(MRM.Model, 'MANOVA') == 1
        MRM.Factors.Within.Number                  = 1;
        MRM.Factors.Within.Factors{1}.Name         = 'DVs';
        MRM.Factors.Within.Factors{1}.LevelsNum    = 1;
        MRM.Factors.Within.Factors{1}.Levels{1}    = '';
    end  
    
    MRM_createDesign();
    
    % For visualising the design
    labels.name = [];
    labels.Y    = [];
    labels.X    = [];
    
    set(handles.fileNumText,  'String', '0 file(s) selected');
    set(handles.conNumText,   'String', '0 contrast(s) specified');
    
    set(handles.testStatMenu,    'Enable', 'off');
    set(handles.exactTestButton, 'Enable', 'off');

% If the number of factors has changes clear out the within structure
else
    switch MRM.Model
        
        case 'Repeated'
            
            if str2double(get(hObject, 'string')) ~= MRM.Factors.Within.Number
                MRM.Factors.Within.Number  = str2double(get(hObject, 'string'));
                MRM.Factors.Within.Factors = [];
                MRM.Design                 = [];
                MRM.Data.Y                 = [];
                MRM.Contrasts              = [];
                MRM.Contrasts.Number       = 0;
                
                labels.name = [];
                labels.Y    = [];
                labels.X    = [];
                
                set(handles.fileNumText,     'String', '0 file(s) selected');
                set(handles.conNumText,      'String', '0 contrast(s) specified');
                set(handles.testStatMenu,    'Enable', 'on');
                set(handles.exactTestButton, 'Enable', 'on');
            end
                
        case 'MANOVA'
            
            if str2double(get(hObject, 'string')) ~= MRM.Factors.Within.Factors{1}.LevelsNum
                MRM.Factors.Within.Number                  = 1;
                MRM.Factors.Within.Factors{1}.Name         = 'DVs';
                MRM.Factors.Within.Factors{1}.LevelsNum    = str2double(get(hObject, 'string'));
                MRM.Factors.Within.Factors{1}.Levels       = [];
                MRM.Design                 = [];
                MRM.Data.Y                 = [];
                MRM.Contrasts              = [];
                MRM.Contrasts.Number       = 0;
                
                labels.name = [];
                labels.Y    = [];
                labels.X    = [];
                
                set(handles.fileNumText,     'String', '0 file(s) selected');
                set(handles.conNumText,      'String', '0 contrast(s) specified');
                
            end
            
            if str2double(get(hObject, 'string')) == 1
                set(handles.testStatMenu,    'Enable', 'off');
                set(handles.exactTestButton, 'Enable', 'off');
            else
                set(handles.testStatMenu,    'Enable', 'on');
                set(handles.exactTestButton, 'Enable', 'on');
            end
                
    end  
end
    
    
    



function WSFactorsNum_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






%-------------------------------------------------------------------------%
% "Specify Levels" button: within-subjects                                %
%-------------------------------------------------------------------------%
function WSSpecifyLevelsButton_Callback(hObject, eventdata, handles)

global MRM

switch MRM.Model
    
    case 'Repeated'

        if str2double(get(handles.WSFactorsNum, 'string')) > 0
            MRM_factorsNew(str2double(get(handles.WSFactorsNum,'String')), 'within');
        else
            msgbox(['You need to specify how many factors before you can set ' ...
                    'the levels.']);    
        end
        
    case 'MANOVA'
        
        if str2double(get(handles.WSFactorsNum, 'string')) > 0
            MRM_labelsNew(1, str2double(get(handles.WSFactorsNum, 'string')), 'within');
        else
            msgbox('You need to specify how many DVs before you can name them.');
        end
end




%-------------------------------------------------------------------------%
% "Number of between-subjects factors" editable text                      %
%-------------------------------------------------------------------------%
function BSFactorsNum_Callback(hObject, eventdata, handles)

global MRM
global labels

if isempty(get(hObject, 'string')) == 1
    set(hObject, 'string', '0')
end

% If 0 is entered clear out the between structure and re-make the design
if str2double(get(hObject, 'string')) == 0
    
    MRM.Factors.Between.Number  = 0;
    MRM.Factors.Between.Factors = [];
    MRM.Design                  = [];
    MRM.Data.Y                  = [];
    MRM.Contrasts               = [];
    MRM.Contrasts.Number        = 0;    
    
    MRM_createDesign();
    
    % These are for drawing the design
    labels.name = [];
    labels.Y    = [];
    labels.X    = [];
    
    set(handles.fileNumText, 'String', '0 file(s) selected');
    set(handles.conNumText,  'String', '0 contrast(s) specified');

% If the number of factors has changes clear out the between structure
elseif str2double(get(hObject, 'string')) ~= MRM.Factors.Between.Number
   
    MRM.Factors.Between.Number  = str2double(get(hObject, 'string'));
    MRM.Factors.Between.Factors = [];
    MRM.Design                  = [];
    MRM.Data.Y                  = [];
    MRM.Contrasts               = [];
    MRM.Contrasts.Number        = 0;
    
    labels.name = [];
    labels.Y    = [];
    labels.X    = [];
    
    set(handles.fileNumText, 'String', '0 file(s) selected');
    set(handles.conNumText, 'String', '0 contrast(s) specified');
    
end




function BSFactorsNum_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%-------------------------------------------------------------------------%
% "Specify Levels" button: between-subjects                               %
%-------------------------------------------------------------------------%
function BSSpecifyLevelsButton_Callback(hObject, eventdata, handles)

if str2double(get(handles.BSFactorsNum, 'string')) > 0
    MRM_factorsNew(str2double(get(handles.BSFactorsNum,'String')), 'between');
else
    msgbox(['You need to specify how many factors before you can set ' ...
            'the levels.']);    
end




%-------------------------------------------------------------------------%
% "Select Data" button                                                    %
%-------------------------------------------------------------------------%
function selectDataButton_Callback(hObject, eventdata, handles)

global MRM

% Check all info is there
specOK = MRM_checkSpecification('Factors', 'both');

if specOK == 1
    
    designOK = MRM_createDesign();
    
    if designOK == 1

        h = MRM_selectDataNew();
        uiwait(h);
        
        fileNum = 0;
        
        % Count the files
        for i = 1:size(MRM.Data.Y,2)
            for j = 1:size(MRM.Data.Y{i}.Cell,2)
                
                fileNum = fileNum + size(MRM.Data.Y{i}.Cell{j}.Scans,2);
                
            end
        end
        
        % Set the file number text
        set(handles.fileNumText, 'String', [num2str(fileNum) ' file(s) selected']);
    end 
else

    msgbox('Some factor information is missing. Make sure all fields are complete', ...
            'Incomplete design specification');
end




%-------------------------------------------------------------------------%
% "View design" button                                                    %
%-------------------------------------------------------------------------%
function viewDesignButton_Callback(hObject, eventdata, handles)

MRM_drawDesign(0);




%=========================================================================%
% "Specify contrasts" button                                              %
%=========================================================================%
function specifyConsButton_Callback(hObject, eventdata, handles)

global MRM

designOK = MRM_checkSpecification('Factors', 'both');

if designOK == 1
    
    h = MRM_contrastsSelector();
    uiwait(h);
    
    % Set the contrasts number text
    set(handles.conNumText, 'String', [num2str(MRM.Contrasts.Number) ...
                                      ' contrast(s) specified']);
    
else

    msgbox('Some factor information is missing. Make sure all fields are complete', ...
            'Incomplete design specification');
end





%=========================================================================%
% "Estimate model" button                                                 %
%=========================================================================%
function estimateModelButton_Callback(hObject, eventdata, handles)

global MRM

%-------------------------------------------------------------------------%
% Check the design matrix is not singular
%-------------------------------------------------------------------------%
if isfield(MRM.Design, 'X') == 1
    if isfield(MRM.Design.X, 'X') == 1
        if isempty(MRM.Design.X.X) ~= 1

            X = MRM.Design.X.X;

            if rcond(X'*X) <= eps
                msgbox(['Your design matrix contains redundant columns. This is likely caused ' ...
                        'by covariates you have added. MRM does not use a pseudoinverse to '    ...
                        'solve for the model parameters. As such inv(X''*X) must exist. '       ...
                        'Please check your design and try again.']);
                return
            end
        end
    end
end


% Make directory (if it doesn't exist)
if exist(MRM.Options.Out, 'dir') ~= 7
   
    [success, msg] = mkdir(MRM.Options.Out);
    
    if success == 0
        msgbox(['There was a problem making the output directory: ', msg]);
        return 
    end
end

% Check for existing MRM.mat in the output directory
if exist([MRM.Options.Out filesep 'MRM.mat'], 'file') ~= 0
    
    cont = questdlg(['Existing analysis found in output folder. ' ...
                     'All files will be deleted. Continue?'], ...
                     'Continue?', 'Yes', 'No', 'No');
                 
    if strcmp(cont, 'Yes')
       
       % Delete the folder and contents 
       [success, msg] = rmdir(MRM.Options.Out, 's');
       if success == 0
           msgbox(['There was a problem deleting the contents of the output directory: ', msg]);
           return
       end
       
       % Remake the folder
       [success, msg] = mkdir(MRM.Options.Out);
       if success == 0
           msgbox(['There was a problem making the output directory: ', msg]);
           return
       end
       
       % Check all the info is there
       designOK = MRM_checkSpecification('All', 'both');
       
       if designOK == 1
           
           if MRM.Contrasts.Number == 0
               
              cont = questdlg(['No contrasts have been specified for this model. ' ...
                               'Do you want to continue?'], ...
                               'No contrasts?', 'Yes', 'No', 'No');
                           
              if strcmp(cont, 'Yes')
                  
                  MRM_estimate();
              else
                  return
              end
              
           else
               
               MRM_estimate();
           end
            
       else
           msgbox(['There appears some information missing. Please check your ' ...
                   'design and try again.']);
       end
        
    end 
    
else
    
    % Check all the info is there
    designOK = MRM_checkSpecification('All', 'both');
    
    if designOK == 1
  
        if MRM.Contrasts.Number == 0
            
            cont = questdlg(['No contrasts have been specified for this model. ' ...
                             'Do you want to continue?'], ...
                             'No contrasts?', 'Yes', 'No', 'No');
            
            if strcmp(cont, 'Yes')
                
                MRM_estimate();
            else
                return
            end
            
        else
            
            MRM_estimate();
        end
         
    else
         msgbox(['There appears some information missing. Please check your ' ...
                 'design and try again.']);
    end
    
end




%=========================================================================%
% "Multivariate test statistic" pop-up menu                               %
%=========================================================================%
function testStatMenu_Callback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Change the choice in the MRM structure depending on what has been
% selected
%--------------------------------------------------------------------------
switch get(hObject, 'Value')
    case 1
        MRM.Options.Stat.Name = 'PT';  
    case 2
        MRM.Options.Stat.Name = 'WL';
    case 3
        MRM.Options.Stat.Name = 'HT';
    case 4 
        MRM.Options.Stat.Name = 'RLR';  
end



%=========================================================================%
% "Thresholding" pop-up menu                                              %
%=========================================================================%
function ThresholdLevelMenu_Callback(hObject, eventdata, handles)

% Use the MRM structures
global MRM
global MRMoptions

% Get current selection
selection         = get(hObject,'Value');
thresholdContents = cellstr(get(hObject,'String'));
threshSelection   = thresholdContents{selection};

% Get current p-value method
pValMethodContents  = cellstr(get(handles.PvalCalcMenu,'String'));
pValMethodSelection = pValMethodContents{get(handles.PvalCalcMenu,'Value')};

% Get current correction method
pValCorrectionContents  = cellstr(get(handles.CorrectionMethodMenu,'String'));
pValCorrectionSelection = pValCorrectionContents{get(handles.CorrectionMethodMenu,'Value')};

if strcmp(threshSelection, MRM.Options.Thresh.Level) == 1
    return
else
    switch selection
        
        %------------------------------------------------------------------
        % None
        %------------------------------------------------------------------
        case 1
            
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.PvalCalcMenu,         'String', {'Permutation', 'Approximate'});
            switch pValMethodSelection    
                case 'Permutation'
                    set(handles.permuNumber,  'Enable', 'On');
                case 'Approximate'
                    set(handles.permuNumber,  'Enable', 'Off');
            end
            set(handles.CorrectionMethodMenu, 'Enable', 'Off');
            set(handles.PThreshold,           'Enable', 'Off');
            set(handles.clusterFormingP,      'Enable', 'Off');
         
        %------------------------------------------------------------------    
        % Uncorrected
        %------------------------------------------------------------------
        case 2
            
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.PvalCalcMenu,         'String', {'Permutation', 'Approximate'});
            set(handles.CorrectionMethodMenu, 'Enable', 'Off');
            switch pValMethodSelection    
                case 'Permutation'
                    set(handles.permuNumber,  'Enable', 'On');
                case 'Approximate'
                    set(handles.permuNumber,  'Enable', 'Off');
            end
            set(handles.PThreshold,           'Enable', 'On');
            set(handles.clusterFormingP,      'Enable', 'Off');
            set(handles.correctedPtext,       'String', 'Uncorrected P threshold');
            set(handles.PThreshold,           'String',  num2str(MRMoptions.UncorrThresh));
        
        %------------------------------------------------------------------
        % Voxel
        %------------------------------------------------------------------
        case 3
            
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.PvalCalcMenu,         'String', {'Permutation', 'Approximate'});
            set(handles.CorrectionMethodMenu, 'Enable', 'On');
            set(handles.PThreshold,           'Enable', 'On');
            set(handles.clusterFormingP,      'Enable', 'Off');
            switch pValMethodSelection
                case 'Permutation'
                    set(handles.permuNumber,  'Enable', 'On');
                    set(handles.CorrectionMethodMenu, 'String', {'FWE', 'FDR'});
                    if strcmp(MRM.Options.Thresh.Method, 'FDR') == 1
                        set(handles.CorrectionMethodMenu, 'Value', 2);
                    end
                case 'Approximate'
                    set(handles.permuNumber,  'Enable', 'Off');
                    set(handles.CorrectionMethodMenu, 'Value', 1);
                    set(handles.CorrectionMethodMenu, 'String', {'FDR'});
            end
            switch pValCorrectionSelection
                case 'FWE'
                    set(handles.PThreshold,   'String',  num2str(MRMoptions.FWEpThresh));
                case 'FDR'
                    set(handles.PThreshold,   'String',  num2str(MRMoptions.FDRqThresh));
            end
            set(handles.correctedPtext, 'String', 'Corrected P threshold');
                         
        %------------------------------------------------------------------
        % Cluster
        %------------------------------------------------------------------
        case 4
            
            set(handles.PvalCalcMenu,         'Value', 1);
            set(handles.PvalCalcMenu,         'String', {'Permutation'});
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.CorrectionMethodMenu, 'Value', 1);
            set(handles.CorrectionMethodMenu, 'String', {'FWE'});
            set(handles.permuNumber,          'Enable', 'On');
            set(handles.PThreshold,           'Enable', 'On');
            set(handles.PThreshold,           'String',  num2str(MRMoptions.FWEpThresh));
            set(handles.correctedPtext,       'String', 'Cluster P threshold');
            set(handles.clusterFormingP,      'Enable', 'On');
            set(handles.clusterFormingP,      'String', num2str(MRMoptions.ClusterP));
    end
    
    % Thresholding method
    thresholdContents = cellstr(get(hObject,'String'));
    threshSelection   = thresholdContents{selection};
    
    % Get current p-value method
    pValMethodContents  = cellstr(get(handles.PvalCalcMenu,'String'));
    pValMethodSelection = pValMethodContents{get(handles.PvalCalcMenu,'Value')};

    % Get current correction method
    pValCorrectionContents  = cellstr(get(handles.CorrectionMethodMenu,'String'));
    pValCorrectionSelection = pValCorrectionContents{get(handles.CorrectionMethodMenu,'Value')};
    
    MRM.Options.Thresh.Level                    = threshSelection;
    MRM.Options.Thresh.Pvals                    = pValMethodSelection;
    MRM.Options.Thresh.Method                   = pValCorrectionSelection;
    MRM.Options.Thresh.nPerms                   = str2num(get(handles.permuNumber,     'String'));
    MRM.Options.Thresh.PThresh                  = str2num(get(handles.PThreshold,      'String'));
    MRM.Options.Thresh.ClustThresh              = str2num(get(handles.clusterFormingP, 'String'));
    
end


%=========================================================================%
% "P-value calculation" pop-up menu                                       %
%=========================================================================%
function PvalCalcMenu_Callback(hObject, eventdata, handles)

global MRM
global MRMoptions

% Get current selection
contents = get(hObject,'String');
selection = contents{get(hObject, 'Value')};

% Thresholding method
thresholdContents = cellstr(get(handles.ThresholdLevelMenu,'String'));
threshSelection   = thresholdContents{get(handles.ThresholdLevelMenu,'Value')};

if strcmp(selection, MRM.Options.Thresh.Pvals) == 1
    return
else
    switch selection
        
        case 'Permutation'
            
            set(handles.permuNumber,              'Enable', 'On');
            set(handles.CorrectionMethodMenu,     'String', {'FWE', 'FDR'});
            if strcmp(MRM.Options.Thresh.Method,  'FDR') == 1
                set(handles.CorrectionMethodMenu, 'Value', 2);
            end
            
        case 'Approximate'
            
            set(handles.permuNumber,          'Enable', 'Off');
            set(handles.CorrectionMethodMenu, 'Value',   1);
            set(handles.CorrectionMethodMenu, 'String', {'FDR'});
            
            if strcmp(threshSelection,  'Uncorrected') ~= 1
                set(handles.PThreshold, 'String',  num2str(MRMoptions.FDRqThresh));
            end
    end
end

% Thresholding method
thresholdContents = cellstr(get(handles.ThresholdLevelMenu,'String'));
threshSelection   = thresholdContents{get(handles.ThresholdLevelMenu,'Value')};

% Get current p-value method
pValMethodContents  = cellstr(get(handles.PvalCalcMenu,'String'));
pValMethodSelection = pValMethodContents{get(handles.PvalCalcMenu,'Value')};

% Get current correction method
pValCorrectionContents  = cellstr(get(handles.CorrectionMethodMenu,'String'));
pValCorrectionSelection = pValCorrectionContents{get(handles.CorrectionMethodMenu,'Value')};

MRM.Options.Thresh.Level                    = threshSelection;
MRM.Options.Thresh.Pvals                    = pValMethodSelection;
MRM.Options.Thresh.Method                   = pValCorrectionSelection;
MRM.Options.Thresh.nPerms                   = str2num(get(handles.permuNumber,     'String'));
MRM.Options.Thresh.PThresh                  = str2num(get(handles.PThreshold,      'String'));
MRM.Options.Thresh.ClustThresh              = str2num(get(handles.clusterFormingP, 'String'));



%=========================================================================%
% "Correction method" pop-up menu                                         %
%=========================================================================%
function CorrectionMethodMenu_Callback(hObject, eventdata, handles)

global MRM
global MRMoptions

% Get current selection
contents  = get(hObject,'String');

if ischar(contents) == 1
    selection = contents;
else
    selection = contents{get(hObject, 'Value')};
end

if strcmp(selection, MRM.Options.Thresh.Method) == 1
    return
else
    switch selection
        
        case 'FWE'
            set(handles.PThreshold, 'String',  num2str(MRMoptions.FWEpThresh));
            
        case 'FDR'
            set(handles.PThreshold, 'String',  num2str(MRMoptions.FDRqThresh));
            
    end
end

% Thresholding method
thresholdContents = cellstr(get(handles.ThresholdLevelMenu,'String'));
threshSelection   = thresholdContents{get(handles.ThresholdLevelMenu,'Value')};

% Get current p-value method
pValMethodContents  = cellstr(get(handles.PvalCalcMenu,'String'));
pValMethodSelection = pValMethodContents{get(handles.PvalCalcMenu,'Value')};

% Get current correction method
pValCorrectionContents  = cellstr(get(handles.CorrectionMethodMenu,'String'));
pValCorrectionSelection = pValCorrectionContents{get(handles.CorrectionMethodMenu,'Value')};

MRM.Options.Thresh.Level                    = threshSelection;
MRM.Options.Thresh.Pvals                    = pValMethodSelection;
MRM.Options.Thresh.Method                   = pValCorrectionSelection;
MRM.Options.Thresh.nPerms                   = str2num(get(handles.permuNumber,     'String'));
MRM.Options.Thresh.PThresh                  = str2num(get(handles.PThreshold,      'String'));
MRM.Options.Thresh.ClustThresh              = str2num(get(handles.clusterFormingP, 'String'));


%=========================================================================%
% Number of permutations edit box
%=========================================================================%
function permuNumber_Callback(hObject, eventdata, handles)

global MRM

if isnan(str2double(get(hObject, 'String'))) == 1
    
   set(hObject, 'String', num2str(MRM.Options.Thresh.nPerms));
   return
     
end

MRM.Options.Thresh.nPerms = str2double(get(hObject, 'String'));



%==========================================================================
% P-value threshold edit text
%==========================================================================
function PThreshold_Callback(hObject, eventdata, handles)

global MRM

if isnan(str2double(get(hObject, 'String'))) == 1
    
   set(hObject, 'String', num2str(MRM.Options.Thresh.PThresh));
   return
   
elseif str2double(get(hObject, 'String')) > 1
    
    msgbox('The threshold cannot be greater than 1');
    set(hObject, 'String', num2str(MRM.Options.Thresh.PThresh));
    return
    
end

MRM.Options.Thresh.PThresh = str2double(get(hObject, 'String'));



%==========================================================================
% Cluster-forming p-value threshold
%==========================================================================
function clusterFormingP_Callback(hObject, eventdata, handles)

global MRM

if isnan(str2double(get(hObject, 'String'))) == 1
    
   set(hObject, 'String', num2str(MRM.Options.Thresh.ClustThresh));
   return

elseif str2double(get(hObject, 'String')) > 1
    
    msgbox('The threshold cannot be greater than 1');
    set(hObject, 'String', num2str(MRM.Options.Thresh.ClustThresh));
    return
    
end

MRM.Options.Thresh.ClustThresh = str2double(get(hObject, 'String'));



%=========================================================================%
% "Save" button                                       %
%=========================================================================%
function saveModelButton_Callback(hObject, eventdata, handles)

global MRM

[FileName,PathName] = uiputfile('*.mat');

% Check the user didn't cancel
if FileName ~= 0
    save([PathName FileName], 'MRM');
end




%=========================================================================%
% Output directory button                                                 %
%=========================================================================%
function outDirButton_Callback(hObject, eventdata, handles)

global MRM

% Select the output directory
outDir = uigetdir();

if outDir == 0
    set(handles.outputDirText, 'String', MRM.Options.Out);
else 
    set(handles.outputDirText, 'String', outDir)
    MRM.Options.Out = get(handles.outputDirText, 'String');
end




%=========================================================================%
% Output directory editable text                                          %
%=========================================================================%
function outputDirText_Callback(hObject, eventdata, handles)

global MRM

if isempty(get(hObject, 'String'))
    set(hObject, 'String', pwd);
    MRM.Options.Out = pwd;
else    
    MRM.Options.Out = get(hObject, 'String');
end



%=========================================================================%
% 'Load' button                                                           %
%=========================================================================%
function loadModelButton_Callback(hObject, eventdata, handles)

global MRM

% Get the filename
[FileName,PathName] = uigetfile('.mat');

% Check for cancel
if isequal(FileName,0) 
    return
end

% Try and load the specified file
try
    temp = load([PathName FileName]);
catch
    msgbox('There was an error loading the file', 'File error') 
end

if isfield(temp, 'MRM')
    
    name = get(handles.figure1, 'Name');
    
    if strcmp(name, 'Repeated measures for neuroimaging') == 1
        
        name = 'Repeated';
        
        if strcmp(name, temp.MRM.Model) ~= 1
            msgbox(['This appears to be a MANOVA model. Please open the ' ...
                    'MANOVA module and try loading it again.']);
            return
        end
    else
        name = 'MANOVA';
        
        if strcmp(name, temp.MRM.Model) ~= 1
            msgbox(['This appears to be a repeated measures model. Please ' ...
                   'open the repeated measures module and try loading it again.']);
            return
        end
    end
   
    MRM = temp.MRM;
    clear temp;
    
    %----------------------------------------------------------------------
    % Fill in the GUI
    %----------------------------------------------------------------------
    
    % Within-subjects number
    %-----------------------
    switch MRM.Model
        
        case 'Repeated'
    
            set(handles.WSFactorsNum, 'String', num2str(MRM.Factors.Within.Number));

            if MRM.Factors.Within.Number == 0
                set(handles.testStatMenu,    'Enable', 'off');
                set(handles.exactTestButton, 'Enable', 'off');
            else
                set(handles.testStatMenu,    'Enable', 'on');
                set(handles.exactTestButton, 'Enable', 'on');
            end
            
        case 'MANOVA'
            
            set(handles.WSFactorsNum, 'String', num2str(MRM.Factors.Within.Factors{1}.LevelsNum));

            if MRM.Factors.Within.Factors{1}.LevelsNum == 0
                set(handles.testStatMenu,    'Enable', 'off');
                set(handles.exactTestButton, 'Enable', 'off');
            else
                set(handles.testStatMenu,    'Enable', 'on');
                set(handles.exactTestButton, 'Enable', 'on');
            end
    end
    
    % Between-subjects number
    %------------------------
    set(handles.BSFactorsNum, 'String', num2str(MRM.Factors.Between.Number));
    
    % File number
    %------------
    if isempty(MRM.Data.Y) ~= 1  
        
        fileNum = 0;
        
        for i = 1:size(MRM.Data.Y,2)
            for j = 1:size(MRM.Data.Y{i}.Cell,2)
                
                fileNum = fileNum + size(MRM.Data.Y{i}.Cell{j}.Scans,2);
                
            end
        end
        
        set(handles.fileNumText, 'String', [num2str(fileNum) ' file(s) selected']);
        
    end
    
    % Covariates number
    %------------------
    set(handles.numCVsEdit, 'String', num2str(MRM.Covariates.Number));
    
    if MRM.Covariates.Number ~= 0
        set(handles.specifyCVsButton, 'Enable', 'on');
    else
        set(handles.specifyCVsButton, 'Enable', 'off');
    end
    
    % Contrast number
    %----------------
    set(handles.conNumText, 'String', [num2str(MRM.Contrasts.Number) ' contrast(s) specified']);

    
    % Multivariate stat
    %------------------
    switch MRM.Options.Stat.Name
        
        case 'PT'
            set(handles.testStatMenu, 'Value', 1);
            
        case 'WL'
            set(handles.testStatMenu, 'Value', 2);
            
        case 'HT'
            set(handles.testStatMenu, 'Value', 3);
            
        case 'RLR'
            set(handles.testStatMenu, 'Value', 4);
    end
    
    
    % P-value method menu
    %--------------------
    switch MRM.Options.Thresh.Level
        
        case 'None'
            set(handles.ThresholdLevelMenu,   'Value',   1);
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.correctedPtext,       'String', 'Corrected P threshold');
            set(handles.PvalCalcMenu,         'String', {'Permutation', 'Approximate'});
            set(handles.CorrectionMethodMenu, 'Enable', 'Off');
            set(handles.permuNumber,          'Enable', 'Off');
            set(handles.PThreshold,           'Enable', 'Off');
            set(handles.clusterFormingP,      'Enable', 'Off');
            set(handles.CorrectionMethodMenu,  'String', {'FWE', 'FDR'});
            switch MRM.Options.Thresh.Pvals
                case 'Permutation'
                    set(handles.PvalCalcMenu, 'Value', 1);
                    set(handles.permuNumber,  'Enable', 'On');
                case 'Approximate'
                    set(handles.PvalCalcMenu, 'Value', 2);
                    set(handles.permuNumber,  'Enable', 'Off');
            end
            
        case 'Uncorrected'   
            set(handles.ThresholdLevelMenu,   'Value',   2);
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.correctedPtext,       'String', 'Uncorrected P threshold');
            set(handles.PvalCalcMenu,         'String', {'Permutation', 'Approximate'});
            set(handles.CorrectionMethodMenu, 'String', {'FWE', 'FDR'});
            set(handles.CorrectionMethodMenu, 'Enable', 'Off');
            set(handles.PThreshold,           'Enable', 'On');
            set(handles.clusterFormingP,      'Enable', 'Off');
            set(handles.correctedPtext,       'String', 'Uncorrected P threshold');
            switch MRM.Options.Thresh.Pvals    
                case 'Permutation'
                    set(handles.permuNumber,  'Enable', 'On');
                    set(handles.PvalCalcMenu, 'Value', 1);
                case 'Approximate'
                    set(handles.permuNumber,  'Enable', 'Off');
                    set(handles.PvalCalcMenu, 'Value', 2);
            end
            
        case 'Voxel'
            set(handles.ThresholdLevelMenu,   'Value',   3);
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.correctedPtext,       'String', 'Corrected P threshold');
            set(handles.PvalCalcMenu,         'String', {'Permutation', 'Approximate'});
            set(handles.CorrectionMethodMenu, 'Enable', 'On');
            set(handles.PThreshold,           'Enable', 'On');
            set(handles.clusterFormingP,      'Enable', 'Off');
            set(handles.CorrectionMethodMenu, 'String', {'FWE', 'FDR'});
            
            switch MRM.Options.Thresh.Pvals
                case 'Permutation'
                    set(handles.PvalCalcMenu, 'Value', 1);
                    set(handles.permuNumber,  'Enable', 'On');
                    
                    switch MRM.Options.Thresh.Method
                        case 'FWE'
                            set(handles.CorrectionMethodMenu, 'Value', 1);
                        case 'FDR'
                            set(handles.CorrectionMethodMenu, 'Value', 2);
                    end
                    
                case 'Approximate'
                    set(handles.PvalCalcMenu, 'Value', 2);
                    set(handles.permuNumber,  'Enable', 'Off');
                    set(handles.CorrectionMethodMenu, 'Value', 1);
                    set(handles.CorrectionMethodMenu, 'String', 'FDR');
            end
            
  
        case 'Cluster'
            set(handles.ThresholdLevelMenu,   'Value',   4);
            set(handles.correctedPtext,       'String', 'Cluster P threshold');
            set(handles.PvalCalcMenu,         'Enable', 'On');
            set(handles.PvalCalcMenu,         'Value',   1);
            set(handles.PvalCalcMenu,         'String', 'Permutation');
            set(handles.CorrectionMethodMenu, 'Enable', 'On');
            set(handles.CorrectionMethodMenu, 'Value',   1);
            set(handles.CorrectionMethodMenu, 'String', 'FWE');
            set(handles.permuNumber,          'Enable', 'On');
            set(handles.PThreshold,           'Enable', 'On');
            set(handles.clusterFormingP,      'Enable', 'On');
            set(handles.CorrectionMethodMenu,  'String', 'FWE');
            
    end
    
    % P-threshold
    %------------
    set(handles.PThreshold,      'String', num2str(MRM.Options.Thresh.PThresh));
    set(handles.permuNumber,     'String', num2str(MRM.Options.Thresh.nPerms));
    set(handles.clusterFormingP, 'String', num2str(MRM.Options.Thresh.ClustThresh));
    
    % Output directory
    %-----------------
    set(handles.outputDirText, 'String', MRM.Options.Out);
    
    % Mask file
    %----------
    set(handles.maskDirText, 'String', MRM.Options.Mask);
    
else
    
    msgbox('The loaded file does not appear to be a valid MRM object', ...
           'Error loading MRM')
       
    clear temp;
    
end




%=========================================================================%
% "Help" button                                                           %
%=========================================================================%
function helpButton_Callback(hObject, eventdata, handles)

% Find the scripts directory
scriptDir = fileparts(mfilename('fullpath'));

% Find out if we are in Windows. If we are use the 'winopen' command,
% otherwise use a 'system' call to open the PDF

if isempty(strfind(computer, 'WIN')) ~= 1
    winopen([scriptDir filesep 'Help' filesep 'Manual.pdf'])
elseif isempty(strfind(computer, 'MAC')) ~= 1
    system(['open ' scriptDir filesep 'Help' filesep 'Manual.pdf']);
else
    system(['xdg-open ' scriptDir filesep 'Help' filesep 'Manual.pdf']);
end




%=========================================================================%
% "Exit" button                                                           %
%=========================================================================%

function exitButton_Callback(hObject, eventdata, handles)

figure1_CloseRequestFcn(handles.figure1, eventdata, handles);




%=========================================================================%
% Mask text field
%=========================================================================%
function maskDirText_Callback(hObject, eventdata, handles)

global MRM

MRM.Options.Mask = get(hObject,'String');



%=========================================================================%
% Mask directory selection button
%=========================================================================%
function maskDirButton_Callback(hObject, eventdata, handles)

global MRM

[filename, pathname, filterindex] = uigetfile('*.nii', 'Pick a mask file');

if filename ~= 0
    MRM.Options.Mask = [pathname, filename];
    set(handles.maskDirText, 'String', MRM.Options.Mask);
end




%==========================================================================
% Close request function
%==========================================================================
function figure1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);




%==========================================================================
% 'Number of covariates' edit text box
%==========================================================================
function numCVsEdit_Callback(hObject, eventdata, handles)

global MRM

% Check a number has been entered
if isnan(str2double(get(hObject, 'String'))) == 1
    
    set(hObject, 'String', num2str(MRM.Covariates.Number));
    
elseif str2double(get(hObject, 'String')) > 0
    
    set(handles.specifyCVsButton, 'Enable', 'on');
    
    % If the number of covariates has changed then empty the CV structure
    % and contrasts and start again
    
    if str2double(get(hObject, 'String')) ~= MRM.Covariates.Number
        
        MRM.Covariates           = [];
        MRM.Covariates.Number    = str2double(get(hObject, 'String'));
         
        % Remove existing contrasts
        MRM.Contrasts         = [];
        MRM.Contrasts.Number  = 0; 
        set(handles.conNumText, 'String', '0 contrast(s) specified');
        
        % Check all the other info has been specified before updating the
        % model
        designOK = MRM_checkSpecification('All', 'both');
        
        if designOK == 1
            MRM_createDesign();
            MRM_buildDesignMatrix();
        end
        
    end
    
else
    
    set(handles.specifyCVsButton, 'Enable', 'off');
    MRM.Covariates           = [];
    MRM.Covariates.Number    = str2double(get(hObject, 'String'));
    
    % Remove existing contrasts
    MRM.Contrasts         = [];
    MRM.Contrasts.Number  = 0; 
    set(handles.conNumText, 'String', '0 contrast(s) specified');
    
    designOK = MRM_checkSpecification('All', 'both');
        
    if designOK == 1
        MRM_createDesign();
        MRM_buildDesignMatrix();
    end
end


%==========================================================================
% 'Specify values' covariates button
%==========================================================================
function specifyCVsButton_Callback(hObject, eventdata, handles)

global MRM

% Check the model has been fully specified first
if isfield(MRM.Design, 'X') == 1
    if isfield(MRM.Design.X, 'X') == 1
        if isempty(MRM.Design.X.X) == 1  
            
            msgbox('You need to select all the data before you can add covariates to the design.');
            set(handles.specifyCVsButton, 'Enable', 'off');
            set(handles.numCVsEdit, 'String', '0');
            MRM.Covariates           = [];
            MRM.Covariates.Number    = 0;
            return 
        end
    else
        
        msgbox('You need to select all the data before you can add covariates to the design.');
        set(handles.specifyCVsButton, 'Enable', 'off');
        set(handles.numCVsEdit, 'String', '0');
        MRM.Covariates           = [];
        MRM.Covariates.Number    = 0;
        return
    end
else
    
    msgbox('You need to select all the data before you can add covariates to the design.');
    set(handles.specifyCVsButton, 'Enable', 'off');
    set(handles.numCVsEdit, 'String', '0');
    MRM.Covariates           = [];
    MRM.Covariates.Number    = 0;
    return
end

MRM_covariates();


%==========================================================================
% Design Orthogonality button
%==========================================================================
function designOrthogButton_Callback(hObject, eventdata, handles)

MRM_designOrthogonality();



%==========================================================================
% Exact test button
%==========================================================================
function exactTestButton_Callback(hObject, eventdata, handles)

global MRM

if MRM.Contrasts.Number == 0
    
    msgbox('Please specify your contrasts before performing the exact test.')
    
else
    
    largestRank = 1;
    X           = MRM.Design.X.X;
    fakeData    = rand(size(MRM.Design.X.X,2), size(MRM.Design.Y.Y,2));

    for i = 1:MRM.Contrasts.Number
        
        A = cell2mat(MRM.Contrasts.Con{i}.A);
        C = cell2mat(MRM.Contrasts.Con{i}.C);

        r = rank((A * fakeData * C')' * inv(A * inv(X'*X) * A') * (A * fakeData * C'));

        if r > largestRank
            largestRank = r;
        end

    end

    if largestRank == 1
        msgbox(['All contrasts are exact. As such all the multivariate statistics '   ...
                'will be identical. If you are using permutation tests then Wilks'' ' ...
                'lambda is recommended as the fastest statistic to compute.']);
    else
        msgbox(['At least one of your contrasts will not be exact and as such the ' ...
                'multivariate statistics will differ.']);
    end

end



%==========================================================================
% CreateFcn
%==========================================================================
function CorrectionMethodMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PvalCalcMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ThresholdLevelMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numCVsEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PThreshold_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clusterFormingP_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function permuNumber_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maskDirText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function outputDirText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function testStatMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
