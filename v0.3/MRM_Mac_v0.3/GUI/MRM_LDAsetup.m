%=========================================================================
% GUI to specify a descriptive linear discriminant analysis
%=========================================================================
% This script creates the GUI used to specify the options for performing an
% dLDA at a specific voxel
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

function varargout = MRM_LDAsetup(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRM_LDAsetup_OpeningFcn, ...
                   'gui_OutputFcn',  @MRM_LDAsetup_OutputFcn, ...
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
function MRM_LDAsetup_OpeningFcn(hObject, eventdata, handles, varargin)

global MRM

% Choose default command line output for MRM_LDAsetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if MRM.Covariates.Number == 0
    set(handles.includeCVsTick, 'Enable', 'off');
end

% Selected voxel
hands = guidata(findobj('type','figure','name','MRM Post-estimation Tools'));

coords = zeros(1,3);

coords(1) = str2num(get(hands.XcoordText, 'String'));
coords(2) = str2num(get(hands.YcoordText, 'String'));
coords(3) = str2num(get(hands.ZcoordText, 'String'));

set(handles.voxCoord, 'String', num2str(coords));

set(hObject, 'Name', 'MRM dLDA');


% --- Outputs from this function are returned to the command line.
function varargout = MRM_LDAsetup_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in designMatGroupButton.
function designMatGroupButton_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.newGroupsButton,    'Value', 0);
    set(handles.defineGroupsButton, 'Enable', 'off');
else
    set(handles.newGroupsButton,    'Value', 1);
    set(handles.defineGroupsButton, 'Enable', 'on');
end


% --- Executes on button press in newGroupsButton.
function newGroupsButton_Callback(hObject, eventdata, handles)

if get(hObject, 'Value') == 1
    set(handles.designMatGroupButton, 'Value', 0);
    set(handles.defineGroupsButton,   'Enable', 'on');
else
    set(handles.designMatGroupButton, 'Value', 1);
    set(handles.defineGroupsButton,   'Enable', 'off');
end


% --- Executes on button press in defineGroupsButton.
function defineGroupsButton_Callback(hObject, eventdata, handles)

MRM_LDAgroups();


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)

figure1_CloseRequestFcn(handles.figure1, eventdata, handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);


%==========================================================================
% Run button
%==========================================================================
function runButton_Callback(hObject, eventdata, handles)

global MRM

%--------------------------------------------------------------------------
% Design matrix definition
%--------------------------------------------------------------------------
if get(handles.designMatGroupButton, 'Value') == 1   
    
    Xlabels = cell(1, size(MRM.Design.X.Cell, 2) + MRM.Covariates.Number);

    for i = 1:size(MRM.Design.X.Cell, 2)
        Xlabels{i} = [MRM.Design.X.Cell{i}.Label ' (' num2str(i) ')'];
    end
else
   if isfield(MRM.Design.X, 'Xlda')
       
       if isempty(MRM.Design.X.Xlda) == 1
           
            msgbox(['You need to define the new groupings if you do not want ' ...
                    'to use the original design matrix.']);
            return 
       else
           
           Xlabels = cell(1,size(MRM.Design.X.Xlda, 2));
           
           for i = 1:size(MRM.Design.X.Xlda, 2)
               Xlabels{i} = num2str(i);
           end
       end
       
   else
       msgbox(['You need to define the new groupings if you do not want ' ...
               'to use the original design matrix.']);
       return
   end
end

%--------------------------------------------------------------------------
% Don't include CVs
%--------------------------------------------------------------------------
if get(handles.includeCVsTick, 'Value') == 0   
    if MRM.Covariates.Number ~= 0   
        if get(handles.designMatGroupButton, 'Value') == 1   
            start    = MRM.Covariates.CV{1}.Col;
            X        = MRM.Design.X.X(:,1:start-1);
            Xgroup   = X;
            GroupNum = size(X,2);
        else  
            X        = MRM.Design.X.Xlda;
            Xgroup   = X;
            GroupNum = size(X,2);
        end 
    else
        if get(handles.designMatGroupButton, 'Value') == 1    
            X        = MRM.Design.X.X;
            Xgroup   = X;
            GroupNum = size(X,2);     
        else   
            X        = MRM.Design.X.Xlda;
            Xgroup   = X;
            GroupNum = size(X,2); 
        end 
    end
%--------------------------------------------------------------------------
% Include CVs
%--------------------------------------------------------------------------    
else
    if MRM.Covariates.Number ~= 0 
        if get(handles.designMatGroupButton, 'Value') == 0  
            Xgroup   = MRM.Design.X.Xlda;
            X        = MRM.Design.X.Xlda;
            GroupNum = size(X,2);
            for i = 1:MRM.Covariates.Number
                if MRM.Covariates.CV{i}.Centering == 1
                    X = [X MRM.Covariates.CV{i}.ValueDM]; 
                else
                    X = [X MRM.Covariates.CV{i}.Value];
                end
            end    
        else
            GroupNum = MRM.Covariates.CV{1}.Col - 1;
            X        = MRM.Design.X.X;
            Xgroup   = X(:,1:GroupNum);
        end    
    else
        if get(handles.designMatGroupButton, 'Value') == 1 
            X        = MRM.Design.X.X;
            Xgroup   = X;
            GroupNum = size(X,2);
        else
            X        = MRM.Design.X.Xlda;
            Xgroup   = X;
            GroupNum = size(X,2); 
        end
    end
end

%--------------------------------------------------------------------------
% Contrast
%--------------------------------------------------------------------------
L = diff(eye(GroupNum));

if MRM.Covariates.Number ~= 0
    if get(handles.includeCVsTick, 'Value') == 1    
        L = [L zeros(size(L,1), MRM.Covariates.Number)];
    end
end

%--------------------------------------------------------------------------
% Get data
%--------------------------------------------------------------------------
MNIcoord = str2num(get(handles.voxCoord, 'String'));

nSubs = size(MRM.Design.Y.Y,1);
nDVs  = size(MRM.Design.Y.Y,2);
cells = size(MRM.Design.X.Cell,2);

files = cell(1,nDVs);

for i = 1:nDVs  
    for j = 1:cells
        files{i} = [files{i} MRM.Data.Y{i}.Cell{j}.Scans];
    end
end

try         
    [matCoords, ~] = MRM_mni2mat(MNIcoord, MRM.Data.Y{1}.Cell{1}.Scans{1});
catch
  msgbox('Cannot find raw data.')
  return
end

Y = zeros(nSubs, nDVs);

for i = 1:nDVs
    for j = 1:nSubs
        Y(j,i) = spm_data_read(spm_data_hdr_read(files{i}{j}), ...
            'xyz',                                             ...
            [matCoords(1) matCoords(2) matCoords(3)]');
    end
end

  
%--------------------------------------------------------------------------
% Run analysis
%--------------------------------------------------------------------------

[Ev, EigVecRaw, UnstandDF, StandDF, varExpl, ~, centroids, functionTests, partialFtests] = MRM_LDA(Y, X, L, GroupNum, ...
      get(handles.returnScoresTick, 'Value'), get(handles.plotScores, 'Value'), ...
      Xgroup, Xlabels, MNIcoord);
  
%--------------------------------------------------------------------------
% Output heading
%--------------------------------------------------------------------------
fprintf('%s\n', ...
    '----------------------------------------------------------------', ...
    'Descriptive Linear Discriminant Analysis'                        , ...
    '----------------------------------------------------------------', ...
    '',                                                                 ...
    ['Results for voxel: ' num2str(MNIcoord)],                          ...
    '');


%--------------------------------------------------------------------------
% Eigenvalues and variance
%--------------------------------------------------------------------------
if get(handles.varExplTick, 'Value') == 1
    fprintf('Eigenvalues and variance explained\n');
    fprintf('----------------------------------\n');
    
    displaytable([Ev varExpl functionTests(:,1) functionTests(:,2) functionTests(:,3) functionTests(:,4)], ...
                 {'Eigenvalues', 'Variance %', 'F', 'df1', 'df2', 'p'}, 11)
    
    fprintf('\n');
    fprintf('\n'); 
end

%--------------------------------------------------------------------------
% Standardised LDF
%--------------------------------------------------------------------------
if get(handles.stdLDTick, 'Value') == 1
   
   fprintf('Standardised discriminant functions\n');
   fprintf('-----------------------------------\n')
    
   labs = cell(1, size(StandDF,2)); 
    
   for i = 1:size(StandDF,2)
       labs{1,i} = ['Function ' num2str(i)];
   end
   
   rows = cell(1, size(MRM.Design.Y.Cell,2));
   
   for i = 1:size(MRM.Design.Y.Cell,2)
       rows{1,i} = MRM.Design.Y.Cell{i}.Label;
   end
   
   displaytable(StandDF, labs, 11, [], rows)
   
   fprintf('\n');
   fprintf('\n');  
end

%--------------------------------------------------------------------------
% Unstandardised LDF
%--------------------------------------------------------------------------
if get(handles.UnstdDFTick, 'Value') == 1
   
   fprintf('Unstandardised discriminant functions\n');
   fprintf('-------------------------------------\n')
    
   labs = cell(1, size(StandDF,2)); 
    
   for i = 1:size(StandDF,2)
       labs{1,i} = ['Function ' num2str(i)];
   end
   
   rows = cell(1, size(MRM.Design.Y.Cell,2));
   
   for i = 1:size(MRM.Design.Y.Cell,2)
       rows{1,i} = MRM.Design.Y.Cell{i}.Label;
   end
   
   rows{1,i+1} = '(Constant)';
   
   displaytable(UnstandDF, labs, 11, [], rows)
   
   fprintf('\n');
   fprintf('\n');   
end

%--------------------------------------------------------------------------
% Group centroids
%--------------------------------------------------------------------------
if get(handles.centroidsTick, 'Value') == 1

   fprintf('Group centroids\n');
   fprintf('---------------\n') 
   
   labs = cell(1, size(StandDF,2)); 
    
   for i = 1:size(StandDF,2)
       labs{1,i} = ['Function ' num2str(i)];
   end
   
    displaytable(centroids, labs, 11, [], Xlabels)

    fprintf('\n');
    fprintf('\n');
end

%--------------------------------------------------------------------------
% Raw eigenvectors
%--------------------------------------------------------------------------
if get(handles.RawEigvecTick, 'Value') == 1
   
   fprintf('Raw eigenvectors\n');
   fprintf('----------------\n')
    
   labs = cell(1, size(StandDF,2)); 
    
   for i = 1:size(StandDF,2)
       labs{1,i} = ['Eigenvector ' num2str(i)];
   end
   
   rows = cell(1, size(MRM.Design.Y.Cell,2));
   
   for i = 1:size(MRM.Design.Y.Cell,2)
       rows{1,i} = MRM.Design.Y.Cell{i}.Label;
   end
   
   displaytable(EigVecRaw, labs, 13, [], rows)
   
   fprintf('\n');
   fprintf('\n');
end

%--------------------------------------------------------------------------
% Partial F-tests
%--------------------------------------------------------------------------

if get(handles.PartialFsTick, 'Value') == 1
    rows = cell(1, size(MRM.Design.Y.Cell,2));

    for i = 1:size(MRM.Design.Y.Cell,2)
        rows{1,i} = MRM.Design.Y.Cell{i}.Label;
    end
    
    fprintf('Partial F-tests for the dependent variables\n');
    fprintf('-------------------------------------------\n');

    displaytable([partialFtests(:,1) partialFtests(:,2) partialFtests(:,3) partialFtests(:,4)], ...
                {'F', 'df1', 'df2', 'p'}, 11, [], rows)
end

fprintf('\n');

msgbox('Results returned in the MATLAB command window');
