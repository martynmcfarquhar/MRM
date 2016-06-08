%==========================================================================
% MRM Contrast Specifier
%==========================================================================
% Function description
%--------------------------------------------------------------------------
% Define the GUI and functions for the contrast specification window. An
% MRM struct as well as an action flag need to be passed. The flag 'action'
% can either be 'plot' or 'con' depending on whether the contrasts are for
% plotting results or for estimation. This only influences where they are
% placed in the MRM struct. If 'con' they go in MRM.Contrasts.Con{i} and if
% they are 'plot' they go in MRM.Contrasts.Plots.Con{i}.
%
%=========================================================================%
% Copyright 2016 Martyn McFarquhar                                        %
%=========================================================================%
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
%=========================================================================%
function MRM = MRM_contrastSpecifier(MRM, action)

screen = get(0, 'ScreenSize');
gui.contrastsWin = figure('Position', [screen(3)/3,screen(4)/3,1/2*screen(3),1/2.5*screen(4)],  ...
                          'DockControls', 'off', 'MenuBar', 'none', 'Visible', 'off',           ...
                          'Name', 'Contrast specifier', 'NumberTitle', 'off',                   ...
                          'CloseRequestFcn', @win_close, 'ResizeFcn', @win_resize);
  
gui.MRM = MRM;
guidata(gcf,gui);

create_widgets();
position_widgets();
setup_weights_tables();

set(gui.contrastsWin, 'Visible', 'on');
uiwait(gcf) % Pause until the add button is pressed, once it is resume the code below

gui = guidata(gcf);

if gui.cancel ~= 1
    switch action
        case 'plot'
            if ~isfield(gui.MRM.Contrasts, 'Plots')
                gui.MRM.Contrasts.Plots = [];
                gui.MRM.Contrasts.Plots.Number = 1;
                gui.MRM.Contrasts.Plots.Con{1}.Name = gui.conNameBox.String;
                gui.MRM.Contrasts.Plots.Con{1}.A    = remove_zero_rows(cell2mat(gui.betweenSubsWeightsTable.Data));
                gui.MRM.Contrasts.Plots.Con{1}.C    = remove_zero_rows(cell2mat(gui.withinSubsWeightsTable.Data));
            else
                num = gui.MRM.Contrasts.Plots.Number+1;
                gui.MRM.Contrasts.Plots.Number        = num;
                gui.MRM.Contrasts.Plots.Con{num}.Name = gui.conNameBox.String;
                gui.MRM.Contrasts.Plots.Con{num}.A    = remove_zero_rows(cell2mat(gui.betweenSubsWeightsTable.Data));
                gui.MRM.Contrasts.Plots.Con{num}.C    = remove_zero_rows(cell2mat(gui.withinSubsWeightsTable.Data));
            end
        case 'con'
            if ~isfield(gui.MRM, 'Contrasts')
                gui.MRM.Contrasts             = [];
                gui.MRM.Contrasts.Number      = 1;
                gui.MRM.Contrasts.Con{1}.Name = gui.conNameBox.String;
                gui.MRM.Contrasts.Con{1}.A    = remove_zero_rows(cell2mat(gui.betweenSubsWeightsTable.Data));
                gui.MRM.Contrasts.Con{1}.C    = remove_zero_rows(cell2mat(gui.withinSubsWeightsTable.Data));
            else
                num = gui.MRM.Contrasts.Number+1;
                gui.MRM.Contrasts.Number        = num;
                gui.MRM.Contrasts.Con{num}.Name = gui.conNameBox.String;
                gui.MRM.Contrasts.Con{num}.A    = remove_zero_rows(cell2mat(gui.betweenSubsWeightsTable.Data));
                gui.MRM.Contrasts.Con{num}.C    = remove_zero_rows(cell2mat(gui.withinSubsWeightsTable.Data));
            end
    end
    MRM = gui.MRM;
else
    MRM = [];
end
delete(gcf);
end


%==========================================================================
% Create widgets
%==========================================================================
function create_widgets()
gui = guidata(gcf);

gui.betweenSubsWeightsTable   = uitable('BackgroundColor', [1 1 1]);
gui.betweenSubsIdentityButton = uicontrol('Style', 'pushbutton', 'String', 'Identity', 'Callback', @ident);
gui.betweenSubsAverageButton  = uicontrol('Style', 'pushbutton', 'String', 'Average',  'Callback', @average);
gui.betweenSubsSyntaxButton   = uicontrol('Style', 'pushbutton', 'String', 'Syntax',   'Callback', @eval_syntax);
gui.withinSubsWeightsTable    = uitable('BackgroundColor', [1 1 1]);
gui.withinSubsIdentityButton  = uicontrol('Style', 'pushbutton', 'String', 'Identity', 'Callback', @ident); 
gui.withinSubsAverageButton   = uicontrol('Style', 'pushbutton', 'String', 'Average',  'Callback', @average);
gui.withinSubsSyntaxButton    = uicontrol('Style', 'pushbutton', 'String', 'Syntax',   'Callback', @eval_syntax);
gui.conNameBox                = uicontrol('Style', 'edit', 'HorizontalAlignment', 'left');
gui.addConButton              = uicontrol('Style', 'pushbutton', 'String', 'Add contrast', 'Callback', @add_con);

guidata(gcf, gui);
end


%==========================================================================
% Position widgets
%==========================================================================
function position_widgets()
gui = guidata(gcf);

mainVBox   = uix.VBox('Parent', gui.contrastsWin, 'Spacing', 10, 'Padding', 10);
topHBox    = uix.HBox('Parent', mainVBox, 'Spacing', 10);
bottomHBox = uix.HBox('Parent', mainVBox);

betweenSubPanel = uipanel('Parent', topHBox, 'Title', 'Between-subject weights');
betweenSubPanelVBox = uix.VBox('Parent', betweenSubPanel, 'Spacing', 5, 'Padding', 10);
gui.betweenSubsWeightsTable.Parent = betweenSubPanelVBox;
betweenButtonsHBox = uix.HBox('Parent', betweenSubPanelVBox);
gui.betweenSubsIdentityButton.Parent = betweenButtonsHBox;
gui.betweenSubsAverageButton.Parent  = betweenButtonsHBox;
gui.betweenSubsSyntaxButton.Parent   = betweenButtonsHBox;
uix.Empty('Parent', betweenButtonsHBox);
betweenButtonsHBox.Widths = [100 100 100 -1];
betweenSubPanelVBox.Heights = [-1 30];

withinSubPanel = uipanel('Parent', topHBox, 'Title', 'Within-subject weights');
withinSubPanelVBox = uix.VBox('Parent', withinSubPanel, 'Spacing', 5, 'Padding', 10);
gui.withinSubsWeightsTable.Parent = withinSubPanelVBox;
withinButtonsHBox = uix.HBox('Parent', withinSubPanelVBox);
gui.withinSubsIdentityButton.Parent = withinButtonsHBox;
gui.withinSubsAverageButton.Parent  = withinButtonsHBox;
gui.withinSubsSyntaxButton.Parent   = withinButtonsHBox;
uix.Empty('Parent', withinButtonsHBox);
withinButtonsHBox.Widths = [100 100 100 -1];
withinSubPanelVBox.Heights = [-1 30];

uicontrol('Parent', bottomHBox, 'Style', 'text', 'String', 'Name')
gui.conNameBox.Parent   = bottomHBox;
gui.addConButton.Parent = bottomHBox;
bottomHBox.Widths       = [50 -1 100];

mainVBox.Heights = [-1 30];

guidata(gcf, gui);
end


%==========================================================================
% Table displays
%==========================================================================
function setup_weights_tables()
gui = guidata(gcf);

gui.betweenSubsWeightsTable.Data    = cell(size(gui.MRM.Design.X.X,2), size(gui.MRM.Design.X.X,2));
gui.betweenSubsWeightsTable.RowName = [];

BSName = cell(size(gui.MRM.Design.X.Cell,2),1);
for i = 1:size(gui.MRM.Design.X.Cell,2)
    BSName{i} = gui.MRM.Design.X.Cell{i}.Label;
end

gui.betweenSubsWeightsTable.ColumnName     = BSName;
gui.betweenSubsWeightsTable.ColumnEditable = logical(ones(1,size(gui.MRM.Design.X.Cell,2))); %#ok<*LOGL>
gui.betweenSubsWeightsTable.ColumnFormat   = repmat({'numeric'},1,size(gui.MRM.Design.X.Cell,2));
gui.betweenSubsWeightsTable.ColumnWidth    = {(gui.betweenSubsWeightsTable.Position(3)/size(gui.MRM.Design.X.Cell,2)-6)};

gui.withinSubsWeightsTable.Data    = cell(size(gui.MRM.Design.Y.Y,2), size(gui.MRM.Design.Y.Y,2));
gui.withinSubsWeightsTable.RowName = [];

if size(gui.MRM.Design.Y.Cell,2) == 1
    gui.withinSubsWeightsTable.ColumnName     = [];
    gui.withinSubsWeightsTable.Data           = {1};
    gui.withinSubsWeightsTable.ColumnEditable = logical(0);
else
    WSName = cell(size(gui.MRM.Design.Y.Cell,2),1);
    for i = 1:size(gui.MRM.Design.Y.Cell,2)
        WSName{i} = gui.MRM.Design.Y.Cell{i}.Label;
    end    
    gui.withinSubsWeightsTable.ColumnName      = WSName;
    gui.withinSubsWeightsTable.ColumnEditable  = logical(ones(1,size(gui.MRM.Design.Y.Cell,2))); %#ok<*LOGL>
    gui.withinSubsWeightsTable.ColumnFormat    = repmat({'numeric'},1,size(gui.MRM.Design.Y.Cell,2));
    gui.withinSubsWeightsTable.ColumnWidth    = {(gui.withinSubsWeightsTable.Position(3)/size(gui.MRM.Design.Y.Cell,2)-6)};

end

guidata(gcf,gui);
end


%==========================================================================
% Make identity matrix
%==========================================================================
function ident(obj,~)
gui = guidata(gcf);
switch obj
    case gui.betweenSubsIdentityButton
         gui.betweenSubsWeightsTable.Data = num2cell(eye(size(gui.MRM.Design.X.Cell,2)));
    case gui.withinSubsIdentityButton
         gui.withinSubsWeightsTable.Data  = num2cell(eye(size(gui.MRM.Design.Y.Cell,2)));
end
guidata(gcf,gui);
end


%==========================================================================
% Make average matrix
%==========================================================================
function average(obj,~)
gui = guidata(gcf);
switch obj
    case gui.betweenSubsAverageButton
        A = repmat(1/size(gui.MRM.Design.X.Cell,2),1,size(gui.MRM.Design.X.Cell,2));
        A = padarray(A,[size(gui.MRM.Design.X.Cell,2)-size(A,1) size(gui.MRM.Design.X.Cell,2)-size(A,2)], 0, 'post');
        gui.betweenSubsWeightsTable.Data = num2cell(A);
    case gui.withinSubsAverageButton
        A = repmat(1/size(gui.MRM.Design.Y.Cell,2),1,size(gui.MRM.Design.Y.Cell,2));
        A = padarray(A,[size(gui.MRM.Design.Y.Cell,2)-size(A,1) size(gui.MRM.Design.Y.Cell,2)-size(A,2)], 0, 'post');
        gui.withinSubsWeightsTable.Data = num2cell(A);
end
guidata(gcf,gui);
end


%==========================================================================
% Evaluate syntax
%==========================================================================
function eval_syntax(obj,~)
gui = guidata(gcf);
syntax = inputdlg('Enter MATLAB syntax', 'Syntax');

if isempty(syntax)
    return
end

try
    A = eval(syntax{:});
catch err
    msgbox(['Error evaluating statement: ' err.message]);
    return
end

switch obj
    case gui.betweenSubsSyntaxButton
        if size(A,1) > size(gui.MRM.Design.X.Cell,2) || size(A,2) > size(gui.MRM.Design.X.Cell,2)
            msgbox('Error: Dimensions of resulting matrix are too large');
            return
        else
            A = padarray(A,[size(gui.MRM.Design.X.Cell,2)-size(A,1) size(gui.MRM.Design.X.Cell,2)-size(A,2)], 0, 'post');
            gui.betweenSubsWeightsTable.Data = num2cell(A);
            guidata(gcf,gui);
        end
    case gui.withinSubsSyntaxButton
        if size(A,1) > size(gui.MRM.Design.Y.Cell,2) || size(A,2) > size(gui.MRM.Design.Y.Cell,2)
            msgbox('Error: Dimensions of resulting matrix are too large');
            return
        else
            A = padarray(A,[size(gui.MRM.Design.Y.Cell,2)-size(A,1) size(gui.MRM.Design.Y.Cell,2)-size(A,2)], 0, 'post');
            gui.withinSubsWeightsTable.Data = num2cell(A);
            guidata(gcf,gui);
        end
end
end


%==========================================================================
% Window close
%==========================================================================
function win_close(~,~)
gui = guidata(gcf);
gui.cancel = 1; % Window X was clicked
guidata(gcf,gui);
uiresume(gcf);
end


%==========================================================================
% Window close
%==========================================================================
function add_con(~,~)
gui = guidata(gcf);

% Check name
if isempty(gui.conNameBox.String)
   msgbox('Please enter a name')
   return
end

% Replace NaN with 0 and then pad with 0
A = cell2mat(gui.betweenSubsWeightsTable.Data);
A(isnan(A)) = 0;
Apadded = zeros(size(gui.MRM.Design.X.Cell,2), size(gui.MRM.Design.X.Cell,2));
Apadded(1:size(A,1), 1:size(A,2)) = A;
gui.betweenSubsWeightsTable.Data = num2cell(Apadded);

B = cell2mat(gui.withinSubsWeightsTable.Data);
B(isnan(B)) = 0;
Bpadded = zeros(size(gui.MRM.Design.Y.Cell,2), size(gui.MRM.Design.Y.Cell,2));
Bpadded(1:size(B,1), 1:size(B,2)) = B;
gui.withinSubsWeightsTable.Data = num2cell(Bpadded);

gui.cancel = 0;
guidata(gcf,gui);
if estimability_check(Apadded, gui.MRM.Design.X.X) == 1
    uiresume(gcf);
else
    msgbox('The between-subjects contrast is not estimable');
end
end


%==========================================================================
% Remove zero rows
%==========================================================================
function mat = remove_zero_rows(mat)
for i = size(mat,1):-1:1
    if sum(mat(i,:) == 0) == size(mat,2)
        mat(i,:) = [];
    end
end
end


%==========================================================================
% Window resize
%==========================================================================
function win_resize(~,~)
gui = guidata(gcf);
gui.betweenSubsWeightsTable.ColumnWidth = {(gui.betweenSubsWeightsTable.Position(3)/size(gui.MRM.Design.X.Cell,2)-6)};
gui.withinSubsWeightsTable.ColumnWidth  = {(gui.withinSubsWeightsTable.Position(3)/size(gui.MRM.Design.Y.Cell,2)-6)};
guidata(gcf,gui);
end


%==========================================================================
% Uses the SAS algorithm detailed at https://support.sas.com/
%==========================================================================
function est = estimability_check(con, X)
H = pinv(X'*X)*X'*X;
con = remove_zero_rows(con);
LH = abs(con - con*H);
C = abs(con);
C(C == 0) = 1;
C = C * 10e-4;
est = 1;
for i = 1:size(LH,1)
    for j = 1:size(LH,2)
        if LH(i,j) > C(i,j)
            est = 0;
            return
        end
    end
end
end
