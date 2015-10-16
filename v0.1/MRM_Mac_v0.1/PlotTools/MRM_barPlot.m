function MRM_barPlot(A, beta, C, Sigma, PlotCI, location, name, save, ...
                     saveWhere, CILevel)
%==========================================================================
% Bar and error plot for a linear combination of parameters in a
% multivariate GLM
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This function creates a bar plot with error bars for the linear
% combination parameters provided in a contrast. It also tries to label the
% bars sensibly given the labels in an MRM structure.
%
%=========================================================================%
% Copyright 2015 Martyn McFarquhar                                        %
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

global MRM

win = get(0, 'ScreenSize');
X   = MRM.Design.X.X;

aDim = size(A,1);
cDim = size(C,1);

means = zeros(1, aDim * cDim);
vars = zeros(1, aDim * cDim);

k = 1;

%--------------------------------------------------------------------------
% Get labels for the plot x-axis
%--------------------------------------------------------------------------

% Design Matrix labels
Xlabels = cell(1, size(MRM.Design.X.Cell, 2) + MRM.Covariates.Number);

for i = 1:size(MRM.Design.X.Cell, 2)
    Xlabels{i} = MRM.Design.X.Cell{i}.Label;
end

if MRM.Covariates.Number ~= 0
    
    if isfield(MRM.Covariates, 'CV')
    
        for i = 1:MRM.Covariates.Number
           Xlabels{size(MRM.Design.X.Cell, 2) + i} = MRM.Covariates.CV{i}.Name; 
        end
    
    end
end

% Column labels Y
Ylabels = cell(1, size(MRM.Design.Y.Cell, 2));

for i = 1:size(MRM.Design.Y.Cell, 2)
    Ylabels{i} = MRM.Design.Y.Cell{i}.Label;
end


%--------------------------------------------------------------------------
% Work through each row combination of A and C, saving the results in
% the vector 'means', 'vars', and 'label'
%--------------------------------------------------------------------------
label = cell(1, aDim * cDim);
invXX = inv(X' * X);

for i = 1:aDim
    for j = 1:cDim
        
        a = A(i, :);
        c = C(j, :);
        
        means(k) = a * beta * c';
        vars(k)  = (c * Sigma * c') * (a * invXX * a');
        
        %------------------------------------------------------------------
        % Between-subjects labels
        %------------------------------------------------------------------
        tempLab = '';
        
        for l = 1:size(a,2)
            if a(:,l) ~= 0
                tempLab = [tempLab ' ' Xlabels{l}];
            end
        end
        
        str        = strsplit(tempLab);
        uniqueLabs = unique(str);
        
        for l = 1:size(uniqueLabs,2)
            label{k} = [label{k} ' ' uniqueLabs{l}];
        end
        
        %------------------------------------------------------------------
        % Within-subjects labels
        %------------------------------------------------------------------
        tempLab = '';
        
        for l = 1:size(c,2)
            if c(:,l) ~= 0
                tempLab = [tempLab ' ' Ylabels{l}];
            end
        end
        
        str        = strsplit(tempLab);
        uniqueLabs = unique(str);
        
        for l = 1:size(uniqueLabs,2)
            label{k} = [label{k} ' ' uniqueLabs{l}];
        end
        
        k = k + 1;
        
    end
end

% Try and construct sensible labels
for i = 1:size(label,2)
   
    str    = strsplit(label{i});
    str(1) = [];
    
    idx = zeros(size(str,2),1);
    
    % Between
    for j = 1:size(str,2)
       
        for k = 1:MRM.Factors.Between.Number
           
            for l = 1:MRM.Factors.Between.Factors{k}.LevelsNum
               
                if strcmp(MRM.Factors.Between.Factors{k}.Levels{l}, str(j)) == 1
                   
                    idx(j) = k;
                    
                end
            end
        end
    end
    
    % Within
    for j = 1:size(str,2)
       
        for k = 1:MRM.Factors.Within.Number
           
            for l = 1:MRM.Factors.Within.Factors{k}.LevelsNum
               
                if strcmp(MRM.Factors.Within.Factors{k}.Levels{l}, str(j)) == 1
                   
                    idx(j) = MRM.Factors.Between.Number + k;
                    
                end
            end
        end
    end
    
    for j = 1:max(idx)

        if size(str(idx == j),2) > 1
            str(idx == j) = [];
            idx(idx == j) = [];
        end
    end
    
    label{i} = '';
    
    for j = 1:size(str,2)
        label{i} = [label{i} ' ' str{j}];
    end
   
end

%--------------------------------------------------------------------------
% Standard error of the contrast
%--------------------------------------------------------------------------
se = sqrt(vars);

%--------------------------------------------------------------------------
% Asymptotic confidence interval
%--------------------------------------------------------------------------
alpha = 1 - CILevel;
zVal  = spm_invNcdf(1 - (alpha/2));
CI    = zVal .* se;

%--------------------------------------------------------------------------
% Are we drawing or are we saving the values?
%--------------------------------------------------------------------------
if save == 0
    
    if isempty(findobj('type','figure','name','MRM Bar Plot')) == 1
    
        h = figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/2], ...
                   'Units',         'Pixels', ...
                   'DockControls',  'off', ...
                   'NumberTitle',   'off', ...
                   'ToolBar',       'none', ...
                   'Name',          'MRM Bar Plot');
    else
        
        h = findobj('type','figure','name','MRM Bar Plot');
        set(0, 'currentfigure', h);
        
    end
       
    bar(means, 'facecolor', [0, 0.75, 0.75]);
    hold on;
    
    if PlotCI == 1
        errorb(means, CI, 'linewidth', 1);
        errorText = [' (' num2str(CILevel * 100) '% confidence intervals)'];
    else
        errorb(means, se, 'linewidth', 1);
        errorText = ' (standard error)';
    end
    
    set(0, 'currentfigure', h);
    
    title(['\bf' name ' at ' num2str(location(1)) ' ' num2str(location(2)) ...
           ' ' num2str(location(3)) errorText], ...
           'FontSize', 12);
      
    set(gca, 'XTick', [1:size(means, 2)]);   
    set(gca, 'XTickLabel', label, 'FontSize', 11);
    rotateXLabels(gca, 45);
       
    hold off;
    
else
    
    nameNoSpace = regexprep(name, ' ', '_');
    
    fileID = fopen([saveWhere filesep nameNoSpace '_plot_values_' ...
                    num2str(location(1)) '_' num2str(location(2)) '_' ...
                    num2str(location(3)) '.txt'], 'w');
    
    outTable = [means', vars', se', (zVal .* se')];
    
    displaytable(outTable, ...
                 {'Mean', 'Variance', 'Standard error', [num2str(CILevel * 100) '% CIs']}, ...
                 [7 8 14 7], ...
                 {'.3f','.3f','.3f','.3f'}, ...
                 [],    ...
                 fileID);
    
    fclose(fileID);
    
end

end