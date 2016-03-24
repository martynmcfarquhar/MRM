function MRM_designOrthogonality(save)
%=========================================================================
% Plot the orthogonality visualisation            
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% This is based entirely on the SPM function of the same purpose, where the
% help indicates:
%
%     For each pair of columns of the design matrix, the
%     orthogonality matrix depicts the magnitude of the cosine of the
%     angle between them, with the range 0 to 1 mapped white to
%     black. Orthogonal vectors (shown in white), have cosine of zero.
%     Colinear vectors (shown in black), have cosine of 1 or -1.
%
%     The cosine of the angle between two vectors a & b is obtained by
%     dividing the dot product of the two vectors by the product of
%     their lengths:
%
%     a' * b / sqrt(sum(a.^2) * sum(b.^2))
%
%     If (and only if) both vectors have zero mean, i.e.
%     sum(a)==sum(b)==0, then the cosine of the angle between the
%     vectors is the same as the correlation between the two variates.
%
% If save == 1 then the design is not shown but is simply saved to the
% folder in MRM.Options.Out
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

if nargin < 1
    save = 0;
end

%--------------------------------------------------------------------------
% Check the design matrix exists
%--------------------------------------------------------------------------
if isfield(MRM.Design, 'X') == 1
    if isfield(MRM.Design.X, 'X') == 1
        if isempty(MRM.Design.X.X) == 1  
            msgbox('You need to select all the data before you can view the design orthogonality.');
            return
        end
    else
        msgbox('You need to select all the data before you can view the design orthogonality.');
        return
    end
else  
    msgbox('You need to select all the data before you can view the design orthogonality.');
    return
end

X = MRM.Design.X.X;
    
%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% Calculate the cosines
%--------------------------------------------------------------------------
cols      = size(X,2);
orthogMat = NaN(cols, cols);

for i = 1:cols
    for j = 1:cols
        
        a = X(:,i);
        b = X(:,j);
        
        orthogMat(i,j) =  a' * b / (sqrt(sum(a.^2) * sum(b.^2)));
        
    end
end

% Use the absolute value for visualisation
orthogMat = abs(orthogMat);


%--------------------------------------------------------------------------
% Draw the lower triangle
%--------------------------------------------------------------------------
win = get(0, 'ScreenSize');

if save == 1
    figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.5], ...
           'Units',         'Pixels', ...
           'Color',         [1 1 1], ...
           'NumberTitle',   'off', ...
           'Visible',       'off');
else
    figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.5], ...
           'Units',         'Pixels',                        ...
           'Color',         [1 1 1],                         ...
           'Name',          'Design orthogonality',          ...
           'MenuBar',       'None',                          ...
           'DockControls',  'off',                           ...
           'NumberTitle',   'off',                           ...
           'ToolBar',       'none',                          ...
           'Resize',        'on');  
end


imagesc(tril(orthogMat));
colormap(flipud(gray));
caxis([0 1]);

set(gca,'xtick',      1:size(MRM.Design.X.Cell, 2) + MRM.Covariates.Number);
set(gca,'xticklabel', Xlabels);
set(gca,'FontSize',   12);

rotateXLabels(gca(), 45);

set(gca,'ytick',      1:size(MRM.Design.X.Cell, 2) + MRM.Covariates.Number);
set(gca,'yticklabel', Xlabels);
set(gca,'FontSize',   12);
set(gca,'box',       'off');

% Add the cosine values only to the lower-triangle
for i = 1:cols
    for j = 1:cols
        if j >= i
            text(i,j,num2str(orthogMat(i,j)), 'Color', 'red', 'FontSize', 12);
        end
    end
end

% Description string
string = {'For each pair of columns in $\mathbf{X}$ this matrix shows the absolute value of the cosine of' ...
          'the angle between the vectors. Orthogonal vectors (i.e. 90$^{\circ}$) are shown in white'       ...
          '(cos($\theta$) = 0), with colinear vectors (i.e. 0$^{\circ}$) shown in black (cos($\theta$) = 1 or -1).' ...
          'Colinearities caused by covariates can often be reduced by mean-centering the' ...
          'covariate in question'};
      
xlabel(string,'interpreter','latex', 'FontSize', 12);

%-------------------------------------------------------------------------%
% Save
%-------------------------------------------------------------------------%
if save == 1
    w = (win(3)/2) / get(0,'ScreenPixelsPerInch');
    h = (win(4)/1.5) / get(0,'ScreenPixelsPerInch');
    
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(gcf, [MRM.Options.Out filesep 'Orthog.jpg'], 'jpg');
    close(gcf);
end

end