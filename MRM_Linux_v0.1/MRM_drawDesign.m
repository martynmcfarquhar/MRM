function MRM_drawDesign(save)
%=========================================================================
% Plot the design visualisation                 
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% This function draws the design matrix X and factorial structure of the
% outcome matrix Y and labels the columns. If the flag to save is given
% (i.e. save = 1) then this is done invisibly and saved to the output
% directory, otherwise a new window is opened containing the image.
%
% This function requires a complete MRM struct available as a global
% variable. This can be checked using MRM_checkSpecification() with the
% 'All' flag. This function also requires the rotateXLabels() function 
% availablefrom: http://www.mathworks.co.uk/matlabcentral/fileexchange/45172
% to provide functionality with earlier versions of MATLAB.
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

% Check the design matrix exists
if isfield(MRM.Design, 'X') == 1
    if isfield(MRM.Design.X, 'X') == 1
        if isempty(MRM.Design.X.X) == 1  
            msgbox('You need to select all the data before you can view the design.');
            return
        end
    else
        msgbox('You need to select all the data before you can view the design.');
        return
    end
else  
    msgbox('You need to select all the data before you can view the design.');
    return
end
    

% Columns labels X
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


win = get(0, 'ScreenSize');

if save == 1
    figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.5], ...
           'Units',         'Pixels', ...
           'Color',         [1 1 1], ...
           'NumberTitle',   'off', ...
           'Visible',       'off');
else
    figure('Position',      [win(3)/4, win(4)/4, win(3)/2.5, win(4)/1.5], ...
           'Units',         'Pixels', ...
           'Color',         [1 1 1], ...
           'Name',          'Design', ...
           'MenuBar',       'None', ...
           'DockControls',  'off', ...
           'NumberTitle',   'off', ...
           'ToolBar',       'none', ...
           'Resize',        'on');    
end



%-------------------------------------------------------------------------%
% Image of Y
%-------------------------------------------------------------------------%
switch MRM.Model
    case 'Repeated'
        
        if MRM.Factors.Within.Number > 0
            subplot(1,2,1);
            imagesc(MRM.Design.Y.Y);
            %caxis([0 size(MRM.Design.Y.Cell, 2)]);
            set(gca,'xtick', 1:size(MRM.Design.Y.Cell, 2));
            set(gca,'xticklabel', Ylabels);
            set(gca,'FontSize', 12);
            rotateXLabels(gca(), 45);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);
            colormap('gray');
            title('Factorial structure of $\mathbf{Y}$', ...
                'interpreter', 'latex', ...
                'FontSize', 15, ...
                'Position', [0.5 1.01 0]);
            
        end
        
    case 'MANOVA'
        
        if MRM.Factors.Within.Factors{1}.LevelsNum > 1
            subplot(1,2,1);
            imagesc(MRM.Design.Y.Y);
            %caxis([0 size(MRM.Design.Y.Cell, 2)]);
            set(gca,'xtick', 1:size(MRM.Design.Y.Cell, 2));
            set(gca,'xticklabel', Ylabels);
            set(gca,'FontSize', 12);
            rotateXLabels(gca(), 45);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);
            colormap('gray');
            title('Multivariate structure of $\mathbf{Y}$', ...
                'interpreter', 'latex', ...
                'FontSize', 15, ...
                'Position', [0.5 1.01 0]);
        end
end



%-------------------------------------------------------------------------%
% Image of X
%-------------------------------------------------------------------------%

X = MRM.Design.X.X;

%Scale covariates so that they visualise better
 if MRM.Covariates.Number > 0
%     
     if isfield(MRM.Covariates, 'CV')
%    
         for i = 1:MRM.Covariates.Number
% 
             colNum      = MRM.Covariates.CV{i}.Col;
%             %v           = X(:,colNum);
%             %values      = v(v ~= 0);
             values      = X(:,colNum);
             values = values ./ max(values); % scale the CV
%             values      = values - min(values(:));
%             values      = (values/(max(values(:)) - min(values(:))) * 2);
%             values      = values - 1;
%             values(values == -1) = 0;
%             %v(v~=0) = values;
%             %X(:,colNum) = v;
             X(:,colNum) = values;
% 
         end
     end
 end

switch MRM.Model
    case 'Repeated'
        
        if MRM.Factors.Within.Number > 0  
            subplot(1,2,2);
        end
        
    case 'MANOVA'
        
        if MRM.Factors.Within.Factors{1}.LevelsNum > 1 
            subplot(1,2,2); 
        end    
end

imagesc(X);
%caxis([-1 1]);
colormap(gray);
set(gca,'xtick', 1:size(MRM.Design.X.Cell, 2) + MRM.Covariates.Number);
set(gca,'xticklabel', Xlabels);
set(gca,'FontSize', 12);
rotateXLabels(gca(), 45);
colormap('gray');

title('Design matrix $\mathbf{X}$', ...
      'interpreter','latex', ...
      'FontSize', 15, ...
      'Position', [0.5 1.005 0]);

  
  
%-------------------------------------------------------------------------%
% Save
%-------------------------------------------------------------------------%
if save == 1
    w = (win(3)/2) / get(0,'ScreenPixelsPerInch');
    h = (win(4)/1.2) / get(0,'ScreenPixelsPerInch');
    
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(gcf, [MRM.Options.Out filesep 'Design.jpg'], 'jpg');
    close(gcf);
end


end