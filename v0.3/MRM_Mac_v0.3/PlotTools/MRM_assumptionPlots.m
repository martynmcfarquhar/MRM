function MRM_assumptionPlots(voxelData, MNIcoords, whichPlot, DVs)
%==========================================================================
% Standard assumption plots for a multivariate GLM
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This function creates a number of plots that can be used to assess how
% well the assumptions of a standard multivariate GLM are met. We provide
% the following:
%   - Histogram of the model residuals - normality
%   - QQ plot of the model resisuals   - normality + outliers
%   - Fitted vs. residual scatter      - homogeneity of variance + outliers
%   - Box plots of the raw data        - normality + outliers
%   - Scatter plots of DVs             - covariance homogeneity
%
% DVs : Flag - 0 if we are plotting all DVs, else it is the number of the
%              DV
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

%==========================================================================
% Re-estimate the model (faster than reading in the images)
%==========================================================================
X         = MRM.Design.X.X;
tempPE    = inv(X' * X) * X' * voxelData;
tempFit   = X * tempPE;
tempResid = voxelData - tempFit;
tempVCOV  = (tempResid' * tempResid) / (size(X,1) - size(X,2));

%==========================================================================
% Residual histograms
%==========================================================================
if whichPlot == 2
   
    if DVs == 0
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.2], ...
               'Units',         'Pixels', ...
               'NumberTitle',   'off',    ...
               'DockControls',  'off', ...
               'ToolBar',       'none');
    else
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.6], ...
               'Units',         'Pixels', ...
               'NumberTitle',   'off',    ...
               'DockControls',  'off', ...
               'ToolBar',       'none');
    end
    
   nPlots = size(tempResid,2) + 1; 
   nCols = 2;

   % Vague rule to try and keep the plots tidy
   if round(nPlots/nCols) > 4
       nCols = 3;
   end
   
   nRows = round(nPlots/nCols);
    
   for k = 1:size(tempResid,2)
      
       %-------------------------------------------------------------------
       % Plot the residual histogram
       %-------------------------------------------------------------------
       if DVs == 0
            subplot(nRows, nCols, k); 
       else
            subplot(size(tempResid,2), 1, k); 
       end
       
       hold on
       
       if license('test', 'statistics_toolbox') == 1 
            % Freedman-Diaconis rule
            bins = (max(tempResid(:,k)) - min(tempResid(:,k))) / (2 * iqr(tempResid(:,k)) / size(tempResid(:,k),1)^(1/3));
            counts = hist(tempResid(:,k), bins);
            hist(tempResid(:,k), bins);
       else
            counts = hist(tempResid(:,k));
            hist(tempResid(:,k));
       end
       h = findobj(gca,'Type','patch');
       set(h,'FaceColor',[1 1 1]); 
       
       %-------------------------------------------------------------------
       % Overlay a normal distribution curve with mean = 0 and variance
       % given by the i,i element of tempVCOV
       %-------------------------------------------------------------------
       samps   = spm_normrnd(0, tempVCOV(k,k), 1000);
       probs   = spm_Npdf(samps, 0, tempVCOV(k,k));
       scaling = max(counts) / max(probs);
       probs   = probs * scaling;
       mat     = [samps' probs'];
       mat     = sortrows(mat, 1);
       plot(mat(:,1), mat(:,2));
       
       %-------------------------------------------------------------------
       % Plot title giving the DV name and estimated variance
       %-------------------------------------------------------------------
       
       if DVs == 0
            title(['\bf ' MRM.Design.Y.Cell{k}.Label ' (Var = ' ...
                   num2str(tempVCOV(k,k)) ')'], ...
                  'FontSize', 12);
       else
           
           title(['\bf ' MRM.Design.Y.Cell{DVs}.Label ' (Var = ' ...
                   num2str(tempVCOV(k,k)) ')'], ...
                  'FontSize', 12); 
       end
       
       hold off;
       
   end
   
   %-----------------------------------------------------------------------
   % One linear combination of the DVs (the sum)
   %-----------------------------------------------------------------------
   if DVs == 0
       
       subplot(nRows, 2, size(tempResid,2) + 1);
       hold on
       
       newResids = sum(tempResid,2);
       
       if license('test', 'statistics_toolbox') == 1   
            % Freedman-Diaconis rule
            bins = (max(newResids) - min(newResids)) / (2 * iqr(newResids) / size(newResids,1)^(1/3));
            counts = hist(newResids, bins);
            hist(newResids, bins);
       else
            counts = hist(tempResid(:,k));
            hist(tempResid(:,k));
       end
       
       h = findobj(gca,'Type','patch');
       set(h,'FaceColor',[1 1 1]);

       v =  (newResids' * newResids) / (size(X,1) - size(X,2));
       
       % This is equivalent to summing all elements in the full covariance
       % matrix given that the variance of a sum is the sum of all variances
       % and covariances

       samps   = spm_normrnd(0, v, 1000);
       probs   = spm_Npdf(samps, 0, v);
       scaling = max(counts) / max(probs);
       probs   = probs * scaling;
       mat     = [samps' probs'];
       mat     = sortrows(mat, 1);
       plot(mat(:,1), mat(:,2));

       title(['\bf Sum of the residuals (Var = ' ...
             num2str(v) ')'], ...
           'FontSize', 12);
   end
   
   hold off;
   
   %-----------------------------------------------------------------------
   % Main plot title
   %-----------------------------------------------------------------------
   axes('Position',[0 0 1 1], ...
       'Xlim',[0 1], ...
       'Ylim',[0 1], ...
       'Box','off', ...
       'Visible','off', ...
       'Units','normalized', ...
       'clipping' , 'off');
   
   text(0.5, 0.98, ...
       ['Residual histogram(s) at ' num2str(MNIcoords(1)) ...
       ' ' num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3))], ...
       'HorizontalAlignment','center', ...
       'VerticalAlignment', 'top');  
end

%==========================================================================
% Residual QQ normal plot
%==========================================================================
if whichPlot == 3
    
    if DVs == 0
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.2], ...
            'Units',         'Pixels', ...
            'NumberTitle',   'off',    ...
            'DockControls',  'off',    ...
            'ToolBar',       'none');
    else
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.6], ...
            'Units',         'Pixels', ...
            'NumberTitle',   'off',    ...
            'DockControls',  'off',    ...
            'ToolBar',       'none');
    end
    
    nPlots = size(tempResid,2) + 1;
    nCols  = 2;
    
    % Vague rule to try and keep the plots tidy
    if round(nPlots/nCols) > 4
        nCols = 3;
    end
   
    nRows  = round(nPlots/nCols);
    
    for k = 1:size(tempResid,2)
        
        %-------------------------------------------------------------------
        % Plot the normal QQ plot
        %-------------------------------------------------------------------
        if DVs == 0
            subplot(nRows, 2, k);
            qqplot(tempResid(:,k));
            title(MRM.Design.Y.Cell{k}.Label);
        else
            subplot(nRows, 1, k);
            qqplot(tempResid(:,k));
            title(MRM.Design.Y.Cell{DVs}.Label);
        end
        
        ylabel('Sample quantiles');
        xlabel('');

    end
    
    %-----------------------------------------------------------------------
    % One linear combination of the DVs (the sum)
    %-----------------------------------------------------------------------
    if DVs == 0
        subplot(nRows, 2, size(tempResid,2) + 1);
        newResids = sum(tempResid,2);
        qqplot(newResids);
        ylabel('Sample quantiles');
        xlabel('Theoretical Quantiles');
        title('Sum of the DVs');
    end
        
    %-----------------------------------------------------------------------
    % Main plot title
    %-----------------------------------------------------------------------
    axes('Position',[0 0 1 1], ...
        'Xlim',[0 1], ...
        'Ylim',[0 1], ...
        'Box','off', ...
        'Visible','off', ...
        'Units','normalized', ...
        'clipping' , 'off');
    
    text(0.5, 0.98, ...
        ['Normal QQ plot(s) at ' num2str(MNIcoords(1)) ...
        ' ' num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3))], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', 'top');
end


%==========================================================================
% Fitted vs. residual scatter plot
%==========================================================================
if whichPlot == 4
   
    if DVs == 0
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.2], ...
            'Units',         'Pixels', ...
            'NumberTitle',   'off',    ...
            'DockControls',  'off',    ...
            'ToolBar',       'none');
    else
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.6], ...
            'Units',         'Pixels', ...
            'NumberTitle',   'off',    ...
            'DockControls',  'off',    ...
            'ToolBar',       'none');
    end
    
    nPlots = size(tempResid,2);
    nCols  = 2;
    
    if round(nPlots/2) > 4
        nCols = 3;
    end
    
    nRows = round(nPlots/nCols);
    
    for k = 1:size(tempResid,2)
       
        if DVs == 0
            subplot(nRows, 2, k);
        else
            subplot(nRows, 1, k);
        end
        
        scatter(tempFit(:,k), tempResid(:,k), 'black');
        
        % Add line at 0
        line('XData', xlim, 'YData', [0 0], 'LineStyle', '-.');
        
        if DVs == 0
            title(MRM.Design.Y.Cell{k}.Label, ...
                 'FontSize', 12);
        else
            title(MRM.Design.Y.Cell{DVs}.Label, ...
                 'FontSize', 12);
        end
        
        ylabel('Residual values', ...
               'FontSize', 12);
           
       if k == size(tempResid,2)
           xlabel('Fitted values', ...
               'FontSize', 12);
       end
        
    end
    
    %-----------------------------------------------------------------------
    % Main plot title
    %-----------------------------------------------------------------------
    axes('Position',[0 0 1 1], ...
        'Xlim',[0 1], ...
        'Ylim',[0 1], ...
        'Box','off', ...
        'Visible','off', ...
        'Units','normalized', ...
        'clipping' , 'off');
    
    text(0.5, 0.98, ...
        ['Fitted vs. residual values at ' num2str(MNIcoords(1)) ...
        ' ' num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3))], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', 'top');
end



%==========================================================================
% Box plots of the groups at each DV
%==========================================================================
if whichPlot == 5
    
    nDVs = size(voxelData,2);
    
    % How many groups?
    groups = size(MRM.Design.X.X,2) - MRM.Covariates.Number;
    
    Xnan            = MRM.Design.X.X(:,1:groups);
    Xnan(Xnan == 0) = NaN;
    
    if DVs == 0
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.2], ...
            'Units',         'Pixels', ...
            'NumberTitle',   'off',    ...
            'DockControls',  'off',    ...
            'ToolBar',       'none');
    else
        figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.6], ...
            'Units',         'Pixels', ...
            'NumberTitle',   'off',    ...
            'DockControls',  'off',    ...
            'ToolBar',       'none');
    end
       
    % Group mapping of the data
    %groupNum = cell(size(voxelData, 1), 1);
    ind = zeros(size(voxelData, 1), 1);
    
    for i = 1:size(voxelData, 1)  
        %groupNum{i} = num2str(find(Xnan(i,:) == 1)); 
        ind(i,1) = find(Xnan(i,:) == 1);
    end
    
    % Labels
    labs = cell(size(MRM.Design.X.Cell,2),1);
    for i = 1:size(MRM.Design.X.Cell,2)
       
        labs{i} = MRM.Design.X.Cell{i}.Label;
        
    end
    
    nCols = 2;
    
    if round(nDVs/2) > 4
        nCols = 3;
    end
    
    nRows = round(nDVs/nCols);
       
    for i = 1:nDVs
        
        if DVs == 0
            subplot(nRows, 2, i);
        else
            subplot(nRows, 1, i);
        end

        % Scatter plot
        Ynew = voxelData(:,i);
        %ind = cell2mat(groupNum);
        %ind = str2double(groupNum{:});
        scatter(ind, Ynew, 'black', 'jitter','on', 'jitterAmount',0.1);
        hold on;
        
        % Box plot
        Ynew = voxelData(:,i);
        %boxplot(Ynew, groupNum);
        boxplot(Ynew, ind);
        
        hold off;
        
        if DVs == 0
            title(['\bf ' MRM.Design.Y.Cell{i}.Label], ...
                   'FontSize', 12);
        else
            title(['\bf ' MRM.Design.Y.Cell{DVs}.Label], ...
                   'FontSize', 12);
        end
        
        set(gca,'xtickmode','auto','xticklabelmode','auto');
        set(gca,'XTick', 1:size(MRM.Design.X.Cell,2));
        set(gca,'XTickLabel', labs);
    end
    
    %----------------------------------------------------------------------
    % Main title
    %----------------------------------------------------------------------
    axes('Position',[0 0 1 1], ...
        'Xlim',[0 1], ...
        'Ylim',[0 1], ...
        'Box','off', ...
        'Visible','off', ...
        'Units','normalized', ...
        'clipping' , 'off');
    
    text(0.5, 0.98, ...
        ['Box and jittered scatter plot(s) at ' num2str(MNIcoords(1)) ...
        ' ' num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3))], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', 'top');
    
end



%==========================================================================
% Scatter plots of DV pairs
%==========================================================================
if whichPlot == 6
   
    nDVs = size(voxelData,2);
    
    if nDVs == 1
        msgbox('Cannot draw scatter plots of DV pairs for only 1 DV.')
        return
    end
    
    % How many columns will we need?
    combos  = flipdim(combnk(1:nDVs,2), 1);
    columns = size(combos,1);
    
    % How many rows (groups)?
    rows = size(MRM.Design.X.X,2) - MRM.Covariates.Number;
    
    Xnan            = X;
    Xnan(Xnan == 0) = NaN;
    
    figure('Position',      [win(3)/4, win(4)/4, win(3)/1.2, win(4)/1.2], ...
           'Units',         'Pixels', ...
           'NumberTitle',   'off',    ...
           'DockControls',  'off',    ...
           'ToolBar',       'none');
    
    ind = 1;
    
    for j = 1:rows
        
        for k = 1:columns
            
            %--------------------------------------------------------------
            % For each grouping separate out the column of the Y matrix and
            % the residuals to form the pair for plotting
            %--------------------------------------------------------------
            Ynew   = [voxelData(:,combos(k,1)) voxelData(:,combos(k,2))];
            ResNew = [tempResid(:,combos(k,1)) tempResid(:,combos(k,2))];
            
            subplot(rows, columns, ind);
            
            Ynew(:,1) = Ynew(:,1) .* Xnan(:,j);
            Ynew(:,2) = Ynew(:,2) .* Xnan(:,j);
            Ynew(any(isnan(Ynew),2),:) = [];
            
            ResNew(:,1) = ResNew(:,1) .* Xnan(:,j);
            ResNew(:,2) = ResNew(:,2) .* Xnan(:,j);
            ResNew(any(isnan(ResNew),2),:) = [];
            
            scatter(Ynew(:,1), Ynew(:,2), '+', 'black');
            
            %axis square
            hold on
            
            %--------------------------------------------------------------
            % Estimate the variance-covariance matrix for this subset
            %--------------------------------------------------------------
            covMat = (ResNew' * ResNew) / (size(ResNew,1) - size(X,2));
            
            elpt = ellipsedata(covMat, ...
                               [mean(Ynew(:,1)), mean(Ynew(:,2))], ...
                               1000, ...
                               2); % 2 standard deviations
            
            plot(elpt(:,1:2:end),elpt(:,2:2:end), 'color', 'black');
            
            %-------------------------------------------------------------------
            % Plot title giving the DV name and estimated variance
            %-------------------------------------------------------------------
            if j == 1
                title(['\bf ' MRM.Design.Y.Cell{combos(k,1)}.Label ' vs ' ...
                    MRM.Design.Y.Cell{combos(k,2)}.Label]);
            end
            
            if k == 1
                ylabel(MRM.Design.X.Cell{j}.Label);
            end
            
            hold off;
            
            ind = ind + 1;
            
        end 
    end
    %-----------------------------------------------------------------------
    % Main plot title
    %-----------------------------------------------------------------------
    axes('Position',[0 0 1 1], ...
        'Xlim',[0 1], ...
        'Ylim',[0 1], ...
        'Box','off', ...
        'Visible','off', ...
        'Units','normalized', ...
        'clipping' , 'off');
    
    text(0.5, 0.98, ...
        ['Scatter plots of DV pairs at ' num2str(MNIcoords(1)) ...
        ' ' num2str(MNIcoords(2)) ' ' num2str(MNIcoords(3)) ...
        '.'], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', 'top');
    
    text(0.5, 0.07, ...
        ['Ellipses represent the countour curve at 2 standard deviations ' ...
        'of the bivariate normal distribution implied by the data.'], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', 'top');
            
end

end