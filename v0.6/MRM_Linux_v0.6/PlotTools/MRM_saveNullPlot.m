function MRM_saveNullPlot(distribution, out, ConName, alpha, observed, clustStat)
%==========================================================================
% Save the histogram of the permutation distribution
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This function save an image of the null distribution constructed via
% permutation testing with labels for the alpha threshold and the observed
% maximum statistic added.
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

if isempty(clustStat) == 1
    lab = 'Maximum F in image';
elseif clustStat == 1
    lab = 'Maximum cluster size in image';
else
    lab = 'Largest mass in image';
end

figure('Visible','off');
    
if exist('histogram', 'builtin') == 0
    
    if license('test', 'statistics_toolbox') == 1
        
        % Freedman-Diaconis rule
        %bins = 2 * iqr(distribution) * size(distribution,1)^(-1/3);
        
        %if isinf(bins) ~= 1
        %    hist(distribution, bins);
        %else
        %    fprintf(1, 'Problem in calculating histogram bins, reverting to the default.\n');
        %    hist(distribution);
        %end
        
        histx(distribution, 'middle');
        
        threshold = quantile(distribution, 1 - alpha);
        lims      = ylim;
        
        line([threshold threshold],                     ...
            [0 lims(2)/2],                             ...
            'LineStyle', '--',                         ...
            'LineWidth', 1);
        text(threshold, (lims(2)/2) + 3, ['Threshold: ' num2str(threshold)], ...
             'Color', 'red');
        line([observed observed],                 ...
            [0 lims(2)/3],                             ...
            'LineStyle', '--',                         ...
            'LineWidth', 1);
        text(observed, (lims(2)/3) + 3, ['Observed: ' num2str(observed)],  ...
             'Color', 'red');
        xlabel(lab);
        ylabel('Count across permutations');
        saveas(gcf, [out 'PermutationContrasts' filesep ConName filesep 'nullDist.jpg'], ...
            'jpeg');
    end
    
else
    histogram(distribution);
    threshold = quantile(distribution,1 - alpha);
    lims      = ylim;
    line([threshold threshold],                              ...
        [0 lims(2)/2],                                      ...
        'LineStyle', '--',                                  ...
        'LineWidth', 1);
    text(threshold, (lims(2)/2) + 3, ['Threshold: ' num2str(threshold)],  ...
             'Color', 'red');
    line([observed observed],    ...
        [0 lims(2)/3],                                      ...
        'LineStyle', '--',                                  ...
        'LineWidth', 1);
    text(observed, (lims(2)/3) + 3, ['Observed: ' num2str(observed)],  ...
             'Color', 'red');
    xlabel(lab);
    ylabel('Count across permutations');
    saveas(gcf, [out 'PermutationContrasts' filesep ConName filesep 'nullDist.jpg'], ...
        'jpeg');
end

end