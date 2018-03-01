function maximaMat = MRM_getMaxima()
%=========================================================================
% Load the matrix of maxima for the currently selected contrast in the MRM 
% Post-estimation Tools.
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% This function basically reads in the results text file for the currently
% selected contrast and saves it in a matrix called maximaMat containing
% the coordinates and stat value.
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

h = findobj('type','figure','name','MRM Post-estimation Tools');

if isempty(h)
    msgbox('MRM window handle not found')
    return
end

handles = guidata(h);

con  = get(handles.contrastExploreMenu, 'Value') - 1;
name = regexprep(MRM.Contrasts.Con{con}.Name, ' ', '_');

if strcmp(MRM.Options.Thresh.Pvals, 'Permutation') == 1
    filepath = [MRM.Options.Out filesep 'PermutationContrasts' filesep name filesep ...
               'MRM_' name '_Results.txt'];
else
    filepath = [MRM.Options.Out filesep 'Contrasts' filesep name filesep ...
               'MRM_' name '_Results.txt'];
end

if strcmp(MRM.Options.Thresh.Level, 'Cluster') == 1

    fileID = fopen(filepath);
    input  = textscan(fileID,'%s %d %d %.3f %d %d %d', ...
                      'Delimiter','|',                 ...
                      'HeaderLines',1);

    maximaMat = zeros(size(input{1},1), 4);

    for i = 1:size(input{1},1)
        
        maximaMat(i,1) = input{3}(i);
        
        ind = 2;
        for j = 5:7
            maximaMat(i,ind) = input{j}(i);
            ind = ind + 1;
        end
    end

else

    fileID = fopen(filepath);
    input  = textscan(fileID,'%s %d %d %.3f %.3f %d %d %d', ...
                     'Delimiter','|', ...
                     'HeaderLines',1);

    maximaMat = zeros(size(input{1},1), 4);

    for i = 1:size(input{1},1)
        
        maximaMat(i,1) = input{4}(i);
        
        ind = 2;
        for j = 6:8
            maximaMat(i,ind) = input{j}(i);
            ind = ind + 1;
        end
    end
end

fclose(fileID);

end