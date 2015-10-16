function designOK = MRM_checkSpecification(checkWhat, betOrWith)
%==========================================================================
% Check design specification
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This function checks whether all the factor specifications in a global MRM
% structure are complete i.e. no missing names or labels. It returns a 1 if
% everything is fine and a 0 otherwise
%
% checkWhat - either 'Factors' or 'All' depending on whether we want to 
%             check the data files as well as the factor structure
% betOrWith - 'between', 'within', or 'both' depending on which factor
%             structure to check.
% 
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

% First check whether a MRM.Factors structure is empty
if isempty(MRM.Factors) == 1
    designOK = 0;
    return
end

if MRM.Factors.Between.Number == 0 && MRM.Factors.Within.Number == 0
    designOK = 0;
    return
end

% Start by checking the between subject specification
if isempty(MRM.Factors.Between) == 0

    facNum = MRM.Factors.Between.Number;

    if facNum ~= 0

        if isempty(MRM.Factors.Between.Factors) == 1
            if strcmp(betOrWith, 'within') ~= 1
                designOK = 0;
            else
                designOK = 1;
            end
            
            return
        end

        if size(MRM.Factors.Between.Factors,2) ~= facNum
            designOK = 0;
            return
        end

        for i = 1:facNum

            % Check a name exists for factor i
            if isfield(MRM.Factors.Between.Factors{i}, 'Name') == 1 
                if isempty(MRM.Factors.Between.Factors{i}.Name) == 1
                    designOK = 0;
                    return
                end
            else
                designOK = 0;
                return
            end

            % Check the number of levels exists for factor i
            if isfield(MRM.Factors.Between.Factors{i}, 'LevelsNum') == 1   
                if isempty(MRM.Factors.Between.Factors{i}.LevelsNum) == 1
                    designOK = 0;
                    return
                end
            else
                designOK = 0;
                return
            end

            % Check a name exists for each level of factor i   
            if isfield(MRM.Factors.Between.Factors{i}, 'Levels') == 1   
                
                a = MRM.Factors.Between.Factors{i}.LevelsNum;
                
                if isempty(MRM.Factors.Between.Factors{i}.Levels) ~= 1
                    
                    for j = 1:a
                        if isempty(MRM.Factors.Between.Factors{i}.Levels{j}) == 1
                            designOK = 0;
                            break
                        end 
                    end
                else
                    designOK = 0;
                    return
                end
            else
                designOK = 0;
                return
            end
        end
    end
end

designOK = 1;

if strcmp(betOrWith, 'within') == 1 || strcmp(betOrWith, 'both')

    % Next check the within subject specification in the same way
    if isempty(MRM.Factors.Within) == 0

        facNum = MRM.Factors.Within.Number;

        if facNum ~= 0

            if isempty(MRM.Factors.Within.Factors) == 1
                if strcmp(betOrWith, 'between') ~= 1
                    designOK = 0;
                else
                    designOK = 1;
                end
                
                return
            end

            if size(MRM.Factors.Within.Factors,2) ~= facNum
                designOK = 0;
                return
            end

            for i = 1:facNum

                % Check a name exists for factor i
                if isfield(MRM.Factors.Within.Factors{i},'Name') == 1 

                    if isempty(MRM.Factors.Within.Factors{i}.Name) == 1

                        designOK = 0;
                        return                    
                    end

                else
                    designOK = 0;
                    return
                end

                % Check the number of levels exists for factor i    
                if isfield(MRM.Factors.Within.Factors{i}, 'LevelsNum') == 1   
                    if isempty(MRM.Factors.Within.Factors{i}.LevelsNum) == 1
                        designOK = 0;
                        return
                    end
                else
                    designOK = 0;
                    return
                end

                % Check a name exists for each level of factor i   
                if isfield(MRM.Factors.Within.Factors{i}, 'Levels') == 1   

                    a = MRM.Factors.Within.Factors{i}.LevelsNum;

                    if isempty(MRM.Factors.Within.Factors{i}.Levels) ~= 1

                        for j = 1:a
                            if isempty(MRM.Factors.Within.Factors{i}.Levels{j}) == 1;
                                designOK = 0;
                                return
                            end 
                        end
                    else
                        designOK = 0;
                        return
                    end
                else
                    designOK = 0;
                    return
                end
            end
        end
    end
end

%=========================================================================%
% Check the specified data files if requested
%=========================================================================%
if strcmp(checkWhat, 'All') == 1
  
    % First check that all the cells have files selected
    numDVs   = size(MRM.Design.Y.Cell, 2);
    numCells = size(MRM.Design.X.Cell, 2);
    
    if isempty(MRM.Data.Y) == 1
        designOK = 0;
        return;
    end
    
    if size(MRM.Data.Y, 2) ~= numDVs
        designOK = 0;
        return
    end
    
    for i = 1:numDVs
        
        for j = 1:numCells
            
            if isempty(MRM.Data.Y{i}.Cell{j}.Scans) == 1
                
                designOK = 0;
                return
            end
        end
    end
    
    % Next check that the file numbers are the same for the cells across
    % the DVs
    
    % Check against the first DV
    fileNum = cell(1, numCells);
    
    for i = 1:numCells 
        fileNum{i} = size(MRM.Data.Y{1}.Cell{i}.Scans, 2);
    end
    
    for i = 2:numDVs
        for j = 1:numCells
            
            if size(MRM.Data.Y{i}.Cell{j}.Scans, 2) ~= fileNum{j}
                
                designOK = 2;
                return
                
            end
        end
    end  
    
end

designOK = 1;

end