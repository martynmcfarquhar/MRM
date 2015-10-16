function designOK = MRM_createDesign
%==========================================================================
% Create the design structure in a global MRM structure
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% Create the MRM.Design.X and MRM.Design.Y structures based on the factors
% specified in MRM.Factors.Between and MRM.Factors.Within. The factor
% structure must therefore be complete (can be checked using
% MRM_checkSpecification).
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

%=========================================================================%
% BETWEEN-SUBJECTS
%=========================================================================%
% Here we are creating the factorial structure of the design matrix.
% Because we use cell-means coding we can define the design matrix
% structure in terms of the cells of the design. These cells are saved in 
% MRM.Design.X.Cell. Within each cell we save the Levels of the factors
% that define that cell as well as a descriptive label of that cell.
%
% e.g. MRM.Design.X.Cell{1}.Levels = [1 1 1]
%      MRM.Design.X.Cell{1}.Label  = 'Old Control Male'
%
% These labels are used in combination with the number of DVs to organise
% the data input.
%-------------------------------------------------------------------------%

if MRM.Factors.Between.Number > 0
    if isempty(MRM.Factors.Between.Factors) ~= 1
    
        %---- Factor info ----%

        between.nCells = 1;

        for i = 1:MRM.Factors.Between.Number

            between.nCells       = between.nCells * MRM.Factors.Between.Factors{i}.LevelsNum;
            between.Indexes{i}   = 1:MRM.Factors.Between.Factors{i}.LevelsNum;

        end   


        %---- Create the factor structure ----%

        outArgs = cell(1, MRM.Factors.Between.Number);

        [outArgs{1:MRM.Factors.Between.Number}] = ndgrid(between.Indexes{:});

        Indexes = [];

        for i = 1:size(outArgs,2)
            Indexes = [Indexes outArgs{i}(:)]; 
        end

        Indexes = sortrows(Indexes, 1:numel(between.Indexes));

        % Error checking
        if size(Indexes, 1) ~= between.nCells
            msgbox('Something went wrong in creating the design', 'createDesign Failure')
            designOK = 0;
            return
        end


        %---- Create the labels for each cell ----%

        for i = 1:size(Indexes, 1)

            row = Indexes(i,:);

            MRM.Design.X.Cell{i}.Levels = row;

            for j = 1:length(row)
                string{j} = MRM.Factors.Between.Factors{j}.Levels{row(j)};
            end

            % Collapse the string elements together to create the label
            MRM.Design.X.Cell{i}.Label = strjoin(string);

        end
        
    else
        
       MRM.Design.X.Cell{1}.Levels = 1;
       MRM.Design.X.Cell{1}.Label = 'Constant'; 
        
    end

else
    
    % If there are 0 between-subjects factors all we need is a constant
    MRM.Design.X.Cell{1}.Levels = 1;
    MRM.Design.X.Cell{1}.Label = 'Constant';
    
end


%=========================================================================%
% WITHIN-SUBJECTS
%=========================================================================%
%
% Here we are creating the factorial structure of the outcome matrix Y. We 
% define column-by-column the structure that Y takes. This is based on the
% highest-order interaction of the within-subject variables
%
%-------------------------------------------------------------------------%

switch MRM.Model
    
    case 'Repeated'

        if MRM.Factors.Within.Number > 0

            if isempty(MRM.Factors.Within.Factors) ~= 1

                %---- Factor info ----%
                within.nCells = 1;

                if isfield(MRM.Factors.Within, 'Factors') == 0
                   return 
                elseif isempty(MRM.Factors.Within.Factors) == 1
                    return
                end

                for i = 1:MRM.Factors.Within.Number

                    within.nCells     = within.nCells * MRM.Factors.Within.Factors{i}.LevelsNum;
                    within.Indexes{i} = 1:MRM.Factors.Within.Factors{i}.LevelsNum;

                end

                %---- Create the factor structure ----%
                outArgs = cell(1, MRM.Factors.Between.Number);

                [outArgs{1:MRM.Factors.Within.Number}] = ndgrid(within.Indexes{:});

                Indexes = [];

                for i = 1:size(outArgs,2)
                    Indexes = [Indexes outArgs{i}(:)]; 
                end

                Indexes = sortrows(Indexes, 1:numel(within.Indexes));

                % Error checking
                if size(Indexes, 1) ~= within.nCells
                    msgbox('Something went wrong in creating the design', 'createDesign Failure')
                    designOK = 0;
                    return
                end

                %---- Create the labels for each cell ----%

                string = [];

                for i = 1:size(Indexes, 1)

                    row = Indexes(i,:);

                    MRM.Design.Y.Cell{i}.Levels = row;

                    for j = 1:length(row)
                        string{j} = MRM.Factors.Within.Factors{j}.Levels{row(j)};
                    end

                    % Collapse the string elements together to create the label
                    MRM.Design.Y.Cell{i}.Label = strjoin(string);

                end
            else

               MRM.Design.Y.Cell{1}.Levels = 1;
               MRM.Design.Y.Cell{1}.Label = ''; 

            end
        else

            % If there are 0 within-subjects factors all we need is a constant
            MRM.Design.Y.Cell{1}.Levels = 1;
            MRM.Design.Y.Cell{1}.Label = '';

        end
        
    case 'MANOVA'
        
        if MRM.Factors.Within.Number > 0

            if isempty(MRM.Factors.Within.Factors) ~= 1

                %---- Factor info ----%
                within.nCells = 1;

                if isfield(MRM.Factors.Within, 'Factors') == 0
                   return 
                elseif isempty(MRM.Factors.Within.Factors) == 1
                    return
                end

                for i = 1:MRM.Factors.Within.Number

                    within.nCells     = within.nCells * MRM.Factors.Within.Factors{i}.LevelsNum;
                    within.Indexes{i} = 1:MRM.Factors.Within.Factors{i}.LevelsNum;

                end

                %---- Create the factor structure ----%
                outArgs = cell(1, MRM.Factors.Between.Number);

                [outArgs{1:MRM.Factors.Within.Number}] = ndgrid(within.Indexes{:});

                Indexes = [];

                for i = 1:size(outArgs,2)
                    Indexes = [Indexes outArgs{i}(:)]; 
                end

                Indexes = sortrows(Indexes, 1:numel(within.Indexes));

                % Error checking
                if size(Indexes, 1) ~= within.nCells
                    msgbox('Something went wrong in creating the design', 'createDesign Failure')
                    designOK = 0;
                    return
                end

                %---- Create the labels for each cell ----%

                string = [];

                for i = 1:size(Indexes, 1)

                    row = Indexes(i,:);

                    MRM.Design.Y.Cell{i}.Levels = row;

                    for j = 1:length(row)
                        string{j} = MRM.Factors.Within.Factors{j}.Levels{row(j)};
                    end

                    % Collapse the string elements together to create the label
                    MRM.Design.Y.Cell{i}.Label = strjoin(string);

                end
            else

               MRM.Design.Y.Cell{1}.Levels = 1;           
               MRM.Design.Y.Cell{1}.Label = ''; 

            end
        else

            % If there are 0 within-subjects factors all we need is a constant
            MRM.Design.Y.Cell{1}.Levels = 1;
            MRM.Design.Y.Cell{1}.Label = '';

        end
        
end
        


designOK = 1;

end