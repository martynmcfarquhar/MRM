function MRM_buildDesignMatrix
%==========================================================================
% Build a cell means design matrix
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This function simply creates a cell means coded design matrix using the
% information in a global MRM object. It also creates the matrix coding the
% factorial structure of Y for use in visualisation.
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

% Design info
cells = size(MRM.Design.X.Cell, 2);
cols  = cell(cells,1);
nSubs = 0;
nDVs  = size(MRM.Design.Y.Cell, 2);
nCVs  = MRM.Covariates.Number;

% Count the number of subjects by counting the files in each cell
for i = 1:cells 
    nSubs = nSubs + size(MRM.Data.Y{1}.Cell{i}.Scans, 2);
end

%--------------------------------------------------------------------------
% Build the design matrix (cell means)
%--------------------------------------------------------------------------
for i = 1:cells
    if i ==1
        n       = size(MRM.Data.Y{1}.Cell{1}.Scans,2);
        cols{i} = [ones(1, n), zeros(1, nSubs - n)];
    else
       n       = size(MRM.Data.Y{1}.Cell{i}.Scans,2);
       temp    = cell2mat(cols(1:i-1));
       nOnes   = sum(temp(:));
       cols{i} = [zeros(1, nOnes), ones(1, n), zeros(1, nSubs - (nOnes + n))];
    end
end

MRM.Design.X.X = cell2mat(cols)';

%--------------------------------------------------------------------------
% Add covariates
%--------------------------------------------------------------------------
if MRM.Covariates.Number ~= 0    
    if isfield(MRM.Covariates, 'CV') == 1
        for i = 1:nCVs
            if MRM.Covariates.CV{i}.Centering == 1
                MRM.Design.X.X           = [MRM.Design.X.X MRM.Covariates.CV{i}.ValueDM];
            else
                MRM.Design.X.X           = [MRM.Design.X.X MRM.Covariates.CV{i}.Value];
            end
            MRM.Covariates.CV{i}.Col = size(MRM.Design.X.X, 2);
        end
    end
end

%--------------------------------------------------------------------------
% Build the factorial matrix for Y (used for plotting only)
%--------------------------------------------------------------------------
MRM.Design.Y.Y = repmat(1:nDVs, nSubs, 1);

end