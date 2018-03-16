function MRM_resultsTable(FimgPath, PimgPath, conFolder, conName, correction)
%=========================================================================
% CREATE AND SAVE RESULTS TABLE FOR A PARTICULAR CONTRAST
%=========================================================================
% Function description
%-------------------------------------------------------------------------
% This function creates the results table for an arbitrary contrast, saving
% the table to a text file. To allow for pretty tables the displaytable() 
% function is used, as written and provided by Matt Caywood:
% http://www.mathworks.com/matlabcentral/fileexchange/27920-displaytable
%
% Inputs are as follows
% FimgPath              - The path to the F-image used to fill the table
%                         and define thge clusters. This is usually the
%                         thresholded F-image (either uncorrected, FDR, or
%                         FWE by permutation)
%
% PimgPath              - The path to the P-image to be displayed
%
% conFolder             - The name of the folder where the contrasts are
%                         being saved. Either 'Contrasts' or
%                         'PermutationContrasts'
%
% conName               - The name of the contrast that the table will be
%                         for. It is assumed that this is also the name of
%                         the folder where the contrast results are saved.
%
% correction            - Either 'none', 'FDR', or 'perm' to indicate the
%                         form of correction used.
%
%=========================================================================
% Copyright 2015 Martyn McFarquhar                                        
%=========================================================================
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
%=========================================================================
      
%--------------------------------------------------------------------------
% Load MRM options in case any settings have been changed
%--------------------------------------------------------------------------
global MRMoptions

MRM_getOptions();

global MRM

Out = [MRM.Options.Out filesep conFolder filesep conName filesep];

% Read in images
Fimg = spm_data_read(spm_data_hdr_read(FimgPath));
Pimg = spm_data_read(spm_data_hdr_read(PimgPath));

[dimX, dimY, dimZ] = size(Fimg);

% Turn the Fimg into a list of values
Fimg(Fimg(:) == 0) = NaN;
Fcollapsed = Fimg(:);

% Create the 3 x n matrix of voxel coordinates
indexes = [];

for i = 1:size(Fimg,3)
    [x, y] = find(Fimg(:,:,i));
    indexes = [indexes; x y repmat(i,size(x,1),1)];
end

indexes = indexes';

% Remove indices to NaN's in the F-image
assignin('base','Fimg',Fimg);
assignin('base','indexes',indexes);

indexes(:,isnan(Fcollapsed))  = [];
Fcollapsed(isnan(Fcollapsed)) = [];

% Use the spm_max function to create list of maxima for the table
[N,Z,M,A,xyz] = spm_max(Fcollapsed,indexes);

if MRMoptions.ClustStat == 0 % Cluster mass
                      
    ClustImg = NaN(dimX, dimY, dimZ);
    Clusters = 1:max(A);
    
    for k = 1:size(Clusters,2)
        
        slices = unique(xyz{k}(3,:));
        
        for l = 1:size(slices,2)
            
            slice  = ClustImg(:,:,slices(l));
            coords = xyz{k}(1:2,xyz{k}(3,:) == slices(l));
            linIND = sub2ind(size(slice), coords(1,:)' , coords(2,:)');
            slice(linIND) = Clusters(k);
            ClustImg(:,:,slices(l)) = slice;
            
        end
    end
    
    ClustSizeRef = unique([A N], 'rows');
    
    % Sum of the F values within each cluster
    for k = 1:size(ClustSizeRef,1)
        ClustSizeRef(k,2) = sum(sum(Fimg(ClustImg == ClustSizeRef(k,1))));
        
        N(A == k) = ClustSizeRef(k,2);
        
    end
end


p = zeros(size(A,1),1);
outTable = [A N Z p M'];

%--------------------------------------------------------------------------
% Create table (if there are results)
%--------------------------------------------------------------------------
if isempty(outTable) ~= 1
    
    % Cluster number ascending F-value decending
    outTable = sortrows(outTable, [1 -3]);
    
    Labels = repmat({''},size(outTable,1),1);
    MNI    = zeros(size(outTable,1), 3);
    
    for i = 1:size(outTable,1)
        
        %-MNI coordinates
        %----------------
        [mniCoords, ~] = MRM_mat2mni([outTable(i,5) outTable(i,6) outTable(i,7)], ...
                                      FimgPath);
        MNI(i,:)       = mniCoords;
        
        %-P-values
        %---------------------
        outTable(i,4) = Pimg(outTable(i,5), outTable(i,6), outTable(i,7));
        
        
        %-Atlas labelling
        %---------------------
        
        Labels{i,1} = MRM_getAtlasLabel(MRMoptions.Atlas, MNI(i,:));
        

    end
    
    outTable = [outTable MNI];
    
    switch correction
        case 'FDR'
            PLabel = 'q-value (FDR - voxel)';
        case 'perm'  
            switch MRM.Options.Thresh.Level
                case 'Voxel'
                    PLabel = 'p-value (FWE - voxel)';
                case 'Cluster'
                    switch MRMoptions.ClustStat
                        case 1
                            PLabel = 'p-value (FWE - cluster)';
                        case 0
                            PLabel = 'p-value (FWE - cluster)';
                    end         
            end
        otherwise
            PLabel = 'p-value (Uncorrected - voxel)';
    end
    
    switch MRMoptions.ClustStat
        case 1
            ClustLabel  = 'Cluster size';
        case 0
            ClustLabel  = 'Cluster mass';
    end
    
    % Remove the matrix coordinates
    outTable(:,5:7) = [];
    
    %----------------------------------------------------------------------
    % Write table to file
    %----------------------------------------------------------------------
    fileID = fopen([Out 'MRM_' conName '_Results.txt'], 'w');
    
    if strcmp(MRM.Options.Thresh.Level, 'Cluster') ~= 1
       
        displaytable(outTable, ...
            {'Cluster', ClustLabel, 'F-value', PLabel, 'X (mm)', 'Y (mm)', 'Z (mm)'}, ...
            [7 15 7 numel(PLabel) 6 6 6], ...
            {'d','.3f','.3f','.3f','d','d','d'}, ...
            Labels, ...
            fileID);
    else
        
        outTable(:,3) = []; 
        
        displaytable(outTable, ...
            {'Cluster', ClustLabel, PLabel, 'X (mm)', 'Y (mm)', 'Z (mm)'}, ...
            [7 15 numel(PLabel) 6 6 6], ...
            {'d','.3f','.3f','d','d','d'}, ...
            Labels, ...
            fileID);
    end
    
    fclose(fileID);
    
else
    
    fileID = fopen([Out 'MRM_' conName '_Results.txt'], 'w');
    
    if strcmp(MRM.Options.Thresh.Level, 'Cluster') == 1
        fprintf(fileID, 'No suprathreshold clusters');
    else
        fprintf(fileID, 'No suprathreshold voxels');
    end
    
    fclose(fileID);
    
end

end