function label = MRM_getAtlasLabel(AtlasString, MNIcoords)
%=========================================================================
% Get the region label from the specified atlas at the specified MNI 
% coordinates             
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% Using the atlases provided in the /Utilies folder the string
% corresponding to the region given by the MNI (mm) coordinates is returned
%
% The atlas string can be either 'TD', 'NM', or 'AAL'. The actual value at
% that coordinate is found using the MRM_mni2mat() function.
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

[pathStr, ~, ~] = fileparts(which('MRM_estimate.m'));

if isempty(pathStr)
    fprintf(1, 'MRM folder not found');
    return
else
    AtlasFolders    = [pathStr filesep 'Utilities' filesep];
end

switch AtlasString
    case 'TD'
        %--------------------------------------------------------------
        % TD Atlas
        %--------------------------------------------------------------      
        fid             = fopen([AtlasFolders 'TDAtlas' filesep 'TD_label.csv']);
        TDLabelText     = textscan(fid,'%d%s','delimiter',',');
        fclose(fid);
        
        fid             = fopen([AtlasFolders 'TDAtlas' filesep 'TD_brodmann.csv']);
        TDBrodText      = textscan(fid,'%d%s','delimiter',',');
        fclose(fid);
        
        % Talaraich labels
        [~, value] = MRM_mni2mat([MNIcoords(1) MNIcoords(2) MNIcoords(3)],      ...
                                 [AtlasFolders 'TDAtlas' filesep 'TD_label.nii']);
        
        TalRow = find(TDLabelText{1} == value);
        
        % Brodmann labels
        [~, value] = MRM_mni2mat([MNIcoords(1) MNIcoords(2) MNIcoords(3)],      ...
                                 [AtlasFolders 'TDAtlas' filesep 'TD_brodmann.nii']);
        
        BrodRow = find(TDBrodText{1} == value);
        
        % Side
        switch sign(MNIcoords(1))
            case -1
                side = ' L';
            case 1
                side = ' R';
            otherwise
                side = ' M';
        end
        
        % Label
        if isempty(TalRow)
            label = '';
        else
            label = [TDLabelText{2}{TalRow} ' '];
        end
        
        if isempty(BrodRow) ~= 1
            label = [label TDBrodText{2}{BrodRow}];
        end
        
        label = [label side];
        
        
    case 'NM'
        %--------------------------------------------------------------
        % Neuromorphometrics
        %--------------------------------------------------------------    
        fid             = fopen([AtlasFolders 'NMAtlas' filesep 'NM_label.csv']);
        NMLabelText     = textscan(fid,'%d%s','delimiter',',');
        fclose(fid);
        
        % Neuromorphometrics labels
        [~, value] = MRM_mni2mat([MNIcoords(1) MNIcoords(2) MNIcoords(3)],  ...
                                 [AtlasFolders 'NMAtlas' filesep 'NM_label.nii']);

        NMRow = find(NMLabelText{1} == value);

        % Label
        if isempty(NMRow)
            switch sign(MNIcoords(1))
                case -1
                    side = ' L';
                case 1
                    side = ' R';
                otherwise
                    side = ' M';
            end
            label = side;
        else
            label = [NMLabelText{2}{NMRow}];
        end
        
    case 'AAL'
        %--------------------------------------------------------------
        % AAL
        %--------------------------------------------------------------
        fid             = fopen([AtlasFolders 'AALAtlas' filesep 'AAL_label.csv']);
        AALLabelText    = textscan(fid,'%d%s','delimiter',',');
        fclose(fid);
        
        % AAL labels
        [~, value] = MRM_mni2mat([MNIcoords(1) MNIcoords(2) MNIcoords(3)],     ...
                                 [AtlasFolders 'AALAtlas' filesep 'AAL_label.nii']);

        AALRow = find(AALLabelText{1} == value);

        %-Row label
        %--------------------
        if isempty(AALRow)
            switch sign(MNIcoords(1))
                case -1
                    side = ' L';
                case 1
                    side = ' R';
                otherwise
                    side = ' M';
            end
            label = side;
        else
            label = [AALLabelText{2}{AALRow}];
        end
        
    otherwise
        fprintf(1, 'Atlas string not recognised. Must be either ''TD'', ''NM'', or ''AAL''.\n');
end


end