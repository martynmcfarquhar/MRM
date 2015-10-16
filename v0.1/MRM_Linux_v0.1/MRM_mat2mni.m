function [mniCoords, value] = MRM_mat2mni(matCoords, imagePath)
%==========================================================================
% Convert matrix indices to MNI (mm) coordinates
%==========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% The mm coordinates can be calculated using the affine transformation
% matrix T in the header using T * [X Y Z 1]'. This is rounded and returned 
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
      
% Load header
hdr = spm_data_hdr_read(imagePath);

% Get the transformation matrix
T = hdr.mat;

% Calculate coordinates
mniCoords = T * [matCoords 1]';
mniCoords(4) = [];
mniCoords = round(mniCoords');

% Read value
value = spm_data_read(hdr, 'xyz', matCoords');
   
% Tidy up
clearvars hdr

end