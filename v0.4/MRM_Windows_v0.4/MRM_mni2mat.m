function [matCoords, value] = MRM_mni2mat(mniCoords, imagePath)
%==========================================================================
% Convert MNI (mm) coordinates to matrix indices
%==========================================================================
% Converting MNI (mm) coordinates to location in the image matrix can be
% achieved using the affine transformation matrix T in the header using
% [X Y Z 1] * inv(T)'
%
%=========================================================================
% Author:      Martyn McFarquhar                                          
% Email:       martyn.mcfarquhar@.manchester.ac.uk                        
% Date:        28/10/2014                                                 
% Institution: The University of Manchester                               
% Version:     0.1                                                        
%=========================================================================


% Load image header
hdr = spm_data_hdr_read(imagePath);

% Load transformation matrix
T = hdr.mat;

% Calculate the MNI coordinates
matCoords = [mniCoords 1] * (inv(T))';
matCoords(4) = [];
matCoords = round(matCoords);

% Make sure we don't have any negative indices (if we do turn them into 1)
for i = 1:3
    if sign(matCoords(i)) == -1
        matCoords(i) = 1;
    end
end

% Get the image value
try
    value = spm_data_read(hdr, 'xyz', matCoords');
catch
    value = NaN;
end
             
clearvars hdr

end