function MRM_defaultOptions
%==========================================================================
% (Re)create the MRM default options .mat file
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% This is used to create the default options .mat file, or to re-create it
% in the event of resetting the options
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

% Find the MRM directory
%-----------------------

if exist('MRM_estimate.m', 'file') == 2
    MRMdir = fileparts(which('MRM_estimate'));
    MRMdir = [MRMdir filesep]; 
else 
    msgbox('Can''t find the MRM directory. Resetting options aborted.');
    return
end

f = filesep;

% Make the options file
%----------------------
MRMoptions.UncorrThresh  = 0.001;
MRMoptions.FDRqThresh    = 0.05;
MRMoptions.FWEpThresh    = 0.05;
MRMoptions.BootPiZero    = 0;
MRMoptions.BootPiResamps = 100;
MRMoptions.PiZeroOne     = 1;
MRMoptions.pFDR          = 0;
MRMoptions.nPerms        = 5000;
MRMoptions.ClusterP      = 0.001;
MRMoptions.checkPerms    = 0;
MRMoptions.ClustStat     = 1;
MRMoptions.Flips         = 1;

MRMoptions.SaveSSCP      = 0;
MRMoptions.SaveResid     = 0;
MRMoptions.SaveMultiStat = 0;

MRMoptions.Template      = [MRMdir 'Utilities' f 'Template' f 'template1.nii'];
MRMoptions.Atlas         = 'TD';

% Save the options file
%----------------------
save([MRMdir 'MRMoptions.mat'], 'MRMoptions');

end