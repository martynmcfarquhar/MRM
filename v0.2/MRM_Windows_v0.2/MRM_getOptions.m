function ok = MRM_getOptions
%=========================================================================
% Load the global options                 
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% This function loads the options stored in the MRMoptions .mat file found
% in the MRM directory. If no .mat file can be found a new one is created
% using the MRM_defaultOptions() function. If the MRM directory cannot be
% found an error is thrown and the function returns a 0;
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

global MRMoptions
ok = 1;

% Find the MRM directory
if exist('MRM_estimate.m', 'file') == 2
        MRMdir = fileparts(which('MRM_estimate'));
        MRMdir = [MRMdir filesep]; 
    else 
        ok = 0;
        return
end   

% If there is no MRMoptions.mat file then make a new one
if exist([MRMdir 'MRMoptions.mat'], 'file') == 0
    MRM_defaultOptions();
end

% Load the options
load([MRMdir 'MRMoptions.mat']);

end