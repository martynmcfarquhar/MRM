function mrm()
%=========================================================================%
% Launcher for Multivariate Repeated Measures (MRM)                       %
%=========================================================================%
% Function description                                                    %
%-------------------------------------------------------------------------%
% Not much in here at the moment, but can be used in future to set initial
% values and options for using the software. At present it is just a
% wrapper for the MRM_main() function that launches the GUI
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

% Check the MATLAB version
v = version('-release');
v = v(1,1:numel(v)-1);

if str2num(v) < 2013  
    fprintf(1, ['MRM requires MATLAB 2013a or higher. You are currently running ' version('-release') ...
                '. Please upgrade to the latest version of MATLAB in order to use the software.\n']);
    return
end

% Check SPM12 is available
if ~exist('spm.m','file')
    fprintf(1, ['SPM12 is required for MRM to function. Please add SPM12 to your path or download it ' ... 
                'from http://www.fil.ion.ucl.ac.uk/spm/\n']);
    return
else
    if ~strcmp(spm('ver'),'SPM12')
        fprintf(1, ['You appear to be using an older version of SPM. Please download SPM12 ' ...
                    'from http://www.fil.ion.ucl.ac.uk/spm/\n']);
        return
    end
end

% If we've got this far then launch the GUI
MRM_launcher();

end