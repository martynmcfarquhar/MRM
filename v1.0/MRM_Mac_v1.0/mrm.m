%=========================================================================%
% Launcher for Multivariate Repeated Measures (MRM)                       %
%=========================================================================%
% Function description                                                    %
%-------------------------------------------------------------------------%
% Not much in here at the moment, but can be used in future to set initial
% values and options for using the software. At present it just checks the
% MATLAB version before calling MRM_launcher
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
function mrm()
% Check the MATLAB version
v = version('-release');
v = v(1,1:numel(v)-1);

if str2double(v) < 2014  
    fprintf(1, ['MRM requires MATLAB 2014b or higher. You are currently running ' version('-release') ...
                '. Please upgrade to the latest version of MATLAB in order to use the software.\n']);
    return
elseif str2double(v) == 2014
    if strcmp(v(1,numel(v)), 'a') % Need 2014b not 2014a
        fprintf(1, ['MRM requires MATLAB 2014b or higher. You are currently running ' version('-release') ...
                '. Please upgrade to the latest version of MATLAB in order to use the software.\n']);
        return
    end
end

% If we've got this far then launch the GUI
MRM_launcher();

end