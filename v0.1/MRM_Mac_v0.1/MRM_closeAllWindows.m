function MRM_closeAllWindows
%==========================================================================
% Close all open windows
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% Convenience function to close all open windows except for the launcher.
% This should make it easier to close windows as new ones are added as
% their definition just needs adding here
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

% Results table
if isempty(findobj('type','figure','name','Results Table')) ~= 1 
    h = findobj('type','figure','name','Results Table');
    delete(h);
end

% Viewer window
if isempty(findobj('type','figure','name','MRMviewer')) ~= 1
    h = findobj('type','figure','name','MRMviewer');
    handles = guidata(h);
    handles.reslicedOverlay = [];
    guidata(h, handles);
    delete(h);
end

% Plotting window
if isempty(findobj('type','figure','name','MRM Bar Plot')) ~= 1
    h = findobj('type','figure','name','MRM Bar Plot');
    delete(h);
end

% Post-estimation tools
if isempty(findobj('type','figure','name','MRM Post-estimation Tools')) ~= 1
    h = findobj('type','figure','name','MRM Post-estimation Tools');
    delete(h);
end

% Repeated measures model
if isempty(findobj('type','figure','name','Repeated measures for fMRI')) ~= 1
    h = findobj('type','figure','name','Repeated measures for fMRI');
    delete(h);
end

% MANOVA model
if isempty(findobj('type','figure','name','MANOVA for fMRI')) ~= 1
    h = findobj('type','figure','name','MANOVA for fMRI');
    delete(h);
end

end