%=========================================================================
% Plot the variance-covariance structure                 
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
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
function MRM_drawVcov(MRM, V, coords)

nSubs = size(MRM.Design.X.X,1);
win   = get(0, 'ScreenSize');  
visV  = kron(eye(nSubs), V);
visV(visV == 0) = NaN;

figure('Position',      [1, 1, win(3)/2.5, win(4)/1.5], ...
       'Units',         'Pixels', ...
       'Color',         [1 1 1],  ...
       'Name',          'Vcov',   ...
       'DockControls',  'off',    ...
       'NumberTitle',   'off',    ...
       'Resize',        'on');   

imagesc(visV);
caxis([(min(V(:)) - 0.1) max(V(:))]);
set(gca,'xtick', []);
set(gca,'xticklabel', []);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
cmap = gray;
cmap(1,:) = [0.85 0.85 0.85];
colormap(cmap);

if size(MRM.Design.Y.Y,2) > 1
    title(['Estimated covariance structure of Vec$(\mathbf{Y})$ for voxel ' num2str(coords(1)) ' ' ...
          num2str(coords(2)) ' ' num2str(coords(3))], 'interpreter', 'latex');
else
    title(['Estimated covariance structure of $\mathbf{Y}$ for voxel ' num2str(coords(1)) ' ' ...
          num2str(coords(2)) ' ' num2str(coords(3))], 'interpreter', 'latex');
end
end