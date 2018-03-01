function MRM_checkUniquePerms(PFmat, conNum, nSubs, nPerm)
%==========================================================================
% Check for unqiue permutations
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% Quick and dirty (i.e. SLOW) method of checking whether all the
% permutations in a nSubs x nSubs x nPerms PFmat matrix are unique
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

reverseStr = '';
j          = 1;

while j <= nPerm
    
    msg = sprintf(['Contrast ' num2str(conNum) ': Checking all permutation matrices are unique: ' ...
                   num2str(round(j/nPerm*100)) ' percent \n']);
    
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    
    idx = 1:nPerm;
    
    k = 1;
    
    while k <= nSubs
        for l = 1:size(idx,2)
            if isnan(idx(l)) ~= 1
                if PFmat(k,:,idx(l)) == PFmat(k,:,j)
                    if idx(l) ~= j
                        idx(l) = l;
                    else
                        idx(l) = NaN;
                    end
                else
                    idx(l) = NaN;
                end
            end
        end
        
        if sum(isnan(idx)) == size(idx,2)
            k = nSubs;
        end
        
        k = k + 1;
        
    end
    
    j = j + 1;
    
    idx(isnan(idx)) = [];
    
    if isempty(idx) ~= 1
        
        for k = 1:size(idx,2)
            
            fprintf(1, ['Contrast ' num2str(conNum) ': Duplicate found, recalculating and resetting checking\n']);
            
            Mat = eye(nSubs);
            r   = randi(2,size(Mat,1),1) - 1;
            
            for l = 1:size(Mat,1)
                if r(l) == 0
                    Mat(l,l) = -1;
                end
            end
            
            PFmat(:,:,idx(k)) = Mat(randperm(size(Mat,1)),:);
            
        end
        
        j    = 1;
        flag = 1;
        
    else
        flag = 0;
    end
    
    if flag == 1
        reverseStr = '';
        flag       = 0;
        j          = 1;
    end
    
end

end