function qValsFull = MRM_qValue(pVals, robust, conNum, doBoot, nBoots)
%==========================================================================
% Calculate the q-values for FDR correction of a list of p-values
%==========================================================================
% Function description                                                    
%--------------------------------------------------------------------------
% Based on the bootstrap procedure for estimating pi_0 described in Storey
% (2002), and as implemented in the q-value software. The typical FDR would
% would be equivalent to using pi_0 = 1.
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
      

%==========================================================================
% Make sure p-values are arranged as a list in a single column
%==========================================================================
if size(pVals,1) == 1  
    pVals = pVals'; 
end

% Final vector to be returned (qvals + NaNs)
qValsFull = NaN(size(pVals,1),1);

% Get the indexes of the non-NaN values
NotNaNind = find(isnan(pVals) == 0);

% Throw out NaNs
pVals(isnan(pVals)) = [];

% Sort the p-values in ascending order
[pVals, index] = sort(pVals);

numP   = size(pVals,1);

%==========================================================================
% Bootstrap estimate of pi0
%==========================================================================
if doBoot == 1
    
    %nBoots      = 100;
    lambdaRange = (0:0.05:0.95)';
    pi0         = zeros(size(lambdaRange,1),1);

    for i = 1:size(lambdaRange,1)  
        pi0(i) = size(pVals(pVals > lambdaRange(i)),1) / ...
                 ((1 - lambdaRange(i)) * size(pVals,1));
    end

    pi0min    = min(pi0);            
    MSE       = zeros(size(lambdaRange,1),1);
    bootPi0   = zeros(size(lambdaRange,1),1);

    reverseStr = '';

    %----------------------------------------------------------------------
    % Resamples
    %----------------------------------------------------------------------
    for i = 1:nBoots

        bootPvals = randsample(pVals,size(pVals,1),true);

        for j = 1:size(lambdaRange,1)

            bootPi0(j) = size(bootPvals(bootPvals > lambdaRange(j)),1) / ...
                         ((1 - lambdaRange(j)) * size(bootPvals,1));

        end

        MSE = MSE + (bootPi0 - pi0min).^2;   

        %------------------------------------------------------------------
        % Print progress
        %------------------------------------------------------------------
        msg = sprintf(['Contrast ' num2str(conNum) ': FDR - Bootstrapping estimate of pi_0: ' ...
                        num2str(round(i/nBoots*100)) ' percent\n']);

        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end

    pi0 = min(pi0(MSE == min(MSE)));
    pi0 = min(pi0,1);
    
    fprintf(1, ['Contrast ' num2str(conNum) ': FDR - pi_0 = ' ...
                num2str(pi0) '\n']);
    
else 
    pi0 = 1;
end

%==========================================================================           
% Calculate q-values
%==========================================================================
   
order = 1:numP;

if robust == 1
    qVals = (pi0 .* numP .* pVals) ./ (order' .* (1 - (1 - pVals).^numP)); 
else
    qVals = (pi0 .* numP .* pVals) ./ order';
end

reverseStr = '';

for i = numP-1:-1:1
         
    qVals(i,1) = min(qVals(i,1), qVals(i+1,1));
    
    prog = num2str(round((numP-i)/(numP-1)*100));

    %----------------------------------------------------------------------
    % Print progress
    %----------------------------------------------------------------------
    msg = sprintf(['Contrast ' num2str(conNum) ': FDR - Calculating q-values: ' ...
                    prog ' percent\n']);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

%==========================================================================
% Return in the original order
%==========================================================================
unsorted = 1:size(pVals,1);
newInd(index) = unsorted;
qVals = qVals(newInd);

%--------------------------------------------------------------------------
% Put the vector back together by inserting the qVals back between the
% original NaN values
%--------------------------------------------------------------------------
qValsFull(NotNaNind) = qVals;

end