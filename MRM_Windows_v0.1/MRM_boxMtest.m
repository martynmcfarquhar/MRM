function [M, F, df1, df2, pval] = MRM_boxMtest(rawData, X)
%==========================================================================
% Box's M Test for the equality of covariance matrices
%==========================================================================
% Function description                                                    %
%-------------------------------------------------------------------------%
% Given a multivariate outcome matrix and a cell means coded design matrix
% this function will compute Box's M Tests for the equality of covariance
% matrices. Here X is the design matrix WITH ADDITIONAL COVARIATES REMOVED. 
% The columns of X therefore exclusively code the cells of the design.
%
% See Rencher & Christensen (2012) for calculation details.
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

% Number of populations
m = size(X,2);

% Number of measures
k = size(rawData,2);

% Individual n
n = zeros(m,1);

for i = 1:m   
    n(i) = size(rawData(X(:,i) == 1,:),1);
end

% Individual covariance matrices
%---------------------------------------------------------
for i = 1:m
    Sind{i} = (rawData(X(:,i) == 1,:) - repmat(mean(rawData(X(:,i) == 1,:)),         ...
                                               size(rawData(X(:,i) == 1,:),1),1))' * ...
              (rawData(X(:,i) == 1,:) - repmat(mean(rawData(X(:,i) == 1,:)),         ...
                                               size(rawData(X(:,i) == 1,:),1),1));
    Sind{i} = Sind{i} / (n(i) - 1);
    
    if rcond(Sind{i}) < eps
       msgbox(['Box''s Test cannot be computed because at least one of ' ...
               'your group-specific covariance matrices is singular.']);
           
       M    = NaN;
       F    = NaN;
       df1  = NaN;
       df2  = NaN;
       pval = NaN;
       
       return
    end
end

% Check for singularities
%---------------------------------------------------------

% Pooled covariance matrix
%---------------------------------------------------------
Ssum = zeros(size(Sind{1},1), size(Sind{1},2));

for i = 1:m
    Ssum = Ssum + ((n(i) - 1) * Sind{i});
end

S = Ssum / (sum(n) - m);

% M stat
%---------------------------------------------------------
Msum = 0;

for i = 1:m
    Msum = Msum + ((n(i) - 1) * log(det(Sind{i})));
end

M = (sum(n) - m) * log(det(S)) - Msum;

% c
%---------------------------------------------------------
cSum = 0;

for i = 1:m
   cSum = cSum + (1/(n(i) - 1)); 
end

cSum = cSum - (1/(sum(n) - m));

c = ((2*k^2 + 3*k - 1) / (6 * (k + 1) * (m - 1))) * cSum;

% df1
%---------------------------------------------------------
df1 = (k*(k+1)*(m-1)) / 2;

% c2
%---------------------------------------------------------
c2sum = 0;

for i = 1:m
   c2sum = c2sum + (1/((n(i) - 1)^2));
end

c2sum = c2sum - (1/((sum(n) - m)^2));

c2 = (((k-1)*(k+2)) / (6*(m-1))) * c2sum;

% df2
%---------------------------------------------------------
df2 = (df1 + 2) / abs(c2 - c^2);

% F calculation
%---------------------------------------------------------
F = ((1 - c - df1/df2) / df1) * M;

% P-value
%---------------------------------------------------------
pval = 1 - spm_Fcdf(F, df1, df2);

end