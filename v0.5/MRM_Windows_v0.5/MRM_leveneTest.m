function [F, df1, df2, pval] = MRM_leveneTest(rawData, X)
%==========================================================================
% Levene's test for homogeneity of variance (center = mean)
%==========================================================================
% Function description                                                    
%--------------------------------------------------------------------------
% Levene's tests is very simple. We take the absolute value of the
% residuals from the full model and compare them in a simple one-way ANOVA
% across the between-subject groups. This is implemented by first computing
% the residuals from the full model (incl. covariates), and then removing
% the covariates from the design-matrix and computing a one-way ANOVA per
% dependent variable, using the absolute values of the residuals in place
% of the raw data.
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


global MRM

beta   = inv(X' * X) * X' * rawData;
resids = rawData - X*beta;

F    = zeros(size(resids,2), 1);
df1  = zeros(size(resids,2), 1);
df2  = zeros(size(resids,2), 1);
pval = zeros(size(resids,2), 1);

X = MRM.Design.X.X(:,1:(size(MRM.Design.X.X,2) - MRM.Covariates.Number));

for i = 1:size(resids,2)
   
    Y    = abs(resids(:,i));
    beta = inv(X' * X) * X' * Y;
    
    con  = diff(eye(size(X,2)));
    
    resid = Y - X*beta;
    
    sigma = sum(resid.^2) / (size(X,1) - size(X,2));
    V     = sigma * inv(X'*X);
    
    F(i) = ((con*beta)'*inv(con*V*con')*(con*beta)) / rank(con*V*con');
    
    df1(i) = rank(con*V*con');
    df2(i) = size(Y,1) - rank(X);
    
    pval(i) = 1 - spm_Fcdf(F(i), df1(i), df2(i));
    
end

end