function [V, F, df1, df2, pval, Qh, Qe, exact] = MRM_stats(Beta, X, A, C, Sigma, stat)
%=========================================================================%
% ESTIMATE THE MULTIVARIATE TEST STATISTIC FOR A HYPOTHESIS OF THE FORM:  %
% ABC' = 0                                                                %
%=========================================================================%
% Function description                                                    %
%-------------------------------------------------------------------------%
% NOTE: This function is not used by MRM, but is used as a reference to the
% calculation of the mutlivariate test statistics outside of an imaging
% context.
%
% This function computes one of four multivariate test statistics 
% associated with a hypothesis test of the form ABC' = 0. Here the matrix B
% is the matrix of parameter estimates from the multivariate linear model 
% Y = XB + E. The matrix A permits the testing of between-subject
% hypotheses and the matrix C permits the testing of within-subject
% hypotheses. Special cases of A and C exist when testing for main effects
% of the between- or within-subject factors. In each case the matrix
% corresponding to tests of no interest contains only 1's. In this sense
% these tests are 'averaged-over' the term of no interest. When testing for
% interactions aross between- and within-subjects factors both matrices
% will contain contrast weights. When more than one between- or within
% subjects factor is present then contrasts will need to be specified with
% care. Both the design-matrix X and the factorial structure in the outcome
% variable Y are presented in cell-means form, and all comparisons should
% be made with this in mind. It would be all too easy to specify a contrast
% that was meaningless.
%
% Whatever the specified contrast all tests are based on two matrices, the
% Sums of Squares and Cross-products (SSCP) matrix for the hypothesis
% (denoted Qh) and the SSCP matrix for the error (denoted Qe). These
% matrices represent the multivariate versions of the numerator and
% denominator in a univariate F-test. All four possible multivariate test
% statistics are based on both the Qh and the Qe matrix.
%
% The following test statistics are available:
%           Pillai's trace     = trace(Qh * inv(Qh + Qe))
%           Wilk's lambda      = det(Qe) / det(Qh + Qe)
%           Hotelling's trace  = trace(Qh * inv(Qe))
%           Roy's largest root = max(eig(inv(Qe) * Qh))
%
% Comparisons between the different test statistics has been conducted and
% each has its own strengths and weaknesses. For each statistic an
% approximation to the univariate F is possible, and under certain
% conditions this approximation is exact and all the test statistics will
% produce the same F value.
%
%-------------------------------------------------------------------------%
% Input arguments                                                         %
%-------------------------------------------------------------------------%
%
% Beta:   The matrix of least-squares parameter estimates
%
% X:      The design-matrix
%
% Sigma:  The estimated variance-covariance matrix from the multivariate
%         model calculated as 1/(n - p) * (e' * e) where n is the number of
%         subjects, p is the number of between-subject parameters (columns 
%         of X), and e is the error matrix.
%
% A:      The left-hand hypothesis matrix in the general hypothesis test 
%         of ABC = 0, where B is the matrix of paprameter estimates. The A 
%         matrix allows for the testing of between-subject factors and is 
%         specified with as many columns as the design matrix. In the case 
%         of no between-subject factors the design matrix contains only an 
%         intercept and A = 1.
%
% C:      The right-hand hypothesis matrix in the general hypothesis test 
%         of ABC' = 0, where B is the matrix of parameter estimates. The C 
%         matrix allows for the testing of within-subject factors and is 
%         specified with as many rows as columns of the outcome matrix Y.
%
% stat:   A string defining which multivariate test statistic to compute. 
%         It should contain one of the following:
%         "PT"  = Pillai's trace
%         "WL"  = Wilk's lambda
%         "HT"  = Hotelling's trace
%         "RLR" = Roy's largest root
%
%-------------------------------------------------------------------------%
% Outputs                                                                 %
%-------------------------------------------------------------------------%
%
% V:    The estimated statistic value for the contrast
%
% F:    The F approximation for the contrast
%
% df1:  The numerator degrees of freedom for the F approximation
%
% df2:  The denominator degrees of freedom for the F approximation
%
% pval: The p-value from the F approximation
%
% Qh:   The sums-of-squares and cross-products matrix for the hypothesis
%
% Qe:   The sums-of-squares and cross-products matrix for the error
%
%-------------------------------------------------------------------------%
% References                                                              %
%-------------------------------------------------------------------------%
%
% (1) Davis, C.S. (2010) Statistical Methods for the Analysis of Repeated 
%     Measurements. London: Springer
% (2) Fox, J. (2008) Applied Regression Analysis and Generalized Linear 
%     Models. London: Sage
% (3) SAS Corp. User's Guide: http://support.sas.com/documentation/
% (4) Crowder, M.J. & Hand, D.J. (1990) Analysis of Repeated Measures. 
%     London: Chapman & Hall
% (5) Christensen, R. (1991). Linear Models for Multivariate, Time Series,
%     and Spatial Data. London: Springer
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
      

n = size(X, 1);     % Number of subjects
h = size(Beta, 1);  % Number of between-subject groups
t = size(Beta, 2);  % Number of within-subject measurements

%-------------------------------------------------------------------------%
% Sums of Squares and Cross Products (SSCP) matrices                      %
%-------------------------------------------------------------------------%

% SSCP matrix for the hypothesis
Qh = (A * Beta * C')' * inv(A * inv(X' * X) * A') * (A * Beta * C');

% SSCP matrix for the residuals
Qe = C * (Sigma * (n - h)) * C';


% Setup info for the test stats
p = rank(Qh + Qe);
q = rank(A * inv(X' * X) * A');
v = n - rank(X);
s = min(p, q);
m = (abs(p - q) - 1) / 2;
o = (v - p - 1) / 2;


% Is this exact?
if rank(Qh) == 1
    exact = 1;
else
    exact = 0;
end

%=========================================================================%
% Test statistics and F approximations 
%=========================================================================%   
% see the SAS documentation for the most comprehensive claculation details 
% 
% http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/
% viewer.htm#statug_introreg_sect012.htm
%
%=========================================================================%
% Pillai's trace                                                          %
%-------------------------------------------------------------------------%

if strcmp(stat, 'PT') == 1

    % Pillai's statistic
    V = trace(Qh * inv(Qh + Qe));

    % F approximation - see SAS documentation
    F = ((2 * o + s + 1)/(2 * m + s + 1)) * (V / (s - V));

    % Degrees of freedom
    df1 = s * (2 * m + s + 1);
    df2 = s * (2 * o + s + 1);

    % P value
    pval = 1 - fcdf(F, df1, df2);

%-------------------------------------------------------------------------%
% Wilk's lambda                                                           %
%-------------------------------------------------------------------------%

elseif strcmp(stat, 'WL') == 1

    % Wilk's lambda (the likelihood ratio statistic)
    V = det(Qe) / det(Qh + Qe);
    
    % F approximation 
    r = v - ((p - q + 1) / 2);
    u = (p * q - 2) / 4;
    
    if p^2 + q^2 - 5 > 0
        
        t = sqrt((p^2 * q^2 - 4) / (p^2 + q^2 - 5));
        
    else
        
        t = 1;
        
    end
    
    F = ((1 - V^(1/t)) / (V^(1/t))) * ((r * t - 2 * u) / (p * q));
    
    % Degrees of freedom
    df1 = p * q;
    df2 = r * t - 2 * u;
    
    % P-value
    pval = 1 - fcdf(F, df1, df2);
    
    % Is this exact? (see Rao, 1973)
    if s <= 2
        exact = 1;
    end

%-------------------------------------------------------------------------%
% Hotelling's trace                                                       %
%-------------------------------------------------------------------------%

elseif strcmp(stat, 'HT') == 1
    
    % Hotelling-Lawley Trace
    V = trace(Qh * inv(Qe));
    
    % F approximation 
    %b = (p + 2 * o) * (q + 2 * o) / (2 * (2 * o + 1) * (o - 1));
    %c = (2 + (p * q + 2) / (b - 1)) / (2 * o);
    
    
    %if o > 0
       
     %   F = (V / c) * ((4 + (p * q + 2) / (b - 1)) / (p * q));
        
        % Degrees of freedom
      %  df1 = p * q;
      %  df2 = 4 + (p * q + 2) / (b - 1);
        
        % P-value
      %  pval = 1 - fcdf(F, df1, df2);
        
   % else
        
        F = (2 * (s * o + 1) * V) / (s^2 * (2 * m + s + 1));
        
        % Degrees of freedom
        df1 = s * (2 * m + s + 1);
        df2 = 2 * (s * o + 1);
        
        % P-value
        pval = 1 - fcdf(F, df1, df2);
        
    %end

%-------------------------------------------------------------------------%
% Roy's largest root                                                      %
%-------------------------------------------------------------------------%

elseif strcmp(stat, 'RLR') == 1

   % Roy's Maximum Root
    V = max(eig(inv(Qe) * Qh));
   
   % F approximation
   r = max(p, q);   
   
   F = V * ((v - r + q) / r);
   
   % Degrees of freedom
   df1 = r;
   df2 = v - r + q;
   
   % P-value
   pval = 1 - fcdf(F, df1, df2);

%-------------------------------------------------------------------------%
% 'stat' can't be anything else so throw an error                         %
%-------------------------------------------------------------------------%

else
 
    error('The argument to ''stat'' is not recognised');
    
end

end