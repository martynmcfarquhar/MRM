function [EvReturn, EigVecRawReturn, UnstandCanonDFreturn, ...
          StandCanonDF, varExpl, discrimVals, centroids] = MRM_LDA(Y, X, L, groupNum, ...
          returnScores, plotScores, Xgroup, groupLabs, coords)
%==========================================================================
% Conduct a descriptive linear discriminants analysis
%==========================================================================
% Function description                                                    
%--------------------------------------------------------------------------
% This script conducts a descriptive linear discriminants analysis with
% inputs as follows
% - Y:         The multivariate outcome matrix used as the predictors
% - X:         The design matrix including any continuous covariates
% - L:         The hypothesis matrix used to construct SSCPH
% - GroupNum:  The number of groups
% - Xgroup:    The design matrix containing only the grouping variables
% - groupLabs: Labels for the groups for plotting
%
% See Rencher & Christensen (2012) and Klecka (1980) for details on how
% the calculations are performed and how the discriminant results are
% interpreted.
%
%-------------------------------------------------------------------------
% References                                                              
%-------------------------------------------------------------------------
% (1) Klecka, W.R., 1980. Discriminant Analysis. Sage, London.
% (2) Rencher, A.C., Christensen, W.F., 2012. Methods of Multivariate 
%     Analysis, 3rd ed. John Wiley & Sons, New York.
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
      
      
%--------------------------------------------------------------------------
% Model estimation
%--------------------------------------------------------------------------
beta  = (X' * X) \ X' * Y;
Resid = Y - X * beta;

% SSCPH and SSCPR
H = (L * beta)' * inv(L * inv(X' * X) * L') * (L * beta);
E = Resid' * Resid;

R = rank(H);

%--------------------------------------------------------------------------
% Eigenvalues of inv(E) * H (in descending order)
%--------------------------------------------------------------------------
Ev       = eig(E \ H);
[Ev, ~]  = sort(Ev, 'descend');

% Only return the first rank(H) eigenvalues
EvReturn = Ev(1:R);

%--------------------------------------------------------------------------
% Cholesky decomposition of E to give the eigenvectors for the
% non-symmetric inv(E) * H
%--------------------------------------------------------------------------

% Because inv(E) * H is not symmetric we use the algorithm for the
% eigenvectors given in Rencher (2012). Briefly, define U = chol(E), then
% the eigenvectors are given by inv(U)*b where b is the corresponding
% eigenvector of inv(U)' * H * inv(U)
U            = chol(E);
[Uvec, Uval] = eig((inv(U))' * H * inv(U));
[~, Uind]    = sort(diag(Uval), 'descend');

% Raw eigenvectors
EigVecRaw       = U \ Uvec;
EigVecRaw       = EigVecRaw(:,Uind);

% Only return the first rank(L) eigenvectors corresponding to the first
% rank(L) eigenvalues
EigVecRawReturn = EigVecRaw(:,1:R);

%--------------------------------------------------------------------------
% Unstandardised discriminant functions
%--------------------------------------------------------------------------

% The unstandardised discriminant functions are the eigenvectors described
% above multiplied by sqrt(N-g) where N is the sample size and g is the
% number of groups
UnstandCanonDF       = EigVecRaw * sqrt(size(X,1) - groupNum);
UnstandCanonDFreturn = UnstandCanonDF(:,1:R);

%--------------------------------------------------------------------------
% Standardised discriminant functions
%--------------------------------------------------------------------------

% The standardised discriminant functions are given by multiplying each
% element of each eigenvector by the square-root of the corresponding
% diagonal element of E
StandCanonDF = zeros(size(Y,2), R);

for i = 1:R    
    StandCanonDF(:,i) = EigVecRaw(:,i) .* sqrt(diag(E));
end

%--------------------------------------------------------------------------
% Variance explained
%--------------------------------------------------------------------------
varExpl = zeros(R,1);

for i = 1:R 
    varExpl(i,1) = Ev(i) / sum(Ev);
end

varExpl = varExpl * 100;

%--------------------------------------------------------------------------
% Constant
%--------------------------------------------------------------------------
cons = zeros(1,R);

for i = 1:R   
    for j = 1:size(UnstandCanonDFreturn,1)  
        cons(1,i) = cons(1,i) + UnstandCanonDFreturn(j,i) * mean(Y(:,j));    
    end 
end

cons = -1 * cons;

UnstandCanonDFreturn = [UnstandCanonDFreturn; cons];

%--------------------------------------------------------------------------
% Discrimination scores
%--------------------------------------------------------------------------
discrimVals = zeros(size(Y,1), R);

Ycons = [Y ones(size(Y,1),1)];

for i = 1:R   
    discrimVals(:,i) = Ycons * UnstandCanonDFreturn(:,i);
end

% Return scores if requested
if returnScores == 1
    assignin('base', 'Ydiscrim', discrimVals);
    fprintf('Discriminant scores returned in variable Ydiscrim\n\n');
end

%--------------------------------------------------------------------------
% Group centroids
%--------------------------------------------------------------------------

% Unstandardised coefficients evaluated at the group means

centroids = zeros(groupNum, R);

for i = 1:R
    
    for j = 1:groupNum
        
        for k = 1:size(UnstandCanonDFreturn,1)-1
            centroids(j,i) = centroids(j,i) + UnstandCanonDFreturn(k,i) * mean(Y(X(:,j) == 1,k));
        end
        
        centroids(j,i) = cons(1,i) + centroids(j,i);
        
    end
end

%--------------------------------------------------------------------------
% Plot scores if requested
%--------------------------------------------------------------------------
if plotScores == 1

    if R == 1
        fprintf('Cannot plot discriminant scores for only one function\n\n');
        return
    end
    
    groupList = zeros(size(Y,1), 1);
    
    for i = 1:size(groupList,1)
        groupList(i,1) = find(Xgroup(i,:));
    end

    ind = 1;

    for i = 1:R  
        for j = 1:R  
            if i < j
                plotnums(ind,1) = i;
                plotnums(ind,2) = j;

                ind = ind + 1;
            end  
        end
    end

    nPlots = size(plotnums,1);
    
    if nPlots == 1
        nCols = 1;
    else
        nCols  = 2;
    end

    % Vague rule to try and keep the plots tidy
    if round(nPlots/nCols) > 4
        nCols = 3;
    end

    nRows = round(nPlots/nCols);

    win = get(0, 'ScreenSize');

    figure('Position',      [win(3)/4, win(4)/4, win(3)/2, win(4)/1.5], ...
        'Units',         'Pixels', ...
        'NumberTitle',   'off',    ...
        'DockControls',  'off', ...
        'ToolBar',       'none');

    for i = 1:size(plotnums,1)
        subplot(nRows, nCols, i);

        gscatter(discrimVals(:,plotnums(i,1)), discrimVals(:,plotnums(i,2)), groupList);
        legend(groupLabs, 'Location', 'best')
        xlabel(['Discrminant function ' num2str(plotnums(i,1))]);
        ylabel(['Discrminant function ' num2str(plotnums(i,2))]);
        
        for j = 1:groupNum
           
            text(centroids(j,plotnums(i,1)), centroids(j, plotnums(i,2)), ...
                 ['+ ' num2str(j)])
            
        end
        
    end
    
    axes('Position',[0 0 1 1], ...
         'Xlim',[0 1], ...
         'Ylim',[0 1], ...
         'Box','off', ...
         'Visible','off', ...
         'Units','normalized', ...
         'clipping' , 'off');
    
    text(0.5, 0.98, ...
        ['Discriminant function scores and group centroids for voxel ' num2str(coords)], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', 'top');
    
end

end