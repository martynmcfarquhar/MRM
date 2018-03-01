function MRM_estimateContrastsPermClust
%=========================================================================
% Estimate contrasts of parameter estimates in a multivariate linear model
% using permutation methods to calculate CLUSTER statistics
%=========================================================================
% Function description
%-------------------------------------------------------------------------
% This function implements permutation tests of contrasts of the form
% ABC' = 0 in the multivariate linear model Y = MB + E.
%
% The form of permutation testing follows largely from the description
% given in Winkler et al. (2014) for univariate fMRI models. In the
% multivariate case we follow the suggestion of Zeng et al. (2011) and
% first define an alternative model of the form
%
% Y* = MB* + E*
%
% where Y* = YC', B* = BC', and E = EC'. This reduces the multivariate
% hypothesis test ABC' = 0 to AB = 0. The model is then partitioned as
%
% X = M * pinv(L')'
%
% and
%
% Z = M - L' * pinv(L')
%
% to give
%
% Y* = XK + ZG + E*
%
% This partitioning is used to allow for permutation in the prescence of
% nuisance covariates, as explained in Winkler et al. (2014). This allows a
% rexpression of the hypothesis of interest as K = AB, so that testing of
% AB can be perfomed by simply testing K = 0. This is constructed in the
% combined model
%
% Y = [X Z] D = E
%
% using an identity matrix of the same dimension as rows in A padded with
% zeros. E.g. if A was 2 x 6 we would test the first 2 rows of D, where
% D = [K;G]
%
% Permutation of the multivariate residuals after fitting the reduced model
%
% Y* = ZG + E**
%
% follows from the scheme of whole-block exchangeability, as discussed by
% Winkler et al. (2014). This is clear if we consider the univariate
% expression of the multivaiate linear model by re-expressing the columns
% from each subject in turn as rows. We can then see that each row of the
% multivaiate data matrix Y corresponds to a per-subject block in the
% univariate case. Exchangeability is then defined on a per-block basis by
% only permuting the rows of E**, not the columns.
%
% Sign-flipping is performed in an identical fashion to the univariate
% case. Because both the permutation matrix P and the sign-flipping matrix
% S are derived from an identity matrix of size equivalent to the rows of
% E** we only ever permute rows, and only ever changes the sign of an
% entire row at once. This is equivalent to the whole-block exchangeability
% scheme described in Winkler et al. (2014; see Figure 3). Sign-flipping is
% a necessity in the multivariate scheme, otherwise permutation of
% rows-only is only effective for hypotheses on the rows of B (the
% between-subject effects) but not the columns of B (the within-subject
% effects).
%
% In terms of implementation we follow the description of the randomize
% algorithm generously provided in Winkler et al. (2014) for general
% univariate fMRI data.
%
% For cluster-based correction we used the spm_max function to define the
% list of maxima and clusters on the thresholded image. The threshold is
% given as a p-value in the GUI. This p-value is converted into a critical
% F-value and used to threshold the image. On each permutation we simply
% threshold the permuted image, calculate the maximum cluster, and use that
% to feed a counter for the number of time a cluster was found that was
% equal to or bigger than the size of the cluster found in the original
% un-permuated image. This counter image is then divided by the number of
% permuations to produce a p-value image for the clusters. Finally, this is
% used to threshold the original F-image for significant clusters.
%
%-------------------------------------------------------------------------
% References
%-------------------------------------------------------------------------
%
% (1) Davis, C.S. (2010) Statistical Methods for the Analysis of Repeated
%     Measurements. London: Springer
% (2) SAS Corp. User's Guide: http://support.sas.com/documentation/
% (3) Crowder, M.J. & Hand, D.J. (1990) Analysis of Repeated Measures.
%     London: Chapman & Hall
% (4) Christensen, R. (1991) Linear Models for Multivariate, Time Series,
%     and Spatial Data. London: Springer
% (5) Winkler, A.M., Ridgway, G.R., Webster, M.A., Smith, S.M., & Nichols,
%     T.E (2014). Permutation inference for the general linear model.
%     NeuroImage, 92, 381-397.
% (6) Zeng, C., Pan, Z., MaWhinney, S., Baron, A.E., & Zerbe, G.O. (2011).
%     Permuation and F distribution of tests in the multivariate general
%     linear model. The American Statistician, 65(1), 31-36.
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

global MRM

%--------------------------------------------------------------------------
% Load MRM options
%--------------------------------------------------------------------------
global MRMoptions

ok = MRM_getOptions();

if ok == 0
    msgbox('There was an error loading the global options.');
    return
end

%--------------------------------------------------------------------------
% Gather Info
%--------------------------------------------------------------------------
M        = MRM.Design.X.X;
outDir   = [MRM.Options.Out filesep];
nDVs     = size(MRM.Design.Y.Cell,2);
cells    = size(MRM.Design.X.Cell,2);
nPerm    = MRM.Options.Thresh.nPerms - 1;
alpha    = MRM.Options.Thresh.PThresh;
checkPerms = MRMoptions.checkPerms;

% Gather together the files into 1 cell per DV
files = cell(1,nDVs);

for i = 1:nDVs
    for k = 1:cells
        files{i} = [files{i} MRM.Data.Y{i}.Cell{k}.Scans];
    end
end

nSubs = size(files{1},2);

% Read in first image to get info
img      = spm_vol(files{1}{1});
dimX     = img(1).dim(1);
dimY     = img(1).dim(2);
dimZ     = img(1).dim(3);
sliceDim = dimX * dimY;
T        = img.mat;

% Read in all the data
data = cell(nSubs,nDVs,1);

for j = 1:nSubs
    for k = 1:nDVs
        data{j,k} = spm_data_read(spm_data_hdr_read(files{k}{j}));
    end
end

% Make the output directory
mkdir([outDir 'PermutationContrasts']);

% Check for mask
if isempty(MRM.Options.Mask) ~= 1
    mask    = spm_data_read(spm_data_hdr_read([outDir 'MRM_mask.nii']));
    useMask = 1; 
else
    % Hacky fix for the moment
    mask    = spm_data_read(spm_data_hdr_read([outDir 'MRM_Covar_1_1.nii']));
    mask(mask > 0) = 1;
    useMask = 1; 
    %useMask = 0;
end

for i = 1:MRM.Contrasts.Number
    
    %======================================================================
    %                             OUTPUT FILES
    %======================================================================
    
    name      = regexprep(MRM.Contrasts.Con{i}.Name, ' ', '_');
    printName = MRM.Contrasts.Con{i}.Name;
    
    mkdir([outDir 'PermutationContrasts' filesep name]);
    
    Out = [outDir 'PermutationContrasts' filesep name filesep];
    
    %-F image
    %----------
    FImg(1) = deal(struct(...
        'fname',   [Out 'MRM_' name '_refF.nii'],                 ...
        'dim',     [dimX, dimY, dimZ],                            ...
        'dt',      [spm_type('float32') spm_platform('bigend')],  ...
        'mat',     T,                                             ...
        'pinfo',   [1 0 0]',                                      ...
        'descrip', 'Reference F'));
    
    FImg = spm_data_hdr_write(FImg);
    
    %-Thresholded F image
    %---------------------
    ThreshFImg(1) = deal(struct(...
        'fname',   [Out 'MRM_' name '_F_FWER_cluster_thresholded.nii'], ...
        'dim',     [dimX, dimY, dimZ],                                  ...
        'dt',      [spm_type('float32') spm_platform('bigend')],        ...
        'mat',     T,                                                   ...
        'pinfo',   [1 0 0]',                                            ...
        'descrip', 'F values thresholded by cluster P values'));
    
    ThreshFImg = spm_data_hdr_write(ThreshFImg);
    
    %-P image
    %---------
    ClustPImg(1) = deal(struct(...
        'fname',   [Out 'MRM_' name '_cluster_P.nii'],           ...
        'dim',     [dimX, dimY, dimZ],                           ...
        'dt',      [spm_type('float32') spm_platform('bigend')], ...
        'mat',     T,                                            ...
        'pinfo',   [1 0 0]',                                     ...
        'descrip', 'Cluster P values'));
    
    ClustPImg = spm_data_hdr_write(ClustPImg);
    
    
    %- 1-P image
    %----------------------------------------------------------------------
    MPImg(1) = deal(struct(...
        'fname',   [outDir filesep 'PermutationContrasts' ...
        filesep name filesep 'MRM_' name ...
        '_1_minus_cluster_P.nii'], ...
        'dim',     [dimX, dimY, dimZ], ...
        'dt',      [spm_type('float32') ...
        spm_platform('bigend')],...
        'mat',     T, ...
        'pinfo',   [1 0 0]',...
        'descrip', ''));
    
    MPImg = spm_data_hdr_write(MPImg);
    
    
    %- Cluster label image
    %----------------------------------------------------------------------
    ClustImg(1) = deal(struct(...
        'fname',   [outDir filesep 'PermutationContrasts' ...
        filesep name filesep 'MRM_' name ...
        '_cluster_labels.nii'], ...
        'dim',     [dimX, dimY, dimZ], ...
        'dt',      [spm_type('float32') ...
        spm_platform('bigend')],...
        'mat',     T, ...
        'pinfo',   [1 0 0]',...
        'descrip', 'Cluster numbers'));
    
    ClustImg = spm_data_hdr_write(ClustImg);
    
    
    %-Multivariate stat image (if requested)
    %----------------------------------------------------------------------
    if MRMoptions.SaveMultiStat == 1
        
        statName = MRM.Options.Stat.Name;
        
        MVImg(1) = deal(struct(...
            'fname',   [outDir filesep 'PermutationContrasts' ...
            filesep name filesep 'MRM_' name ...
            '_ref' statName ...
            '.nii'], ...
            'dim',     [dimX, dimY, dimZ], ...
            'dt',      [spm_type('float32') ...
            spm_platform('bigend')],...
            'mat',     T, ...
            'pinfo',   [1 0 0]',...
            'descrip', ''));
        
        MVImg = spm_data_hdr_write(MVImg);
        
    end
    
    A    = cell2mat(MRM.Contrasts.Con{i}.A);
    C    = cell2mat(MRM.Contrasts.Con{i}.C);
    cols = size(zeros(nSubs, nDVs) * C',2);
    
    %======================================================================
    %                             PRINT INFO
    %======================================================================
    fprintf(1, '\n');
    fprintf('%s', repmat('=',1,72));
    fprintf(1, '\n');
    fprintf(1, ['Contrast ' num2str(i) '\n']);
    fprintf('%s', repmat('=',1,72));
    fprintf(1, '\n');
    fprintf(1, ['Name:' '       ' printName '\n']);
    fprintf(1, ['P-values:' '   Permutation\n']);
    
    if cols > 1
        switch MRM.Options.Stat.Name
            case 'PT'
                fprintf(1, ['Statistic:' '  Pillai''s trace\n']);
            case 'WL'
                fprintf(1, ['Statistic:' '  Wilks'' lambda\n']);
            case 'HT'
                fprintf(1, ['Statistic:' '  Hotelling-Lawley trace\n']);
            case 'RLR'
                fprintf(1, ['Statistic:' '  Roy''s largest root\n']);
        end
    end
    
    if MRMoptions.ClustStat == 1
        fprintf(1,  'Correction: FWE cluster-level (cluster size)\n');
    else
        fprintf(1,  'Correction: FWE cluster-level (cluster mass)\n');
    end
    fprintf(1, ['Threshold:' '  pFWE < ' num2str(alpha) '\n']);
    fprintf('%s', repmat('-',1,72));
    fprintf(1, '\n');
    
    
    %======================================================================
    %                   REDUCE THE MODEL TO Y* = XB* + E*
    %======================================================================
    
    %     A = cell2mat(MRM.Contrasts.Con{i}.A);
    %     C = cell2mat(MRM.Contrasts.Con{i}.C);
    
    Y     = zeros(nSubs, nDVs, sliceDim);
    Ystar = NaN(nSubs, size(zeros(nSubs, nDVs) * C',2), sliceDim, dimZ);
    
    fprintf(1, ['Contrast ' num2str(i) ': Reducing model ...\n']);
    
    % Y* = Y * C'
    %----------------------------------------------------------------------
    for j = 1:dimZ
        
        for k = 1:nSubs
            for l = 1:nDVs
                slice = data{k,l}(:,:,j);
                Y(k,l,:) = reshape(slice, 1, 1, sliceDim);
            end
        end
        
        for k = 1:sliceDim
            if sum(sum(isnan(Y(:,:,k)))) ~= nSubs * nDVs
                Ystar(:,:,k,j) = Y(:,:,k) * C';
            end
        end
    end
    
    %======================================================================
    %                      REPARAMETERISE THE MODEL
    %======================================================================
    
    % Model partition after Ridgway (2009) - see Winkler et al. (2014)
    X     = M * pinv(A);
    Z     = M - (M * A' * pinv(A'));
    [Z,~] = svd(Z);
    Z     = Z(:,1:rank(M) - rank(A'));
    Rz    = eye(size(Z,1)) - Z * pinv(Z);    % (I - ZZ*)
    X     = Rz * X;                          % Orthogonalise X wrt Z
    
    % New contrast for testing K = 0 in Y* = XK + ZG + E*, specified in
    % terms of the combined model Y* = [X Z]D = E*
    Con  = padarray(eye(size(A,1)), [0 (size(X,2) + size(Z,2)) - size(A,1)], 0, 'post');
    
    pXX  = pinv([X Z]' * [X Z]);            % (M'M)*
    CXXC = pinv(Con * pXX * Con');          % (C (M'M)* C')*
    Rz   = eye(size(Z,1)) - Z * pinv(Z);    % (I - ZZ*)
    pM   = pinv([X Z]);                     %  M*
    Ry   = eye(size([X Z],1)) - [X Z] * pM; % (I - MM*)
    % * denotes pseudoinverse
    
    %======================================================================
    %                  REFERENCE F, P, V, AND 1-P IMAGES
    %======================================================================
    rows = size(X,2) + size(Z,2);
    
    ResidZ = NaN(nSubs, cols, sliceDim, dimZ);
    PEs    = NaN(rows, cols, sliceDim);
    F      = NaN(sliceDim, dimZ);
    Resids = NaN(nSubs, cols, sliceDim);
    SSCPH  = NaN(cols, cols, sliceDim);
    SSCPR  = NaN(cols, cols, sliceDim);
    
    fprintf(1, ['Contrast ' num2str(i) ': Calculating reference images as permutation 1 ...\n']);
    
    for j = 1:dimZ
        
        for k = 1:cols
            
            ResidZ(:,k,:,j) = Rz * squeeze(Ystar(:,k,:,j));
            PEs(:,k,:)      = pM * squeeze(ResidZ(:,k,:,j));
            Resids(:,k,:)   = Ry * squeeze(ResidZ(:,k,:,j));
            
        end
        
        for k = 1:sliceDim
            
            SSCPH(:,:,k) = (Con * PEs(:,:,k))' * CXXC * (Con * PEs(:,:,k));
            SSCPR(:,:,k) = Resids(:,:,k)' * Resids(:,:,k);
            
        end
        
        %------------------------------------------------------------------
        % Find the first voxel in this slice where SSCPH and SSCPR are 
        % not NaN. Use this to calculate p = rank(SSCPH + SSCPR). Check for
        % SSCPH and SSCPR where spme values are essentially 0 as well. This
        % will help prevent scenarios in the multivariate comparisons where
        % there is a voxel in one modality with values that are empty in the
        % other
        %------------------------------------------------------------------
        for k = 1:sliceDim
            
            nanSSCPH = isnan(SSCPH(:,:,k));
            nanSSCPR = isnan(SSCPR(:,:,k));
            
            if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0
                
                p = rank(SSCPH(:,:,k) + SSCPR(:,:,k));
                
                if min(min(abs(SSCPH(:,:,k)))) > eps && min(min(abs(SSCPR(:,:,k)))) > eps && p >= cols
                    break
                else
                    p = NaN;
                end
                
            else
                p = NaN;
            end
        end
        
        % Calculate the parameters we need if there is data in this slice
        if ~isnan(p)
            q   = rank(Con * pinv([X Z]' * [X Z]) * Con');
            v   = nSubs - rank([X Z]);
            s   = min(p,q);
            m   = (abs(p - q) - 1) / 2;
            o   = (v - p - 1) / 2;
            df1 = s * (2 * m + s + 1);
            df2 = s * (2 * o + s + 1);
            LHS = ((2 * o + s + 1)/(2 * m + s + 1));
        end
        
        if strcmp(MRM.Options.Stat.Name, 'PT') == 1
            
            %--------------------------------------------------------------
            % Pillai's trace
            %--------------------------------------------------------------
            if isnan(p) ~= 1
                ClustP = MRM.Options.Thresh.ClustThresh;
                critF  = spm_invFcdf(1 - ClustP, [df1 df2]);
            end
            
            V = NaN(1,1,sliceDim);
            P = NaN(1,1,sliceDim);
            
            for k = 1:sliceDim
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0 && ~isnan(p)
                    
                    V(1,1,k) = sum(diag(SSCPH(:,:,k) / (SSCPH(:,:,k) + SSCPR(:,:,k))));
                    F(k,j)   = LHS * (V(1,1,k) / (s - V(1,1,k)));
                    
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:)     = 1 - spm_Fcdf(F(:,j), [df1 df2]);
                P(P(:) == 1) = NaN;
            end
            
        elseif strcmp(MRM.Options.Stat.Name, 'WL') == 1
            
            %--------------------------------------------------------------
            % Wilk's lambda
            %--------------------------------------------------------------
            if isnan(p) ~= 1
                r = v - ((p - q + 1) / 2);
                u = (p * q - 2) / 4;
                
                if p^2 + q^2 - 5 > 0
                    t = sqrt((p^2 * q^2 - 4) / (p^2 + q^2 - 5));
                else
                    t = 1;
                end
                
                df1 = p * q;
                df2 = r * t - 2 * u;
                RHS = df2/df1;
                
                ClustP = MRM.Options.Thresh.ClustThresh;
                critF  = spm_invFcdf(1 - ClustP, [df1 df2]);
            end
            
            V = NaN(1,1,sliceDim);
            P = NaN(1,1,sliceDim);
            
            for k = 1:sliceDim
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0 && ~isnan(p)
                    
                    V(1,1,k) = det(SSCPR(:,:,k)) / det(SSCPH(:,:,k) + SSCPR(:,:,k));
                    F(k,j)   = ((1 - V(1,1,k)^(1/t)) / (V(1,1,k)^(1/t))) * RHS;
                    
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:)     = 1 - spm_Fcdf(F(:,j), [df1 df2]);
                P(P(:) == 1) = NaN;
            end
            
        elseif strcmp(MRM.Options.Stat.Name, 'HT') == 1
            
            %--------------------------------------------------------------
            % Hotelling-Lawley trace
            %--------------------------------------------------------------
            
            if isnan(p) ~= 1
                
                df1 = s * (2 * m + s + 1);
                df2 = 2 * (s * o + 1);
                
                ClustP = MRM.Options.Thresh.ClustThresh;
                critF  = spm_invFcdf(1 - ClustP, [df1 df2]);
            end
            
            V = NaN(1,1,sliceDim);
            P = NaN(1,1,sliceDim);
            
            for k = 1:sliceDim
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0 && ~isnan(p)
                    
                    V(1,1,k) = sum(diag(SSCPR(:,:,k) \ (SSCPH(:,:,k))));
                    
                    F(k,j) = (df2 * V(1,1,k)) / (s^2 * (2 * m + s + 1));
                    
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:) = 1 - spm_Fcdf(F(:,j), [df1 df2]);
                P(P(:) == 1) = NaN;
            end
            
        elseif strcmp(MRM.Options.Stat.Name, 'RLR') == 1
            
            %--------------------------------------------------------------
            % Roy's largest root
            %--------------------------------------------------------------
            if isnan(p) ~= 1
                z   = max(p, q);
                df1 = z;
                df2 = v - z + q;
            end
            
            V = NaN(1,1,sliceDim);
            P = NaN(1,1,sliceDim);
            
            for k = 1:sliceDim
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0 && ~isnan(p)
                    
                    V(1,1,k) = max(eig(SSCPR(:,:,k) \ SSCPH(:,:,k)));
                    F(k,j) = V(1,1,k) * (v - z + q) / z;
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:) = 1 - spm_Fcdf(F(:,j), [df1 df2]);
                P(P(:) == 1) = NaN;
            end
            
        end
        
        %------------------------------------------------------------------
        % Save slice
        %------------------------------------------------------------------
        if useMask == 1
            
            out = reshape(F(:,j), dimX, dimY);
            out(isnan(mask(:,:,j))) = NaN;
            
            if isnan(p) ~= 1
                out(out <= critF) = NaN; % Cluster threshold
            end
            
            FImg(1).private.dat(:,:,j) = out;
            
        else
            
            out = reshape(F(:,j), dimX, dimY);
            
            if isnan(p) ~= 1
                out(out <= critF) = NaN; % Cluster threshold
            end
            
            FImg(1).private.dat(:,:,j) = out;
            
        end
        
        if MRMoptions.SaveMultiStat == 1
            
            if useMask == 1
                out = reshape(V(1,1,:), dimX, dimY);
                out(isnan(mask(:,:,j))) = NaN;
                
                if isnan(p) ~= 1
                    out(out <= critF) = NaN; % Cluster threshold
                end
                
                MVImg(1).private.dat(:,:,j) = out;
            else
                
                out = reshape(V(1,1,:), dimX, dimY);
                
                if isnan(p) ~= 1
                    out(out <= critF) = NaN; % Cluster threshold
                end
                
                MVImg(1).private.dat(:,:,j) = out;
            end
            
        end
    end
    
    %----------------------------------------------------------------------
    % Calculate cluster stats
    %----------------------------------------------------------------------
    
    % Turn the Fimg into a list of values
    Fcollapsed = FImg(1).private.dat(:);
    
    % Create the 3 x n matrix of voxel coordinates
    indexes = [];
    
    for mm = 1:size(FImg(1).private.dat,3)
        [x, y] = find(FImg(1).private.dat(:,:,mm));
        indexes = [indexes; x y repmat(mm,size(x,1),1)];
    end
    
    indexes = indexes';
    
    % Remove indices to NaN's in the F-image
    indexes(:,isnan(Fcollapsed))  = [];
    Fcollapsed(isnan(Fcollapsed)) = [];
    
    % Use the spm_max function to create list of maxima for the table
    [N,~,~,A,xyz] = spm_max(Fcollapsed,indexes);
    
    %----------------------------------------------------------------------
    % Create the image of cluster numbers
    %----------------------------------------------------------------------
    clustImg = NaN(dimX, dimY, dimZ);
    clusters = 1:max(A);
    
    for k = 1:size(clusters,2)
        
        slices = unique(xyz{k}(3,:));
        
        for l = 1:size(slices,2)
            
            slice  = clustImg(:,:,slices(l));
            coords = xyz{k}(1:2,xyz{k}(3,:) == slices(l));
            linIND = sub2ind(size(slice), coords(1,:)' , coords(2,:)');
            slice(linIND) = clusters(k);
            clustImg(:,:,slices(l)) = slice;
            
        end
    end
    
    ClustImg(1).private.dat(:,:,:) = clustImg(:,:,:);
    
    ClustSizeRef = unique([A N], 'rows');
    
    % Cluster mass
    if MRMoptions.ClustStat == 0
        
        Fstats = FImg(1).private.dat(:,:,:);
        
        for k = 1:size(ClustSizeRef,1)
            ClustSizeRef(k,2) = sum(sum(Fstats(clustImg == ClustSizeRef(k,1))));
        end
        
    end
    
    
    % Skip the permutations if there are no clusters found after
    % thresholding
    if isempty(ClustSizeRef) == 1
        fprintf(1, ['Contrast ' num2str(i) ': No clusters found - skipping contrast\n \n']);
        
        % delete the results folder
        rmdir([outDir 'PermutationContrasts' filesep name], 's');
        
        clearvars SSCPH SSCPR SSCPHImg SSCPRImg FImg ClustPImg MPImg ...
            UFWER ThreshFImg FWERPImg Fstar FWERMPImg ClustImg;
        
        continue
    end
    
    % Counter for the cluster sizes
    ClustCounter = zeros(dimX, dimY, dimZ);
    ClustCounter(isnan(FImg(1).private.dat(:,:,:))) = NaN;
    
    clearvars SSCPH SSCPR nanSSCPH nanSSCPR V P PEs Resids
    
    
    %======================================================================
    %                            PERMUTATIONS
    %======================================================================
    
    Fperm           = NaN(sliceDim, dimZ);
    PEsPerm         = NaN(rows, cols, sliceDim);
    ResidsPerm      = NaN(nSubs, cols, sliceDim);
    PFmat           = zeros(nSubs, nSubs, nPerm);
    
    analysisMask           = zeros(sliceDim, dimZ);
    analysisMask(isnan(F)) = NaN;
    
    
    % If there is a mask use it to mask ResidZ and F
    if useMask == 1
        
        analysisMask(isnan(reshape(mask,sliceDim, dimZ))) = NaN;
        
        temp = NaN(dimX, dimY, dimZ);
        
        for j = 1:nSubs
            for k = 1:cols
                temp(:,:,:) = reshape(ResidZ(j,k,:,:), dimX, dimY, dimZ);
                temp(isnan(mask)) = NaN;
                ResidZ(j,k,:,:) = reshape(temp(:,:,:), sliceDim, dimZ);
            end
        end
    end
    
    %----------------------------------------------------------------------
    % Permutation and sign-flipping matrices
    %----------------------------------------------------------------------
    
    fprintf(1, ['Contrast ' num2str(i) ': Calculating permutation matrices ...\n']);
    
    for j = 1:nPerm
        
        Mat = eye(nSubs);
        
        if strcmp(MRM.Model, 'Repeated') == 1
            
            r   = randi(2,size(Mat,1),1) - 1;
            
            for k = 1:size(Mat,1)
                if r(k) == 0
                    Mat(k,k) = -1;
                end
            end
            
        end
        
        PFmat(:,:,j) = Mat(randperm(size(Mat,1)),:);
    end
    
    %----------------------------------------------------------------------
    % If requested check the permutation matrices for uniqueness
    %----------------------------------------------------------------------
    if checkPerms == 1
        
        MRM_checkUniquePerms(PFmat, i, nSubs, nPerm)
        
    end
    
    %----------------------------------------------------------------------
    % Start permutations
    %----------------------------------------------------------------------
    reverseStr = '';
    
    maxDist = zeros(nPerm + 1, 1);
    
    % If this is a multivariate contrast initialise the matrices we need
    if cols ~= 1
        SSCPH    = NaN(cols, cols, sliceDim);
        SSCPR    = NaN(cols, cols, sliceDim);
        A        = NaN(size(Con,1), sliceDim, cols);
        L        = NaN(sliceDim, size(Con,1), cols);
    end
    
    
    for j = 1:nPerm
        
        tic
        
        %------------------------------------------------------------------
        % Univariate contrast
        %------------------------------------------------------------------
        if cols == 1
            
            for k = 1:dimZ
                
                if sum(sum(isnan(F(:,k)))) ~= sliceDim
                    
                    residZresh        = reshape(ResidZ(:,1,:,k), nSubs, sliceDim);
                    PEsPerm(:,1,:)    = pM * (PFmat(:,:,j) * residZresh);
                    ResidsPerm(:,1,:) = Ry * (PFmat(:,:,j) * residZresh);
                    
                    A           = Con * reshape(PEsPerm, rows, sliceDim);
                    L           = A' * CXXC;
                    L           = L' .* A;
                    
                    if size(L,1) > 1
                        L = sum(L);
                    end
                    
                    residsReshape = reshape(ResidsPerm,nSubs,sliceDim);
                    sig           = residsReshape .* residsReshape;
                    sig           = sum(sig);
                    V             = L ./ (L + sig);
                    Fperm(:,k)    = (LHS .* (V ./ (s - V)))';
                    
                end
                
            end
            
            %------------------------------------------------------------------
            % Multivariate contrast
            %------------------------------------------------------------------
        else
            
            switch MRM.Options.Stat.Name
                
                %----------------------------------------------------------
                % Pillai's trace
                %----------------------------------------------------------
                case 'PT'
                    
                    for k = 1:dimZ
                        
                        if sum(sum(isnan(analysisMask(:,k)))) ~= sliceDim
                            
                            for l = 1:cols
                                residZresh        = reshape(ResidZ(:,l,:,k), nSubs, sliceDim);
                                PEsPerm(:,l,:)    = pM * (PFmat(:,:,j) * residZresh);
                                ResidsPerm(:,l,:) = Ry * (PFmat(:,:,j) * residZresh);
                            end
                            
                            V = NaN(cols, sliceDim);
                            
                            for l = 1:cols
                                A(:,:,l) = Con * reshape(PEsPerm(:,l,:), rows, sliceDim);
                                L(:,:,l) = A(:,:,l)' * CXXC;
                            end
                            
                            for l = 1:cols
                                for x = 1:cols
                                    
                                    temp = L(:,:,l)' .* A(:,:,x);
                                    
                                    if size(temp,1) > 1
                                        SSCPH(l,x,:)    = sum(temp,1);
                                    else
                                        SSCPH(l,x,:)    = temp;
                                    end
                                end
                            end
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    SSCPR(:,:,l) = ResidsPerm(:,:,l)' * ResidsPerm(:,:,l);
                                    
                                end
                            end
                            
                            SSCPR = SSCPH + SSCPR;
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    V(:,l) = diag(SSCPH(:,:,l) / SSCPR(:,:,l));
                                    
                                end
                            end
                            
                            Fperm(:,k) = LHS * (sum(V) ./ (s - sum(V)));
                            
                        end
                        
                    end
                    
                    %--------------------------------------------------------------
                    % Wilks' lambda
                    %--------------------------------------------------------------
                case 'WL'
                    
                    for k = 1:dimZ
                        
                        if sum(sum(isnan(analysisMask(:,k)))) ~= sliceDim
                            
                            for l = 1:cols
                                residZresh        = reshape(ResidZ(:,l,:,k), nSubs, sliceDim);
                                PEsPerm(:,l,:)    = pM * (PFmat(:,:,j) * residZresh);
                                ResidsPerm(:,l,:) = Ry * (PFmat(:,:,j) * residZresh);
                            end
                            
                            V = NaN(1, sliceDim);
                            
                            for l = 1:cols
                                A(:,:,l) = Con * reshape(PEsPerm(:,l,:), rows, sliceDim);
                                L(:,:,l) = A(:,:,l)' * CXXC;
                            end
                            
                            for l = 1:cols
                                for x = 1:cols
                                    
                                    temp = L(:,:,l)' .* A(:,:,x);
                                    
                                    if size(temp,1) > 1
                                        SSCPH(l,x,:)    = sum(temp,1);
                                    else
                                        SSCPH(l,x,:)    = temp;
                                    end
                                end
                            end
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    SSCPR(:,:,l) = ResidsPerm(:,:,l)' * ResidsPerm(:,:,l);
                                    
                                end
                            end
                            
                            SSCPH = SSCPH + SSCPR;
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    V(1,l) = det(SSCPR(:,:,l)) / det(SSCPH(:,:,l));
                                    
                                end
                            end
                            
                            Fperm(:,k) = ((1 - V.^(1/t)) ./ (V.^(1/t))) * RHS;
                            
                        end
                        
                    end
                    
                    
                    %--------------------------------------------------------------
                    % Hotelling-Lawley trace
                    %--------------------------------------------------------------
                case 'HT'
                    
                    for k = 1:dimZ
                        
                        if sum(sum(isnan(analysisMask(:,k)))) ~= sliceDim
                            
                            for l = 1:cols
                                residZresh        = reshape(ResidZ(:,l,:,k), nSubs, sliceDim);
                                PEsPerm(:,l,:)    = pM * (PFmat(:,:,j) * residZresh);
                                ResidsPerm(:,l,:) = Ry * (PFmat(:,:,j) * residZresh);
                            end
                            
                            V = NaN(cols, sliceDim);
                            
                            for l = 1:cols
                                A(:,:,l) = Con * reshape(PEsPerm(:,l,:), rows, sliceDim);
                                L(:,:,l) = A(:,:,l)' * CXXC;
                            end
                            
                            for l = 1:cols
                                for x = 1:cols
                                    
                                    temp = L(:,:,l)' .* A(:,:,x);
                                    
                                    if size(temp,1) > 1
                                        SSCPH(l,x,:)    = sum(temp,1);
                                    else
                                        SSCPH(l,x,:)    = temp;
                                    end
                                end
                            end
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    SSCPR(:,:,l) = ResidsPerm(:,:,l)' * ResidsPerm(:,:,l);
                                    
                                end
                            end
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    V(:,l) = diag(SSCPR(:,:,l) \ SSCPH(:,:,l));
                                    
                                end
                            end
                            
                            Fperm(:,k) = (df2 * sum(V)) ./ (s^2 * (2 * m + s + 1));
                            
                        end
                    end
                    
                    %--------------------------------------------------------------
                    % Roy's largest root
                    %--------------------------------------------------------------
                case 'RLR'
                    
                    for k = 1:dimZ
                        
                        if sum(sum(isnan(analysisMask(:,k)))) ~= sliceDim
                            
                            for l = 1:cols
                                residZresh        = reshape(ResidZ(:,l,:,k), nSubs, sliceDim);
                                PEsPerm(:,l,:)    = pM * (PFmat(:,:,j) * residZresh);
                                ResidsPerm(:,l,:) = Ry * (PFmat(:,:,j) * residZresh);
                            end
                            
                            V = NaN(1, sliceDim);
                            
                            for l = 1:cols
                                A(:,:,l) = Con * reshape(PEsPerm(:,l,:), rows, sliceDim);
                                L(:,:,l) = A(:,:,l)' * CXXC;
                            end
                            
                            for l = 1:cols
                                for x = 1:cols
                                    
                                    temp = L(:,:,l)' .* A(:,:,x);
                                    
                                    if size(temp,1) > 1
                                        SSCPH(l,x,:)    = sum(temp,1);
                                    else
                                        SSCPH(l,x,:)    = temp;
                                    end
                                end
                            end
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    SSCPR(:,:,l) = ResidsPerm(:,:,l)' * ResidsPerm(:,:,l);
                                    
                                end
                            end
                            
                            
                            for l = 1:sliceDim
                                if isnan(analysisMask(l,k)) ~= 1
                                    
                                    V(:,l) = max(eig(SSCPR(:,:,l) \ SSCPH(:,:,l)));
                                    
                                end
                            end
                            
                            Fperm(:,k) = V .* ((v - z + q) / z);
                            
                        end
                    end
            end
        end
            
            %------------------------------------------------------------------
            % Calculate cluster stats
            %------------------------------------------------------------------
            
            % Turn the threholded Fimg into a list of values
            FpermReshaped = reshape(Fperm, dimX, dimY, dimZ);
            FpermReshaped(FpermReshaped <= critF) = NaN; % Cluster thresh
            Fcollapsed = FpermReshaped(:);
            
            % Create the 3 x n matrix of voxel coordinates
            indexes = [];
            
            for mm = 1:dimZ
                [x, y] = find(FpermReshaped(:,:,mm));
                indexes = [indexes; x y repmat(mm,size(x,1),1)];
            end
            
            indexes = indexes';
            
            % Remove indices to NaN's in the F-image
            indexes(:,isnan(Fcollapsed))  = [];
            Fcollapsed(isnan(Fcollapsed)) = [];
            
            %------------------------------------------------------------------
            % Cluster mass
            %------------------------------------------------------------------
            if MRMoptions.ClustStat == 0
                
                [N,~,~,count,xyz] = spm_max(Fcollapsed,indexes);
                
                PermClustImg = NaN(dimX, dimY, dimZ);
                PermClusters = 1:max(count);
                
                for k = 1:size(PermClusters,2)
                    
                    slices = unique(xyz{k}(3,:));
                    
                    for l = 1:size(slices,2)
                        
                        slice  = PermClustImg(:,:,slices(l));
                        coords = xyz{k}(1:2,xyz{k}(3,:) == slices(l));
                        linIND = sub2ind(size(slice), coords(1,:)' , coords(2,:)');
                        slice(linIND) = PermClusters(k);
                        PermClustImg(:,:,slices(l)) = slice;
                        
                    end
                end
                
                PermClustSizeRef = unique([count N], 'rows');
                
                % Sum of the F values within each cluster
                for k = 1:size(PermClustSizeRef,1)
                    PermClustSizeRef(k,2) = sum(sum(FpermReshaped(PermClustImg == PermClustSizeRef(k,1))));
                end
                
                if isempty(PermClustSizeRef) == 1
                    maxClust = [];
                else
                    maxClust = max(PermClustSizeRef(:,2));
                end
                
                %------------------------------------------------------------------
                % Cluster size
                %------------------------------------------------------------------
            else
                [N,~,~,~,~] = spm_max(Fcollapsed,indexes);
                maxClust = max(N);
            end
            
            for mm = 1:size(ClustSizeRef,1)
                
                if maxClust >= ClustSizeRef(mm,2)
                    
                    ClustCounter(clustImg == ClustSizeRef(mm,1)) =        ...
                        ClustCounter(clustImg == ClustSizeRef(mm,1)) + 1;
                end
            end
            
            % Save the distribution of the maximum cluster (either size or mass)
            if isempty(maxClust) == 1
                maxDist(j) = 0;
            else
                maxDist(j) = maxClust;
            end
            
            time = toc;
            
            msg = sprintf(['Contrast ' num2str(i) ': Permutation ' num2str(j + 1) ...
                '/' num2str(nPerm + 1) ' in ' num2str(time) ' seconds\n']);
            
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    %----------------------------------------------------------------------
    % Save plot of the null
    %----------------------------------------------------------------------
    maxDist(j + 1) = max(ClustSizeRef(:,2)); % Add in the original maximum cluster size
    
    MRM_saveNullPlot(maxDist, outDir, name, alpha, max(ClustSizeRef(:,2)), MRMoptions.ClustStat);
    
    
    %----------------------------------------------------------------------
    % Save FWER cluster p-values and thresholded F image
    %----------------------------------------------------------------------
    ClustPImg(1).private.dat(:,:,:) = (ClustCounter(:,:,:) + 1) ./ (nPerm + 1);
    MPImg(1).private.dat(:,:,:)     = 1 - ((ClustCounter(:,:,:) + 1) ./ (nPerm + 1));
    
    ThreshF = FImg(1).private.dat(:,:,:);
    ClustP  = ClustPImg(1).private.dat(:,:,:);
    ThreshF(ClustP >= alpha) = NaN;
    
    ThreshFImg(1).private.dat(:,:,:) = ThreshF;
    
    
    %----------------------------------------------------------------------
    % Results file
    %----------------------------------------------------------------------
    fprintf(1, ['Contrast ' num2str(i) ': Writing results table ...\n']);
    fprintf(1, '\n');
    
    MRM_resultsTable(ThreshFImg(1).fname, ClustPImg(1).fname, ...
        'PermutationContrasts', name, 'perm');
    
    
    %----------------------------------------------------------------------
    % Tidy up
    %----------------------------------------------------------------------
    clearvars SSCPH SSCPR SSCPHImg SSCPRImg FImg ClustPImg MPImg ...
        UFWER ThreshFImg FWERPImg Fstar FWERMPImg ClustImg;
    
    
end