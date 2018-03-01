function MRM_estimateContrastsPerm
%=========================================================================
% Estimate contrasts of parameter estimates in a multivariate linear model
% using permutation methods
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
% The flag below will cause only one permutation to run using an identity
% matrix. This is useful for checking that the test statistic is being
% calculated correctly during the permutations. This will be saved in a new
% image in the output folder called 'permTestStat.nii' and should be
% identical to the '...refF.nii' image
%--------------------------------------------------------------------------
testPerms = 0;

%--------------------------------------------------------------------------
% Gather Info
%--------------------------------------------------------------------------
M          = MRM.Design.X.X;
outDir     = [MRM.Options.Out filesep];
nDVs       = size(MRM.Design.Y.Cell,2);
cells      = size(MRM.Design.X.Cell,2);
nPerm      = MRM.Options.Thresh.nPerms - 1;
alpha      = MRM.Options.Thresh.PThresh;
checkPerms = MRMoptions.checkPerms;
flips      = MRMoptions.Flips;

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
    %---------
    FImg(1) = deal(struct(...
                        'fname',   [Out 'MRM_' name '_refF.nii'],                ...
                        'dim',     [dimX, dimY, dimZ],                           ...
                        'dt',      [spm_type('float64') spm_platform('bigend')], ...
                        'mat',     T,                                            ...
                        'pinfo',   [1 0 0]',                                     ...
                        'descrip', 'Reference F'));
    
    FImg = spm_data_hdr_write(FImg);
    
    %-P image
    %---------
    PImg(1) = deal(struct(...
                      'fname',   [Out 'MRM_' name '_refP_UC.nii'],     ...
                      'dim',     [dimX, dimY, dimZ],                            ...
                      'dt',      [spm_type('float64') spm_platform('bigend')],  ...
                      'mat',     T,                                             ...
                      'pinfo',   [1 0 0]',                                      ...
                      'descrip', 'Reference P values'));
    
    PImg = spm_data_hdr_write(PImg);
    
    %- 1-P image
    %------------
    MPImg(1) = deal(struct(...
                      'fname',   [Out 'MRM_' name '_1_minus_refP_UC.nii'], ...
                      'dim',     [dimX, dimY, dimZ],                                ...
                      'dt',      [spm_type('float64') spm_platform('bigend')],      ...
                      'mat',     T,                                                 ...
                      'pinfo',   [1 0 0]',                                          ...
                      'descrip', 'Reference 1-P values'));
    
    MPImg = spm_data_hdr_write(MPImg);
    
    %-Permuted P image
    %----------------------------------------------------------------------
    PermPImg(1) = deal(struct(...
                            'fname',   [Out 'MRM_' name '_permutedP_UC.nii'],        ...
                            'dim',     [dimX, dimY, dimZ],                           ...
                            'dt',      [spm_type('float64') spm_platform('bigend')], ...
                            'mat',     T,                                            ...
                            'pinfo',   [1 0 0]',                                     ...
                            'descrip', 'Permutation P values'));
    
    PermPImg = spm_data_hdr_write(PermPImg);
    
    %- 1-Permuted P image
    %----------------------------------------------------------------------
    MPermPImg(1) = deal(struct(...
                            'fname',   [Out 'MRM_' name '_1_minus_permutedP_UC.nii'], ...
                            'dim',     [dimX, dimY, dimZ],                            ...
                            'dt',      [spm_type('float64') spm_platform('bigend')],  ...
                            'mat',     T,                                             ...
                            'pinfo',   [1 0 0]',                                      ...
                            'descrip', 'Permutation 1-P values'));
    
    MPermPImg = spm_data_hdr_write(MPermPImg);
    
    %----------------------------------------------------------------------
    % Images for the thresholding options
    %----------------------------------------------------------------------
    switch MRM.Options.Thresh.Level
        
        case 'Voxel'
            
            switch MRM.Options.Thresh.Method
                
                %----------------------------------------------------------
                % Images for Familywise error correction
                %----------------------------------------------------------
                case 'FWE'
    
                    %-FWE thresholded F image
                    %-------------------------
                    FWEFImg(1) = deal(struct(...
                                            'fname',   [Out 'MRM_' name '_F_FWER_thresholded.nii'],  ...
                                            'dim',     [dimX, dimY, dimZ],                           ...
                                            'dt',      [spm_type('float64') spm_platform('bigend')], ...
                                            'mat',     T,                                            ...
                                            'pinfo',   [1 0 0]',                                     ...
                                            'descrip', 'Familywise error thresholded F'));
                    
                    FWEFImg = spm_data_hdr_write(FWEFImg);

                    %-FWE P image
                    %--------------
                    FWEPImg(1) = deal(struct(...
                                            'fname',   [Out 'MRM_' name '_permutedP_FWER.nii'],      ...
                                            'dim',     [dimX, dimY, dimZ],                           ...
                                            'dt',      [spm_type('float64') spm_platform('bigend')], ...
                                            'mat',     T,                                            ...
                                            'pinfo',   [1 0 0]',                                     ...
                                            'descrip', 'Familywise error adjusted P values'));
                    
                    FWEPImg = spm_data_hdr_write(FWEPImg);

                    %-FWE 1-P image
                    %-----------------------
                    FWEMPImg(1) = deal(struct(...
                                            'fname',   [Out 'MRM_' name '_1_minus_permutedP_FWER.nii'], ...
                                            'dim',     [dimX, dimY, dimZ],                              ...
                                            'dt',      [spm_type('float64') spm_platform('bigend')],    ...
                                            'mat',     T,                                               ...
                                            'pinfo',   [1 0 0]',                                        ...
                                            'descrip', 'Familywise error adjusted 1-P values'));
                    
                    FWEMPImg = spm_data_hdr_write(FWEMPImg);
                 
                %----------------------------------------------------------
                % Images for False Discovery Rate Correction
                %----------------------------------------------------------    
                case 'FDR'
                    
                    %-FDR thresholded F image
                    %-------------------------
                    FDRFImg(1) = deal(struct(...
                                        'fname',   [Out 'MRM_' name '_F_FDRperm_thresholded.nii'], ...
                                        'dim',     [dimX, dimY, dimZ],                             ...
                                        'dt',      [spm_type('float64') spm_platform('bigend')],   ...
                                        'mat',     T,                                              ...
                                        'pinfo',   [1 0 0]',                                       ...
                                        'descrip', 'FDR thresholded F'));
                    
                    FDRFImg = spm_data_hdr_write(FDRFImg);
                    
                    %-FDR Q image
                    %--------------------------------------------------------------
                    QImg(1) = deal(struct(...
                                        'fname',   [Out 'MRM_' name '_FDR_Qvalues_Perm.nii'],    ...
                                        'dim',     [dimX, dimY, dimZ],                           ...
                                        'dt',      [spm_type('float64') spm_platform('bigend')], ...
                                        'mat',     T,                                            ...
                                        'pinfo',   [1 0 0]',                                     ...
                                        'descrip', 'FDR Q values'));

                    QImg = spm_data_hdr_write(QImg);
                    
                    %-FDR 1-Q image
                    %--------------------------------------------------------------
                    MQImg(1) = deal(struct(...
                                        'fname',   [Out 'MRM_' name '_FDR_1_minus_Qvalues_Perm.nii'], ...
                                        'dim',     [dimX, dimY, dimZ],                                ...
                                        'dt',      [spm_type('float64') spm_platform('bigend')],      ...
                                        'mat',     T,                                                 ...
                                        'pinfo',   [1 0 0]',                                          ...
                                        'descrip', 'FDR 1-Q values'));
                    
                    MQImg = spm_data_hdr_write(MQImg);
                    
            end
        
        %------------------------------------------------------------------
        % Images for Uncorrected thresholding
        %------------------------------------------------------------------
        case 'Uncorrected'
            
            ThreshFImg(1) = deal(struct(...
                                     'fname',   [Out 'MRM_' name '_F_Thresholded_UCperm.nii'], ...
                                     'dim',     [dimX, dimY, dimZ],                            ...
                                     'dt',      [spm_type('float64') spm_platform('bigend')],  ...
                                     'mat',     T,                                             ...
                                     'pinfo',   [1 0 0]',                                      ...
                                     'descrip', 'Uncorrected thresholded F'));

            ThreshFImg = spm_data_hdr_write(ThreshFImg);
            
    end
    
    if testPerms == 1
        
        %-Test Perms image
        %----------------------------------------------------------------------
        TestFImg(1) = deal(struct(...
                              'fname',   [outDir 'PermutationContrasts' ...
                                          filesep name filesep          ...
                                          'permTestStat.nii'],          ...
                              'dim',     [dimX, dimY, dimZ],            ...
                              'dt',      [spm_type('float64')           ...
                                          spm_platform('bigend')],      ...
                              'mat',     T, ...
                              'pinfo',   [1 0 0]',...
                              'descrip', ''));

        TestFImg = spm_data_hdr_write(TestFImg);
    end

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
                          'dt',      [spm_type('float64') ...
                                      spm_platform('bigend')],...
                          'mat',     T, ...
                          'pinfo',   [1 0 0]',...
                          'descrip', ''));
    
        MVImg = spm_data_hdr_write(MVImg);
    
    end
    
    
    A    = cell2mat(MRM.Contrasts.Con{i}.A);
    C    = cell2mat(MRM.Contrasts.Con{i}.C);
    cols = size(ones(nSubs, nDVs) * C',2);
    
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
    
    switch MRM.Options.Thresh.Level
        case 'Voxel'
            switch MRM.Options.Thresh.Method
                case 'FWE'
                    fprintf(1,  'Correction: FWE voxel-level\n');
                    fprintf(1, ['Threshold:' '  pFWE < ' num2str(alpha) '\n']);                 
                case 'FDR'
                    fprintf(1,  'Correction: FDR voxel-level\n');
                    fprintf(1, ['Threshold:' '  qFDR < ' num2str(alpha) '\n']);
            end
        case 'Uncorrected'
            fprintf(1,  'Correction: None\n');
            fprintf(1, ['Threshold:' '  p < ' num2str(alpha) '\n']);
        case 'None'
            fprintf(1,  'Correction: None\n');
            fprintf(1, ['Threshold:' '  None\n']);
    end
    fprintf('%s', repmat('-',1,72));
    fprintf(1, '\n');
    
    
    %======================================================================
    %                   REDUCE THE MODEL TO Y* = XB* + E*
    %======================================================================
    
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
            %if isnan(p) ~= 1
            %    df1 = s * (2 * m + s + 1);
            %    df2 = s * (2 * o + s + 1);
            %    LHS = ((2 * o + s + 1)/(2 * m + s + 1));
            %end
            
            
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
            % Wilks' lambda
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
            end
            
            V = NaN(1,1,sliceDim);
            P = NaN(1,1,sliceDim);

            for k = 1:sliceDim
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));

                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0 && ~isnan(p)

                    V(1,1,k) = sum(diag(SSCPH(:,:,k) / SSCPR(:,:,k)));
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
            
            out                         = reshape(F(:,j), dimX, dimY);
            out(isnan(mask(:,:,j)))     = NaN;
            FImg(1).private.dat(:,:,j)  = out;
            
            out                         = reshape(P(1,1,:), dimX, dimY);
            out(isnan(mask(:,:,j)))     = NaN;
            PImg(1).private.dat(:,:,j)  = out;
            
            MPImg(1).private.dat(:,:,j) = 1 - out;
            
        else
            
            FImg(1).private.dat(:,:,j)  = reshape(F(:,j), dimX, dimY);
            PImg(1).private.dat(:,:,j)  = reshape(P(1,1,:), dimX, dimY);
            MPImg(1).private.dat(:,:,j) = 1 - reshape(P(1,1,:), dimX, dimY);
            
        end

        if MRMoptions.SaveMultiStat == 1
            
            if useMask == 1
                out                         = reshape(V(1,1,:), dimX, dimY);
                out(isnan(mask(:,:,j)))     = NaN;
                MVImg(1).private.dat(:,:,j) = out;
            else
                MVImg(1).private.dat(:,:,j) = reshape(V(1,1,:), dimX, dimY);
            end
            
        end
    end
    
    clearvars SSCPH SSCPR nanSSCPH nanSSCPR V P PEs Resids
    
    
    %======================================================================
    %                            PERMUTATIONS
    %====================================================================== 
    
    U                      = zeros(sliceDim, dimZ); % Counter for uncorrected permuted p-value
    U(isnan(F))            = NaN;
    UFWER                  = zeros(sliceDim, dimZ); % Counter for FWER corrected p-values
    UFWER(isnan(F))        = NaN;                   
    analysisMask           = zeros(sliceDim, dimZ);
    analysisMask(isnan(F)) = NaN;
    Fref                   = reshape(FImg(1).private.dat(:,:,:), sliceDim, dimZ);
    Fperm                  = NaN(sliceDim, dimZ);
    PEsPerm                = NaN(rows, cols, sliceDim);
    ResidsPerm             = NaN(nSubs, cols, sliceDim);
    PFmat                  = zeros(nSubs, nSubs, nPerm);
    
    % If there is a mask use it
    if useMask == 1
        
        temp = NaN(dimX, dimY, dimZ);
        
        for j = 1:nSubs
            for k = 1:cols
                temp(:,:,:)       = reshape(ResidZ(j,k,:,:), dimX, dimY, dimZ);
                temp(isnan(mask)) = NaN;
                ResidZ(j,k,:,:)   = reshape(temp(:,:,:), sliceDim, dimZ);
            end
        end
        
        temp(:,:,:)       = reshape(UFWER(:,:), dimX, dimY, dimZ);
        temp(isnan(mask)) = NaN;
        
        UFWER             = reshape(temp(:,:,:), sliceDim, dimZ);
        U                 = reshape(temp(:,:,:), sliceDim, dimZ);
        analysisMask      = reshape(temp(:,:,:), sliceDim, dimZ);
        
        F(isnan(analysisMask)) = NaN;
    end
    
    %----------------------------------------------------------------------
    % Permutation and sign-flipping matrices
    %----------------------------------------------------------------------
    
    fprintf(1, ['Contrast ' num2str(i) ': Creating permutation matrices ...\n']);
    
    if testPerms == 1      
        PFmat(:,:,1) = eye(nSubs);
    else
        for j = 1:nPerm

            Mat = eye(nSubs);
            
            if flips == 1
            
                r   = randi(2,size(Mat,1),1) - 1;

                for k = 1:size(Mat,1)
                    if r(k) == 0
                        Mat(k,k) = -1;
                    end
                end

            end

            PFmat(:,:,j) = Mat(randperm(size(Mat,1)),:); 
        end
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
    
    if testPerms == 1
        nPerm = 1;
    end
    
    % If this is a multivariate contrast initialise the matrices we need
    if cols > 1
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
                
                U(Fperm(:,k) >= Fref(:,k),k) = U(Fperm(:,k) >= Fref(:,k),k) + 1;

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
                                if ~isnan(analysisMask(l,k))
                                    V(:,l) = diag(SSCPH(:,:,l) / SSCPR(:,:,l));
                                end
                            end
                            
                            Fperm(:,k) = LHS * (sum(V) ./ (s - sum(V)));
                            
                        end
                        
                        U(Fperm(:,k) >= Fref(:,k),k) = U(Fperm(:,k) >= Fref(:,k),k) + 1;
                        
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
                                
                                A(:,:,l)          = Con * reshape(PEsPerm(:,l,:), rows, sliceDim);
                                L(:,:,l)          = A(:,:,l)' * CXXC;
                            end
                            
                            V = NaN(1, sliceDim);
                            
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
                        
                        U(Fperm(:,k) >= Fref(:,k),k) = U(Fperm(:,k) >= Fref(:,k),k) + 1;
                        
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
                                    
                                     V(:,l) = diag(SSCPH(:,:,l) / SSCPR(:,:,l));
                                    
                                end
                            end

                            Fperm(:,k) = (df2 * sum(V)) ./ (s^2 * (2 * m + s + 1));
                                
                        end
                        
                        U(Fperm(:,k) >= Fref(:,k),k) = U(Fperm(:,k) >= Fref(:,k),k) + 1;
                        
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
                        
                        U(Fperm(:,k) >= Fref(:,k),k) = U(Fperm(:,k) >= Fref(:,k),k) + 1;
                        
                    end
            end
        end
        
        Fmax  = max(max(Fperm));
        UFWER = UFWER + (Fmax >= F);

         % Save the distribution of the maximum F
        if isempty(Fmax) == 1
            maxDist(j) = 0;
        else
            maxDist(j) = Fmax;
        end
        
        %evalin('base', ['times(' num2str(j) ') = ' num2str(toc) ';']);
        
        time = toc;        
        msg  = sprintf(['Contrast ' num2str(i) ': Permutation ' num2str(j + 1) ...
                        '/' num2str(nPerm + 1) ' in ' num2str(time) ' seconds\n']);
                   
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    %----------------------------------------------------------------------
    % Save the permutation p-values
    %----------------------------------------------------------------------
    U = (U + 1) ./ (nPerm + 1); % The 1 is for the original test
    
    for j = 1:dimZ
        PermPImg(1).private.dat(:,:,j)  = reshape(U(:,j), dimX, dimY);
        MPermPImg(1).private.dat(:,:,j) = 1 - reshape(U(:,j), dimX, dimY);
    end
    
    
    %======================================================================
    %                         CORRECTION OPTIONS
    %======================================================================
    switch MRM.Options.Thresh.Level   
        
        case 'Voxel'
            
            switch MRM.Options.Thresh.Method
                
                case 'FWE'
    
                    %----------------------------------------------------------------------
                    % Family-wise error correction
                    %----------------------------------------------------------------------
                    
                    % Save plot of the null
                    maxDist(j + 1) = max(max(F)); % Add in the original maximum F-value
    
                    if testPerms == 0   
                        MRM_saveNullPlot(maxDist, outDir, name, alpha, max(max(F)), []);
                    end
    
                    % Save the FWER corrected p-values
                    UFWER(UFWER == 0) = 1;
                    UFWER             = UFWER ./ (nPerm + 1);

                    for j = 1:dimZ
                        FWEPImg(1).private.dat(:,:,j)  = reshape(UFWER(:,j), dimX, dimY);
                        FWEMPImg(1).private.dat(:,:,j) = 1 - reshape(UFWER(:,j), dimX, dimY);
                    end

                    % Threshold the F image using the corrected p-values
                    UFWER(UFWER >= alpha) = NaN;
                    Fstar                 = F(:,:);
                    Fstar(isnan(UFWER))   = NaN;

                    % Save the thresholded F image
                    for j = 1:dimZ  
                        FWEFImg(1).private.dat(:,:,j) = reshape(Fstar(:,j), dimX, dimY);  
                    end
                    
                    % Results file
                    fprintf(1, ['Contrast ' num2str(i) ': Writing results table ...\n']);
                    fprintf(1, '\n');
                    MRM_resultsTable(FWEFImg(1).fname, FWEPImg(1).fname, ...
                                    'PermutationContrasts', name, 'perm');
                    
                    clearvars FWEFImg FWEPImg Fstar FWEMPImg
                    
                case 'FDR'
                    
                    %----------------------------------------------------------------------
                    % FDR correction using the permuted p-values
                    %----------------------------------------------------------------------
 
                    %Pvals  = reshape(PermPImg(1).private.dat(:,:,:), sliceDim * dimZ, 1);
                    Pvals  = PermPImg(1).private.dat(:);
                    Qvals  = MRM_qValue(Pvals,MRMoptions.pFDR,i,MRMoptions.BootPiZero,    ...
                                        MRMoptions.BootPiResamps);
                    Qvals = reshape(Qvals, sliceDim, dimZ, 1);
                    QImg(1).private.dat(:,:,:)  = reshape(Qvals, dimX, dimY, dimZ);
                    MQImg(1).private.dat(:,:,:) = 1 - reshape(Qvals, dimX, dimY, dimZ);

                    % Use the q-values to threshold the F image
                    FDRF                   = FImg(1).private.dat(:,:,:);
                    QvalImg                = QImg(1).private.dat(:,:,:);
                    FDRF(QvalImg >= alpha) = NaN;

                    % Save to disc
                    FDRFImg(1).private.dat(:,:,:) = FDRF;
                    
                    % Results file
                    fprintf(1, ['Contrast ' num2str(i) ': Writing results table ...\n']);
                    fprintf(1, '\n');
                    MRM_resultsTable(FDRFImg(1).fname, QImg(1).fname, ...
                                    'PermutationContrasts', name, 'FDR');
                    
                    clearvars Qvals Pvals FDRF QvalImg QImg FDRFImg MQImg    
            end
            
        case 'Uncorrected'
            
            %----------------------------------------------------------------------
            % Uncorrected thresholding using the permuted p-values
            %----------------------------------------------------------------------
            
            % Threshold the permutation p-values
            UncorrectedP                        = PermPImg(1).private.dat(:,:,:);
            UncorrectedP(UncorrectedP >= alpha) = NaN;

            % Use the p-values to threshold the F image
            ThreshF                      = FImg(1).private.dat(:,:,:);
            ThreshF(isnan(UncorrectedP)) = NaN;

            % Save to disc
            ThreshFImg(1).private.dat(:,:,:) = ThreshF;
            
            % Results file
            fprintf(1, ['Contrast ' num2str(i) ': Writing results table ...\n']);
            fprintf(1, '\n');
            MRM_resultsTable(ThreshFImg(1).fname, PermPImg(1).fname, ...
                             'PermutationContrasts', name, 'none');
            
            clearvars UncorrectedP ThreshF ThreshFImg
    end
    
    %----------------------------------------------------------------------
    % Test mode
    %----------------------------------------------------------------------
    if testPerms == 1
        for j = 1:dimZ
            TestFImg(1).private.dat(:,:,j) = reshape(Fperm(:,j), dimX, dimY);
        end
    end

    
    %----------------------------------------------------------------------
    % Tidy up
    %----------------------------------------------------------------------
    clearvars SSCPH SSCPR SSCPHImg SSCPRImg FImg PermPImg PImg MPImg UFWER ...
              Fstar maxDist MPermPImg;

    if MRMoptions.SaveMultiStat == 1        
        clearvars MVImg       
    end
    
    if testPerms == 1
        clearvars TestFImg
    end
    
end

end