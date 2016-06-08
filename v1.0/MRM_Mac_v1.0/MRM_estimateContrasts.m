function MRM_estimateContrasts
%=========================================================================
% Estimate contrasts of parameter estimates in a multivariate linear model                
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% The general linear hypothesis testing scheme in a multivariate linear
% model is ABC' = 0, where A is the contrast matrix for the between-subject
% structure and C is the contrast matrix for the within-subject structure.
% B is the matrix of estimated parameters. For cases where size(B,2) = 1
% (i.e. a univariate model) C = 1 and the scheme is identical to the
% univariate case.
%
% Four different test statistics are available, with their
% F-approximations, and include
%   - Pillai's trace
%   - Wilks' lambda
%   - Hotelling-Lawley trace
%   - Roy's largest root
% See any of the references for details on the calculation of these values.
%
% Thresholded images are saved using either an uncorrected p-value scheme,
% or the False Discovery Rate (FDR) correction.
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
% (5) Rencher, A., and Christensen, (2014). Methods of Multivariate
%     Analysis.
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
global MRMoptions

% Gather the parameter estimates together
%----------------------------------------
rows  = size(MRM.Design.X.X, 2);
cols  = size(MRM.Data.Y, 2);

img = spm_vol([MRM.Options.Out filesep 'MRM_PE_1_1.nii']);

dimX = img.dim(1);
dimY = img.dim(2);
dimZ = img.dim(3);
T    = img.mat;

nVox = dimX * dimY; % Voxels in a slice

PEs = NaN(rows, cols, dimX, dimY, dimZ, 'double');

for i = 1:rows  
    for j = 1:cols 
        PEs(i,j,:,:,:) = spm_data_read(spm_data_hdr_read([MRM.Options.Out filesep 'MRM_PE_' num2str(i) '_' num2str(j) '.nii']));
    end
end

% Gather the variance-covariance matrix together
%-----------------------------------------------
Vcov = NaN(cols, cols, dimX, dimY, dimZ, 'double');

for i = 1:cols
    for j = 1:cols
        Vcov(i,j,:,:,:) = spm_data_read(spm_data_hdr_read([MRM.Options.Out filesep 'MRM_Covar_' num2str(i) '_' num2str(j) '.nii'])); 
    end
end

% Unscale
X    = MRM.Design.X.X;
Vcov = Vcov .* (size(X, 1) - size(X, 2));


% Calculate the inverse to save doing it for every contrast
%----------------------------------------------------------
XX   = inv(X' * X);
n    = size(X,1);


% Correction?
%----------------
switch MRM.Options.Thresh.Level  
    case 'None'
        doCorrection = 0;
    case 'Uncorrected'   
        doCorrection = 1; 
    case 'Voxel'  
        doCorrection = 1;
end

alpha = MRM.Options.Thresh.PThresh;


% Make folder
%------------
mkdir([MRM.Options.Out filesep 'Contrasts']);

%==========================================================================
%                           CONTRASTS LOOPS
%==========================================================================

for i = 1:MRM.Contrasts.Number
    
    % Contrast name
    name      = regexprep(MRM.Contrasts.Con{i}.Name, ' ', '_');
    printName = MRM.Contrasts.Con{i}.Name;
    
    % Contrast specific sub-folder
    mkdir([MRM.Options.Out filesep 'Contrasts' filesep name]);
    Out = [MRM.Options.Out filesep 'Contrasts' filesep name filesep];

    %======================================================================
    % OUTPUT FILES
    %======================================================================
    CRows = size(cell2mat(MRM.Contrasts.Con{i}.C),1);
    
    if MRMoptions.SaveSSCP == 1
        
        % SSCP folder
        %-------------------
        mkdir([Out 'SSCP']);
        
        %-SSCPH
        %----------------------------------------------------------------------
        SSCPHImg(1:CRows^2) = deal(struct(...
                                          'fname',   [], ...
                                          'dim',     [dimX, dimY, dimZ], ...
                                          'dt',      [spm_type('float64') ...
                                                      spm_platform('bigend')],...
                                          'mat',     T, ...
                                          'pinfo',   [1 0 0]',...
                                          'descrip', ''));

        ind = 1;
        for j = 1:CRows
            for k = 1:CRows 
               SSCPHImg(ind).fname    = [Out 'SSCP' filesep 'MRM_' name ...
                                         '_SSCPH_' num2str(j) '_' ...
                                         num2str(k) '.nii'];
               ind = ind + 1;
            end
        end

        SSCPHImg = spm_data_hdr_write(SSCPHImg);

        %-SSCPR
        %----------------------------------------------------------------------
        SSCPRImg(1:CRows^2) = deal(struct(...
                                          'fname',   [], ...
                                          'dim',     [dimX, dimY, dimZ], ...
                                          'dt',      [spm_type('float64') ...
                                                      spm_platform('bigend')],...
                                          'mat',     T, ...
                                          'pinfo',   [1 0 0]',...
                                          'descrip', ''));

        ind = 1;
        for j = 1:CRows
            for k = 1:CRows 
               SSCPRImg(ind).fname    = [Out 'SSCP' filesep 'MRM_' name ...
                                         '_SSCPR_' num2str(j) '_' ...
                                         num2str(k) '.nii'];
               ind = ind + 1;
            end
        end

        SSCPRImg = spm_data_hdr_write(SSCPRImg);
    end
    
    %-F image
    %----------------------------------------------------------------------
    FImg(1) = deal(struct(...
                          'fname',   [Out 'MRM_' name '_F.nii'], ...
                          'dim',     [dimX, dimY, dimZ], ...
                          'dt',      [spm_type('float64') ...
                                      spm_platform('bigend')],...
                          'mat',     T, ...
                          'pinfo',   [1 0 0]',...
                          'descrip', ''));
    
    MRM.Contrasts.Con{i}.File = [Out 'MRM_' name '_F.nii'];
    FImg = spm_data_hdr_write(FImg);
    
    switch MRM.Options.Thresh.Level
        
        %-FDR F image
        %------------------------------------------------------------------
        case 'Voxel'
        
            FDRFImg(1) = deal(struct(...
                                     'fname',   [Out 'MRM_' name '_F_Thresholded_FDR.nii'], ...
                                     'dim',     [dimX, dimY, dimZ], ...
                                     'dt',      [spm_type('float64') ...
                                                 spm_platform('bigend')],...
                                     'mat',     T, ...
                                     'pinfo',   [1 0 0]',...
                                     'descrip', ''));

            MRM.Contrasts.Con{i}.FileThresh = [Out 'MRM_' name '_F_Thresholded_FDR.nii'];
            FDRFImg = spm_data_hdr_write(FDRFImg);

            %-FDR Q image
            %--------------------------------------------------------------
            QImg(1) = deal(struct(...
                                      'fname',   [Out 'MRM_' name '_FDR_Qvalues.nii'], ...
                                      'dim',     [dimX, dimY, dimZ], ...
                                      'dt',      [spm_type('float64') ...
                                                  spm_platform('bigend')],...
                                      'mat',     T, ...
                                      'pinfo',   [1 0 0]',...
                                      'descrip', ''));

             QImg = spm_data_hdr_write(QImg);
             
             %-FDR 1-Q image
             %--------------------------------------------------------------
             MQImg(1) = deal(struct(...
                                      'fname',   [Out 'MRM_' name '_FDR_1_minus_Qvalues.nii'], ...
                                      'dim',     [dimX, dimY, dimZ], ...
                                      'dt',      [spm_type('float64') ...
                                                  spm_platform('bigend')],...
                                      'mat',     T, ...
                                      'pinfo',   [1 0 0]',...
                                      'descrip', ''));

             MQImg = spm_data_hdr_write(MQImg);
        
        %-Thresholded F image
        %------------------------------------------------------------------
        case 'Uncorrected'
        
            ThreshFImg(1) = deal(struct(...
                                     'fname',   [Out 'MRM_' name '_F_Thresholded_UC.nii'], ...
                                     'dim',     [dimX, dimY, dimZ], ...
                                     'dt',      [spm_type('float64') ...
                                                 spm_platform('bigend')],...
                                     'mat',     T, ...
                                     'pinfo',   [1 0 0]',...
                                     'descrip', ''));

            MRM.Contrasts.Con{i}.FileThresh = [Out 'MRM_' name '_F_Thresholded_UC.nii'];
            ThreshFImg = spm_data_hdr_write(ThreshFImg);
        
    end
    
    
    %-P image
    %----------------------------------------------------------------------
    PImg(1) = deal(struct(...
                          'fname',   [Out 'MRM_' name '_P_UC_.nii'], ...
                          'dim',     [dimX, dimY, dimZ], ...
                          'dt',      [spm_type('float64') ...
                                      spm_platform('bigend')],...
                          'mat',     T, ...
                          'pinfo',   [1 0 0]',...
                          'descrip', ''));
    
    PImg = spm_data_hdr_write(PImg);
    
    %- 1-P image
    %----------------------------------------------------------------------
    MPImg(1) = deal(struct(...
                          'fname',   [Out 'MRM_' name '_1_minus_P_UC.nii'], ...
                          'dim',     [dimX, dimY, dimZ], ...
                          'dt',      [spm_type('float64') ...
                                      spm_platform('bigend')],...
                          'mat',     T, ...
                          'pinfo',   [1 0 0]',...
                          'descrip', ''));
    
    MPImg = spm_data_hdr_write(MPImg);
    
    %-Multivariate stat image (if requested)
    %----------------------------------------------------------------------
    if MRMoptions.SaveMultiStat  == 1
        
        statName = MRM.Options.Stat.Name;
        
        MVImg(1) = deal(struct(...
                          'fname',   [Out 'MRM_' name '_' statName '.nii'], ...
                          'dim',     [dimX, dimY, dimZ], ...
                          'dt',      [spm_type('float32') ...
                                      spm_platform('bigend')],...
                          'mat',     T, ...
                          'pinfo',   [1 0 0]',...
                          'descrip', ''));
    
        MVImg = spm_data_hdr_write(MVImg);
    
    end

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
    fprintf(1, ['P-values:' '   Approximate\n']);
    
    if CRows > 1
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
            fprintf(1,  'Correction: FDR voxel-level\n');
            fprintf(1, ['Threshold:' '  qFDR < ' num2str(alpha) '\n']);
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
    %                             CONTRASTS
    %======================================================================   
    A = cell2mat(MRM.Contrasts.Con{i}.A);
    C = cell2mat(MRM.Contrasts.Con{i}.C);
    
    AXXA = inv(A * XX * A');

    SSCPH = NaN(CRows, CRows, nVox, 'double');
    SSCPR = NaN(CRows, CRows, nVox, 'double');
    
    reverseStr = '';

    for j = 1:dimZ
        
        PEsReshaped  = reshape(PEs(:,:,:,:,j),  [size(PEs, 1), size(PEs, 2), nVox]);
        VcovReshaped = reshape(Vcov(:,:,:,:,j), [size(Vcov, 1), size(Vcov, 2), nVox]);
        
        
        %==================================================================
        % SSCPH and SSCPR for this slice
        %==================================================================
        for k = 1:nVox
            
            SSCPH(:,:,k) = (A * PEsReshaped(:,:,k) * C')' * AXXA * (A * PEsReshaped(:,:,k) * C');
            SSCPR(:,:,k) = C * VcovReshaped(:,:,k) * C';
            
        end
        
        % Find the first voxel in this slice for which SSCPH and SSCPR are
        % not NaN. Use this to calculate p = rank(SSCPH + SSCPR)
        for k = 1:nVox
            
            nanSSCPH = isnan(SSCPH(:,:,k));
            nanSSCPR = isnan(SSCPR(:,:,k));
            
            if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0 
                
                p = rank(SSCPH(:,:,k) + SSCPR(:,:,k));
                
                if p > 0
                    
                    if sum(sum(abs(SSCPH(:,:,k)) <= eps)) == 0 && sum(sum(abs(SSCPR(:,:,k)) <= eps)) == 0
                        break
                    else
                        p = NaN;
                    end
                    
                else
                    p = NaN;
                end
            else 
                p = NaN;
            end
        end
        
        
        % Save slice if requested
        if MRMoptions.SaveSSCP == 1
            ind = 1;
            for k = 1:CRows
                for l = 1:CRows
                    
                    SSCPHImg(ind).private.dat(:,:,j) = reshape(SSCPH(k,l,:), dimX, dimY);
                    SSCPRImg(ind).private.dat(:,:,j) = reshape(SSCPR(k,l,:), dimX, dimY);
                    
                    ind = ind + 1;
                    
                end
            end
        end
        
        
        %==================================================================
        % Multivariate stat, F approximation, and p-values
        %==================================================================
        % See the SAS documentation:
        %
        % http://support.sas.com/documentation/cdl/en/statug/63033/HTML/
        % default/viewer.htm#statug_introreg_sect012.htm
        %------------------------------------------------------------------
        
        % Calculate the parameters we need if there is data in this slice
        if isnan(p) ~= 1
            q = rank(A * inv(X' * X) * A');
            v = n - rank(X);
            s = min(p,q);
            m = (abs(p - q) - 1) / 2;
            o = (v - p - 1) / 2;
        end
        
        if strcmp(MRM.Options.Stat.Name, 'PT') == 1
            
            %--------------------------------------------------------------
            % Pillai's trace
            %--------------------------------------------------------------
            
            if isnan(p) ~= 1
                df1 = s * (2 * m + s + 1);
                df2 = s * (2 * o + s + 1);
                LHS = ((2 * o + s + 1)/(2 * m + s + 1));
            end
            
            V = NaN(1,1,nVox);
            F = NaN(1,1,nVox);
            P = NaN(1,1,nVox);
            
            for k = 1:nVox
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0   

                        V(1,1,k) = sum(diag(SSCPH(:,:,k) / (SSCPH(:,:,k) + SSCPR(:,:,k))));
                        F(1,1,k) = LHS * (V(1,k) / (s - V(1,k)));

                end
                
            end
            
            if isnan(p) ~= 1
                P(1,1,:) = 1 - spm_Fcdf(F(1,1,:), [df1 df2]);
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
            
            V = NaN(1,1,nVox);
            F = NaN(1,1,nVox);
            P = NaN(1,1,nVox);
            
            for k = 1:nVox
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 && sum(nanSSCPR(:)) == 0
                        V(1,1,k) = det(SSCPR(:,:,k)) / det(SSCPH(:,:,k) + SSCPR(:,:,k));
                        F(1,1,k) = ((1 - V(1,1,k)^(1/t)) / (V(1,1,k)^(1/t))) * RHS;
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:) = 1 - spm_Fcdf(F(1,1,:), [df1 df2]);
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
            
            V = NaN(1,1,nVox);
            F = NaN(1,1,nVox);
            P = NaN(1,1,nVox);
            
            for k = 1:nVox
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 || sum(nanSSCPR(:)) == 0
                    
                    V(1,1,k) = sum(diag(SSCPR(:,:,k) \ SSCPH(:,:,k)));
                    F(1,1,k) = (df2 * V(1,1,k)) / (s^2 * (2 * m + s + 1));
                    
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:) = 1 - spm_Fcdf(F(1,1,:), [df1 df2]);
                P(P(:) == 1) = NaN;
            end
            
            
        elseif strcmp(MRM.Options.Stat.Name, 'RLR') == 1
            
            %--------------------------------------------------------------
            % Roy's largest root
            %--------------------------------------------------------------
            
            if isnan(p) ~= 1
                r   = max(p, q);
                df1 = r;
                df2 = v - r + q;
            end
            
            V = NaN(1,1,nVox);
            F = NaN(1,1,nVox);
            P = NaN(1,1,nVox);
            
            for k = 1:nVox
                
                nanSSCPH = isnan(SSCPH(:,:,k));
                nanSSCPR = isnan(SSCPR(:,:,k));
                
                if sum(nanSSCPH(:)) == 0 || sum(nanSSCPR(:)) == 0
                    
                    V(1,1,k) = max(eig(SSCPR(:,:,k) \ SSCPH(:,:,k)));
                    F(1,1,k) = V(1,1,k) * (v - r + q) / r;
                end
            end
            
            if isnan(p) ~= 1
                P(1,1,:) = 1 - spm_Fcdf(F(1,1,:), [df1 df2]);
                P(P(:) == 1) = NaN;
            end
            
        end
            
        %------------------------------------------------------------------
        % Save slice
        %------------------------------------------------------------------
         FImg(1).private.dat(:,:,j) = reshape(F(1,1,:), dimX, dimY);
         PImg(1).private.dat(:,:,j) = reshape(P(1,1,:), dimX, dimY);
        MPImg(1).private.dat(:,:,j) = 1 - reshape(P(1,1,:), dimX, dimY);

        if MRMoptions.SaveMultiStat == 1
            MVImg(1).private.dat(:,:,j) = reshape(V(1,1,:), dimX, dimY);
        end
        
        msg = sprintf(['Contrast ' num2str(i) ': Estimating slice: ' ...
                        num2str(j) '/' num2str(dimZ) '\n']);
                   
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
       
    switch MRM.Options.Thresh.Level
        
        case 'Voxel'
            %==============================================================
            % FDR correction
            %==============================================================
       
            Pvals  = reshape(PImg(1).private.dat(:,:,:), nVox * dimZ, 1);
            Qvals = MRM_qValue(Pvals,MRMoptions.pFDR,i,MRMoptions.BootPiZero, ...
                               MRMoptions.BootPiResamps);
            Qvals = reshape(Qvals, nVox, dimZ, 1);
            QImg(1).private.dat(:,:,:) = reshape(Qvals, dimX, dimY, dimZ);
            MQImg(1).private.dat(:,:,:) = 1 - reshape(Qvals, dimX, dimY, dimZ);

            % Use the q-values to threshold the F image
            FDRF    = FImg(1).private.dat(:,:,:);
            QvalImg = QImg(1).private.dat(:,:,:);
            FDRF(QvalImg >= alpha) = NaN;

            % Save to disc
            FDRFImg(1).private.dat(:,:,:) = FDRF;
            
        case 'Uncorrected'
            %==============================================================    
            % Uncorrected thresholding
            %==============================================================
        
            % Threshold the p-values
            UncorrectedP = PImg(1).private.dat(:,:,:);
            UncorrectedP(UncorrectedP >= alpha) = NaN;

            % Use the p-values to threshold the F image
            ThreshF = FImg(1).private.dat(:,:,:);
            ThreshF(isnan(UncorrectedP)) = NaN;

            % Save to disc
            ThreshFImg(1).private.dat(:,:,:) = ThreshF;
    end
    
    %======================================================================
    % Save DF to text file
    %======================================================================
    fileID = fopen([Out 'MRM_' name '_DF.txt'], 'w');
    fprintf(fileID, '%s\n', num2str(df1));
    fprintf(fileID, '%s\n', num2str(df2));
    fclose(fileID);
    
    
    %======================================================================
    % Results file
    %======================================================================
    switch MRM.Options.Thresh.Level
        case 'Voxel'
            fprintf(1, ['Contrast ' num2str(i) ': Writing results table ...\n']);
            fprintf(1, '\n');
            MRM_resultsTable(FDRFImg(1).fname, QImg(1).fname, 'Contrasts', name, 'FDR');
        case 'Uncorrected'
            fprintf(1, ['Contrast ' num2str(i) ': Writing results table ...\n']);
            fprintf(1, '\n');
            MRM_resultsTable(ThreshFImg(1).fname, PImg(1).fname, 'Contrasts', name, 'none');
        case 'None'
            fprintf(1, '\n');
    end
   
        
    %======================================================================
    % Tidy up
    %======================================================================
    clearvars SSCPHImg SSCPRImg FImg PImg MPImg;

    if MRMoptions.SaveMultiStat == 1        
        clearvars MVImg       
    end
    
    switch MRM.Options.Thresh.Level
        case 'Voxel'
            clearvars FDRFImg FDRP FDRF Pvals thresh QvalImg QImg MQImg
        case 'Uncorrected'
            clearvars ThreshFImg UncorrectedP ThreshF 
    end
    
end

save([MRM.Options.Out filesep 'MRM.mat'], 'MRM');

end