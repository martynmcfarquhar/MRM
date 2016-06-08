function MRM_estimate
%=========================================================================
% Estimate the parameters in a multivariate linear model                 
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% This function estimates the parameters, residuals, and covariance 
% structure in the multivariate linear model given by Y = XB + E, where 
%
% Y is n x k
% X is n x m 
% B is m x k
% E is n x k
%
% Given a set of n images for each of k dependent variables (DVs), and a
% between-subject grouping structure encoded in the design matrix X, this
% function will apply the least-squares estimates inv(X'*X)*X'*Y at each
% voxel. This always has a unique solution so long as inv(X'*X) exists. In
% MRM cell-means coded design matrices are used exclusively to guarentee
% estimability.
%
% A useful feature of multivariate linear models is that the parameter
% estimates are identical to multiple univariate models. In addition, the
% variance estimates in the multivariate covariance matrix are also the
% same as the residual mean square from a univariate model. The results 
% obtained with this function should therefore be identical to multiple 
% univariate models fit with other analysis packages (e.g. SPM, FSL) using 
% a cell-means coded design matrix and assuming i.i.d. errors.
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

%--------------------------------------------------------------------------
% Load MRM options in case any settings have been changed
%--------------------------------------------------------------------------
global MRMoptions

ok = MRM_getOptions();

if ok == 0
    msgbox('There was an error loading the global options.');
    return
end

%--------------------------------------------------------------------------

global MRM

% Clear out any old options if an MRM object from a previous run is still
% loaded
if MRM.Estimated == 1
    if isfield(MRM.Design.X, 'Xlda') == 1
        MRM.Design.X = rmfield(MRM.Design.X, 'Xlda');
    end
    
    if isfield(MRM.Design.X, 'XldaVec') == 1
        MRM.Design.X = rmfield(MRM.Design.X, 'XldaVec');
    end
    
    if isfield(MRM.Contrasts, 'Plots') == 1
        MRM.Contrasts = rmfield(MRM.Contrasts, 'Plots');
    end
end

% Get info
outDir   = [MRM.Options.Out filesep];
residOut = MRMoptions.SaveResid;
allScans = 1;
nDVs     = size(MRM.Design.Y.Cell,2);
cells    = size(MRM.Design.X.Cell,2);
useMask  = 0;

% Gather together the files into 1 cell per DV
files = cell(1,nDVs);

for i = 1:nDVs  
    for j = 1:cells
        files{i} = [files{i} MRM.Data.Y{i}.Cell{j}.Scans];
    end
end

nSubs = size(files{1},2); % Number of subjects is taken from the first list 
                          % of files
X    = MRM.Design.X.X;
nPEs = size(X,2);

invXX = inv(X' * X);

%=========================================================================%
% File setup
%=========================================================================%

clc;
fprintf(1, 'Initialising files ...\n \n');

%-------------------------------------------------------------------------%
% Image info
%-------------------------------------------------------------------------%

% Read in first image to get info
img      = spm_vol(files{1}{1});
dimX     = img(1).dim(1);
dimY     = img(1).dim(2);
dimZ     = img(1).dim(3);
sliceDim = dimX * dimY;


%=========================================================================%
%                               FILE CHECK
%=========================================================================%
 for i = 1:nDVs
   
    for j = 2:nSubs
       
        img = spm_data_hdr_read(files{i}{j});
        
        % Check the image dimensions
        if img(1).dim(1) ~= dimX || img.dim(2) ~= dimY || img.dim(3) ~= dimZ
            [~,name,ext] = fileparts(files{i}{j}); 
            msgbox(['Error: The file ', name, ext, ' has different dimensions from the first image.']);
            return
        end
        
    end
 end

% Transformation matrix
T = img(1).mat;

%=========================================================================%
%                               MASK
%=========================================================================%
if isempty(MRM.Options.Mask) ~= 1
   
    try
        mask = spm_vol(cellstr(MRM.Options.Mask));
    catch error
        msgbox(['There was a problem loading the mask file: ' error.message]);
        return
    end
    
    useMask = 1;
    
    fprintf(1, 'Reslicing mask ...\n \n');
    
    % Reslice the mask to the same dimensions as the images
    P = {files{1}{1} MRM.Options.Mask};
    
    flags.mask      = true;
    flags.mean      = false;
    flags.interp    = 0;
    flags.which     = 1;
    flags.wrap      = [0 0 0];
    flags.prefix    = 'resliced_';
    
    spm_reslice(P, flags);
    
    % Load and save the mask data
    [path, name, ext]      = fileparts(MRM.Options.Mask);         
    mask                   = spm_vol([path filesep 'resliced_' name ext]);
    mask                   = mask.private.dat(:,:,:);
    mask(abs(mask) <= eps) = NaN;
    mask(isnan(mask) == 0) = 1; 
    
    % Save image
    MaskImg(1) = deal(struct(...
                             'fname',   [outDir 'MRM_mask.nii'], ...
                             'dim',     [dimX, dimY, dimZ],      ...
                             'dt',      [spm_type('float32')        ...
                                         spm_platform('bigend')],   ...
                             'mat',     T, ...
                             'pinfo',   [1 0 0]',...
                             'descrip', ''));
    
    MaskImg = spm_data_hdr_write(MaskImg);
    
    MaskImg.private.dat(:,:,:) = mask;
    
    delete([path filesep 'resliced_' name ext]);
    
end


%=========================================================================%
%                           OUTPUT FILES
%=========================================================================%


%-------------------------- Parameter Estimates --------------------------%
outPE(1:size(X,2) * nDVs) = deal(struct(...
                                'fname',   [], ...
                                'dim',     [dimX, dimY, dimZ],          ...
                                'dt',      [spm_type('float64')         ...
                                            spm_platform('bigend')],    ...
                                'mat',     T,                           ...
                                'pinfo',   [1 0 0]',                    ...
                                'descrip', 'MRM - Parameter Estimates'));

k = 1;                        
for i = 1:nDVs
    for j = 1:nPEs
        outPE(k).fname    = [outDir 'MRM_PE_' num2str(j) '_' num2str(i) ...
                            '.nii'];
        outPE(k).descript = ['MRM - Parameter Estimate: row ' num2str(j) ...
                            ' col ' num2str(i)];
        k = k + 1;
    end
end

% Write header
outPE = spm_data_hdr_write(outPE);

                                         
%------------------------------- Residuals -------------------------------%
if residOut == 1
   outResid(1:nSubs * nDVs) = deal(struct(...
                                   'fname',   [], ...
                                   'dim',     [dimX, dimY, dimZ], ...
                                   'dt',      [spm_type('float64') ...
                                               spm_platform('bigend')],...
                                   'mat',     T, ...
                                   'pinfo',   [1 0 0]',...
                                   'descrip', 'MRM - Residuals'));

    k = 1;                        
    for i = 1:nDVs
        for j = 1:nSubs 
            outResid(k).fname    = [outDir 'MRM_Resid_' num2str(j) '_' ...
                                   num2str(i) '.nii'];
            outResid(k).descript = ['MRM - Residuals: row ' num2str(j) ...
                                   ' col ' num2str(i)];
            k = k + 1;
        end
    end
    
    % Write header
    outResid = spm_data_hdr_write(outResid);
end


                                                  
%---------------------- Variance-covariance matrix -----------------------%
outCovar(1:nDVs*nDVs) = deal(struct(...
                            'fname',   [], ...
                            'dim',     [dimX, dimY, dimZ], ...
                            'dt',      [spm_type('float64') ...
                                        spm_platform('bigend')],...
                            'mat',     T, ...
                            'pinfo',   [1 0 0]',...
                            'descrip', 'MRM - Variance-covariance'));
k = 1;                
for i = 1:nDVs
    for j = 1:nDVs
        outCovar(k).fname    = [outDir 'MRM_Covar_' num2str(j), '_' ...
                                num2str(i), '.nii'];
        outCovar(k).descript = ['MRM - Variance-covariance: row ' ... 
                                num2str(j) ' col ' num2str(i)];
        k = k + 1; 
    end
end                

% Write header
outCovar = spm_data_hdr_write(outCovar);


%=========================================================================%
%                             PRINT INFO
%=========================================================================%
fprintf('%s', repmat('=',1,72));
fprintf(1, '\n');
fprintf(1, 'Model estimation\n');
fprintf('%s', repmat('=',1,72));
fprintf(1, '\n');


%=========================================================================%
%                           INITIALISE MATRICES
%=========================================================================%

paramMat   = NaN(dimX, dimY, dimZ, size(X,2),   'double');
outcomeMat = NaN(size(files, 2), sliceDim,      'double');
residMat   = NaN(dimX, dimY, dimZ, nSubs, nDVs, 'double');
estimates  = NaN(size(X,2), nDVs, 'double');

%=========================================================================%
%                               ESTIMATION
%=========================================================================%

p = 0; %-index for the PEs across DVs

for m = 1:nDVs

    reverseStr = '';

    data = cell(nSubs,1);
    
    for i = 1:nSubs
        data{i} = spm_data_read(spm_data_hdr_read(files{m}{i}));
    end

    for i = 1:dimZ

        % Gather slices from across subjects and arrange to form a matrix
        % of subjects x voxels
        
        for j = 1:nSubs    
            slice             = data{j}(:,:,i);
            slice(slice == 0) = NaN;            % Replace zeros with NaN
            outcomeMat(j,:)   = reshape(slice, 1, sliceDim);
        end
        
        % Estimate parameters and residuals      
        if sum(isnan(outcomeMat(:))) ~= sliceDim
            estimates = invXX * X' * outcomeMat;
            resid     = outcomeMat - X * estimates; 
        end

        % Reshape the matrix of estimated parameters into a separate slice
        % per column of X. Save these to the output image(s).
        for k = 1:size(X,2)
            
            if sum(isnan(outcomeMat(:))) ~= sliceDim  
                paramMat(:,:,i,k) = reshape(estimates(k,:), dimX, dimY); 
            else
                paramMat(:,:,i,k) = NaN(dimX, dimY, 'double');
            end
            
            % Mask PEs if provided
            if useMask == 1   
                out                             = paramMat(:,:,i,k);
                out(isnan(mask(:,:,i)))         = NaN;
                outPE(p + k).private.dat(:,:,i) = out;  
            else   
                outPE(p + k).private.dat(:,:,i) = paramMat(:,:,i,k); 
            end
        end
        
        
        % Save residuals for this slice into a 5D matrix of X,Y,Z,S,DV,
        % where S is the number of subjects and DV is the number of DVs
        
        for l = 1:nSubs
            
            if sum(isnan(outcomeMat(:))) ~= sliceDim  
                residMat(:,:,i,l,m) = reshape(resid(l,:), dimX, dimY);
            else
                residMat(:,:,i,l,m) = nan(dimX, dimY);
            end
        end    
       
        msg = sprintf(['Estimating parameters and saving images for DV ' ...
                       num2str(m) ' slice ' num2str(i) '/' num2str(dimZ) '\n']);
                   
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
    end
    
p = p + k;    

end

fprintf(1, '\n');

% Tidy-up
if allScans == 1  
    clear data;
end


%=========================================================================%
%                               SAVE RESIDUALS
%=========================================================================%

if residOut == 1

    p = 0;
    
    for i = 1:nDVs
        
        reverseStr = '';
        
        for j = 1:dimZ
            
            for k = 1:nSubs
                
                if useMask == 1
                
                    out                                = residMat(:,:,j,k,i);
                    out(isnan(mask(:,:,j)))            = NaN;
                    outResid(p + k).private.dat(:,:,j) = out;
                
                else
                    
                    outResid(p + k).private.dat(:,:,j) = residMat(:,:,j,k,i);
                
                end
                
            end

            msg = sprintf(['Saving residuals for DV: ' num2str(i) ...
                            ' slice: ' num2str(j) '/' num2str(dimZ) '\n']);
                   
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            
        end
        
       p = p + k;
       
    end
    
    fprintf(1, '\n');
end



%=========================================================================%
%                          VARIANCE-COVARIANCE MATRIX
%=========================================================================%
%
% The covariance structure is estimated using the residuals from the
% multivariate linear model Y = XB + E. This estimation takes the form of
% (E' * E) / (n - m).
%
%=========================================================================%
 
reverseStr = '';

for k = 1:dimZ

    resid = zeros(nSubs, nDVs, sliceDim);
    
    % Reshape residuals so they are organised as subs x resid x voxels
    for i = 1:nSubs
        
        for j = 1:nDVs
            
            resid(i,j,:) = reshape(residMat(:,:,k,i,j), 1, sliceDim);
            
        end
    end
    
    sliceCov = zeros(nDVs, nDVs, size(resid,3));
    
    % Estimate covariance matrix for this slice 
    for i = 1:size(resid,3)
        
        sliceCov(:,:,i) = (resid(:,:,i)' * resid(:,:,i)) / (nSubs - nPEs);
        
    end
    
    % Work through each element of the covariance matrix, saving each as a
    % slice in a separate image
    l = 1;
    for i = 1:size(sliceCov,1)
        
        for j = 1:size(sliceCov,2)
            
            % Mask the output if requested
            if useMask == 1
                
                out = reshape(sliceCov(i,j,:), dimX, dimY);
                out(isnan(mask(:,:,k))) = NaN;
                outCovar(l).private.dat(:,:,k) = out;
                
            else
                
                outCovar(l).private.dat(:,:,k) = reshape(sliceCov(i,j,:), ...
                                                     dimX, dimY);
            end
            
            l = l + 1;
            
        end
    end
    
    msg = sprintf(['Estimating and writing covariance structure for slice ' ...
                   num2str(k) '/' num2str(dimZ) '\n']);
                   
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
end

clear residMat
clear resid
clear sliceCov

fprintf(1, '\n');

%=========================================================================%
%               SAVE DESIGN VISUALISATION AND MRM STRUCTURE
%=========================================================================%

fprintf(1, 'Saving design images ...\n \n');
MRM_drawDesign(1);
MRM_designOrthogonality(1);

% Save MRM structure
MRM.Estimated = 1;
save([outDir 'MRM.mat'], 'MRM');

%=========================================================================%
%                               CONTRASTS
%=========================================================================%
if MRM.Contrasts.Number > 0 
    switch MRM.Options.Thresh.Level
        case 'Cluster'
            MRM_estimateContrastsPermClust();
        otherwise
            switch MRM.Options.Thresh.Pvals
                case 'Permutation'
                    MRM_estimateContrastsPerm();
                case 'Approximate'
                    MRM_estimateContrasts();
            end
    end
end

fprintf(1, 'Done!\n');

end