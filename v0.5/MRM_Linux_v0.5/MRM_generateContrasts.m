function MRM_generateContrasts
%=========================================================================
% Generation of standard ANOVA contrasts             
%=========================================================================
% Function description                                                    
%-------------------------------------------------------------------------
% The methods used here are a bit brute-force as the contrast for each
% scenario is generated separately. I dare say there is a better way to do
% this programmatically but, generally speaking, for up to a 4-way 
% interaction tests the contrasts can be generated using
%
% kron(kron(kron(a,b),c),d)
%
% where a,b,c,d are either the contrasts of the differences of factors 1-4
% given by the MATLAB diff() command, or a row vector of zeros the size of
% the number of levels of the factors 1-4. If the contrast of interest
% contains a particular factor then the diff() contrast is used, otherwise
% zeros are used to pad the expression.
%
% For example, the main effect of Factor 2 contains only factor 2 and so in
% the expression above b = diff(b), with a = zeros(a,1), c = zeros(c,1)
% etc. For interactions you just include more than one diff() command. The
% maximum 4-way interaction is therefore formed by each a:d being a diff()
% command.
%
% This function was inspired by the algorithm found in the 
% spm_make_contrasts() function written by Will Penny at the FIL.
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

global MRM

warning = 0;

if MRM.Factors.Between.Number + MRM.Factors.Within.Number > 4
    warning = 1;
end

nCon = MRM.Contrasts.Number;

k = [];

if MRM.Factors.Between.Number == 0   
    k = 1;
elseif MRM.Factors.Between.Number > 4
    warning = 1;
    for i = 1:4  
        k = [k MRM.Factors.Between.Factors{i}.LevelsNum]; 
    end
else
    for i = 1:MRM.Factors.Between.Number  
        k = [k MRM.Factors.Between.Factors{i}.LevelsNum]; 
    end
end

k  = [k(:)' 1 1 1];

C1a = ones(k(1),1);
C2a = ones(k(2),1);
C3a = ones(k(3),1);
C4a = ones(k(4),1);
D1a = -diff(eye(k(1)))';
D2a = -diff(eye(k(2)))';
D3a = -diff(eye(k(3)))';
D4a = -diff(eye(k(4)))';

k = [];

if MRM.Factors.Within.Number == 0   
    k = 1;
elseif MRM.Factors.Within.Number > 4
    warning = 1;
    for i = 1:4  
        k = [k MRM.Factors.Within.Factors{i}.LevelsNum]; 
    end
else
    for i = 1:MRM.Factors.Within.Number  
        k = [k MRM.Factors.Within.Factors{i}.LevelsNum]; 
    end
end

k  = [k(:)' 1 1 1];

C1c = ones(k(1),1);
C2c = ones(k(2),1);
C3c = ones(k(3),1);
C4c = ones(k(4),1);
D1c = -diff(eye(k(1)))';
D2c = -diff(eye(k(2)))';
D3c = -diff(eye(k(3)))';
D4c = -diff(eye(k(4)))';

if warning == 1
    msgbox(['The auto generation feature is designed for a maximum '       ...
            '4-way interaction across the between and within-subject '     ...
            'factors. The maximum number of contrasts will be generated, ' ...
            'however, a number of the highest-order effects will need to ' ...
            'be entered manually due to the large number of contrasts '    ...
            'necessary to test all the effects.']);
end

%==========================================================================
% Effect of Task
%==========================================================================
MRM.Contrasts.Con{nCon + 1}.Name = 'Task';
MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');

if strcmp(MRM.Model, 'Repeated') == 1
    MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
else
    MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
end

MRM.Contrasts.Number = MRM.Contrasts.Number + 1;
nCon                 = MRM.Contrasts.Number;

%==========================================================================
% Between-subjects Contrasts (up to a 4-way interaction)
%==========================================================================
if MRM.Factors.Between.Number > 0
    
    %----------------------------------------------------------------------
    % Main effect of BS Factor 1
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Between.Factors{1}.Name;
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),C3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;
    nCon                 = MRM.Contrasts.Number; 
    
end

if MRM.Factors.Between.Number > 1  
    
    %----------------------------------------------------------------------
    % Main effect of BS Factor 2
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Between.Factors{2}.Name;
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;
    nCon                 = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % Interaction of BS Factors 1 and 2
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{2}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;
    nCon                 = MRM.Contrasts.Number; 
end

if MRM.Factors.Between.Number > 2
    
    %----------------------------------------------------------------------
    % Factor 3
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Between.Factors{3}.Name;
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),D3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;   
    nCon                 = MRM.Contrasts.Number;    

    %----------------------------------------------------------------------
    % BS Factor 1 x BS Factor 3
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{3}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),D3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;       
    nCon                 = MRM.Contrasts.Number;

    %----------------------------------------------------------------------
    % BS Factor 2 x BS Factor 3
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Between.Factors{3}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),D3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;     
    nCon                 = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factors 1 x BS Factor 2 x BS Factor 3
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Between.Factors{3}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),D3a),C4a)');
    
    if strcmp(MRM.Model, 'Repeated') == 1
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    else
        MRM.Contrasts.Con{nCon + 1}.C = num2cell(eye(MRM.Factors.Within.Factors{1}.LevelsNum));
    end
    
    MRM.Contrasts.Number = MRM.Contrasts.Number + 1;
    nCon                 = MRM.Contrasts.Number;

end

if MRM.Factors.Between.Number > 3
   
    %----------------------------------------------------------------------
    % BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Between.Factors{4}.Name;
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 1 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),C3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;             
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 2 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;                 
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 3 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{3}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),D3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;                
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 1 x BS Factor 2 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name       ...
                                        ' x ' MRM.Factors.Between.Factors{2}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name ];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;              
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 1 x BS Factor 3 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name       ...
                                        ' x ' MRM.Factors.Between.Factors{3}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name ];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),D3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;               
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 2 x BS Factor 3 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name       ...
                                        ' x ' MRM.Factors.Between.Factors{3}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name ];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),D3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;               
    
    nCon = MRM.Contrasts.Number;
    
    %----------------------------------------------------------------------
    % BS Factor 1 x BS Factor 2 x BS Factor 3 x BS Factor 4
    %----------------------------------------------------------------------
    MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name       ...
                                        ' x ' MRM.Factors.Between.Factors{2}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{3}.Name ...
                                        ' x ' MRM.Factors.Between.Factors{4}.Name ];
    MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),D3a),D4a)');
    MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
    MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;                
    
    nCon = MRM.Contrasts.Number;
end


if strcmp(MRM.Model, 'Repeated') == 1

    %==========================================================================
    % Within-subject Contrasts (up to 4-way interaction)
    %==========================================================================
    if MRM.Factors.Within.Number > 0

        % Main effect of WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Within.Factors{1}.Name;
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number; 
    end

    if MRM.Factors.Within.Number > 1  

        % Main effect of WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Within.Factors{2}.Name;
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % Interaction of WS Factors 1 and 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number; 
    end

    if MRM.Factors.Within.Number > 2

        % Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Within.Factors{3}.Name;
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),D3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;   
        nCon = MRM.Contrasts.Number;    

        % WS Factor 1 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),D3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1; 
        nCon = MRM.Contrasts.Number;

        % WS Factor 2 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),D3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factors 1 x WS Factor 2 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),D3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

    end

    if MRM.Factors.Within.Number > 3

        % WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = MRM.Factors.Within.Factors{4}.Name;
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 1 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{4}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 2 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{4}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),C3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 3 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{4}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),D3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 1 x WS Factor 2 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{4}.Name ];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),C3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 1 x WS Factor 3 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{4}.Name ];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),D3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 2 x WS Factor 3 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{4}.Name ];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),D3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % WS Factor 1 x WS Factor 2 x WS Factor 3 x WS Factor 4
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{4}.Name ];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),D3c),D4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end

    %==========================================================================
    % Between-subject x Within-subject Contrasts
    %==========================================================================
    %
    % We generate a maximum 4-way interaction across these terms.  This is done
    % purely for simplicity. Any higher and the number of contrasts necessary 
    % becomes restrictively large
    %
    %==========================================================================

    if MRM.Factors.Between.Number == 1 && MRM.Factors.Within.Number == 1

        % BS Factor 1 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,C2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end

    if MRM.Factors.Between.Number == 1 && MRM.Factors.Within.Number == 2

        % BS Factor 1 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,C2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(C1c,D2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,D2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end

    if MRM.Factors.Between.Number == 1 && MRM.Factors.Within.Number == 3

        % BS Factor 1 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,C2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(C1c,D2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(C1c,C2c),D3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,D2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 1 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,C2c),D3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 2 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(C1c,D2c),D3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 1 x WS Factor 2 x WS Factor 3
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{3}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,D2c),D3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end

    if MRM.Factors.Between.Number == 2 && MRM.Factors.Within.Number == 1

        % BS Factor 1 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(D1a,C2a),C3a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(D1c,C2c),C3c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 2 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x BS Factor 2 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end

    if MRM.Factors.Between.Number == 2 && MRM.Factors.Within.Number == 2

        % BS Factor 1 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 2 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 2 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x WS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 2 x WS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x BS Factor 2 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x BS Factor 2 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;  

        % BS Factor 1 x BS Factor 2 x WS Factor 1 x WS Factor 2
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{2}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end

    if MRM.Factors.Between.Number == 3 && MRM.Factors.Within.Number == 1

        % BS Factor 1 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 2 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 3 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,C2a),D3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x BS Factor 2 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),C3a),C4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x BS Factor 3 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,C2a),D3a),D4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 2 x BS Factor 3 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Between.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(C1a,D2a),D3a),D4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;

        % BS Factor 1 x BS Factor 2 x BS Factor 3 x WS Factor 1
        MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Between.Factors{1}.Name ' x ' MRM.Factors.Between.Factors{2}.Name ' x ' MRM.Factors.Between.Factors{3}.Name ' x ' MRM.Factors.Within.Factors{1}.Name];
        MRM.Contrasts.Con{nCon + 1}.A    = num2cell(kron(kron(kron(D1a,D2a),D3a),D4a)');
        MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
        MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
        nCon = MRM.Contrasts.Number;
    end
end

%==========================================================================
% Covariate contrasts
%==========================================================================

if MRM.Covariates.Number > 0
   
   % Pad existing contrasts with zeros for the covariate(s)  
   for i = 1:size(MRM.Contrasts.Con,2)
      
       A     = cell2mat(MRM.Contrasts.Con{i}.A);
       pad   = zeros(size(A,1),1);
       for j = 1:MRM.Covariates.Number    
           A = [A pad];  
       end
       MRM.Contrasts.Con{i}.A = num2cell(A);
   end
    
   sizeX = size(MRM.Design.X.X, 2);
   
   for i = 1:MRM.Covariates.Number
       
       %-------------------------------------------------------------------
       % Main effect of the CV
       %-------------------------------------------------------------------
       MRM.Contrasts.Con{nCon + 1}.Name = MRM.Covariates.CV{i}.Name;
       MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                    1 ...
                                                    zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
       MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),C3c),C4c)');
       MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
       nCon = MRM.Contrasts.Number;
       
       if strcmp(MRM.Model, 'Repeated') == 1
      
           if MRM.Factors.Within.Number > 0

               %---------------------------------------------------------------
               % WS Factor 1 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),C3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
               nCon = MRM.Contrasts.Number;
           end

           if MRM.Factors.Within.Number > 1

               %---------------------------------------------------------------
               % WS Factor 2 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{2}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),C3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
               nCon = MRM.Contrasts.Number;

               %---------------------------------------------------------------
               % WS Factor 1 x WS Factor 2 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ...
                                                   ' x ' ...
                                                   MRM.Factors.Within.Factors{2}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),C3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
               nCon = MRM.Contrasts.Number;

           end

           if MRM.Factors.Within.Number > 2

               %---------------------------------------------------------------
               % WS Factor 3 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{3}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,C2c),D3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
               nCon = MRM.Contrasts.Number;

               %---------------------------------------------------------------
               % WS Factor 1 x WS Factor 3 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ...
                                                   ' x ' ...
                                                   MRM.Factors.Within.Factors{3}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,C2c),D3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;
               nCon = MRM.Contrasts.Number;

               %---------------------------------------------------------------
               % WS Factor 2 x WS Factor 3 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{2}.Name ...
                                                   ' x ' ...
                                                   MRM.Factors.Within.Factors{3}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(C1c,D2c),D3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;

               %---------------------------------------------------------------
               % WS Factor 1 x WS Factor 2 x WS Factor 3 x CV
               %---------------------------------------------------------------
               MRM.Contrasts.Con{nCon + 1}.Name = [MRM.Factors.Within.Factors{1}.Name ...
                                                   ' x ' ...
                                                   MRM.Factors.Within.Factors{2}.Name ...
                                                   ' x ' ...
                                                   MRM.Factors.Within.Factors{3}.Name ...
                                                   ' x ' ...
                                                   MRM.Covariates.CV{i}.Name];
               MRM.Contrasts.Con{nCon + 1}.A    = num2cell([zeros(1,MRM.Covariates.CV{i}.Col - 1) ... 
                                                            1 ...
                                                            zeros(1,sizeX - MRM.Covariates.CV{i}.Col)]);
               MRM.Contrasts.Con{nCon + 1}.C    = num2cell(kron(kron(kron(D1c,D2c),D3c),C4c)');
               MRM.Contrasts.Number             = MRM.Contrasts.Number + 1;

           end

           if MRM.Factors.Within.Number > 3

               msgbox(['Not all possible contrasts have been generated for ' ...
                       'the covariate ' MRM.Covariates.CV{i}.Name]);

           end
       end
   end
end

