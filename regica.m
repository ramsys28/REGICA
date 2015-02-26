% regica() - Reject ocular artifacts from EEG signals using the reg-ica
%            methodology
%  
%
% Usage:
%   >> Y = regica(X,opt); 
%
%
% Inputs
%   X              - Input data matrix, d x L (d=#sensors, L=#sample points)
%   opt.icatype    - ['runica'|'binica'|'jader'|'acsobiro'|'sobi'] ICA algorithm 
%                    to use for the ICA decomposition. (DEFAULT: runica)
%   opt.extica     - ['yes'|'no'] Used the extended version of ICA or
%                    BINICA (DEFAULT:'yes')
%   opt.regtype    - ['LMS'|'CRLS'|'SRLS'|'H_inf_TV'|'H_inf_EW'|'Schlogl_algo'] 
%                    Regression algorithm for the decontamination of the
%                    ICA components. (DEFAULT: SRLS)
%   opt.crittype   - ['fractal dimension' | 'correlation'] Criterion
%                    algorithm to choose the artifactual ICs
%   opt.corthr     - Correlation threshold to filter the artifactual ICs
%   opt.EOG        - Input EOG data matrix k x L (k=#EOG Signals, L=#sample
%                    points)
%   opt.M          - Order of the adaptive filter (DEFAULT: 3). Used on
%                    LMS,CRLS,SRLS,H_inf_TV,H_inf_EW
%   opt.mu         - Learning rate (DEFAULT: 1e-6). Used only in LMS
%   opt.lambda     - Forgetting factor (DEFAULT: 0.9999). Used only in CRLS and
%                    H_inf_EW
%   opt.sigma      - Parameter described in [1] (DEFAULT: 0.01). Used only in
%                    CRLS and SRLS
%   opt.prec       - Precision (in bits) to use for the computations (DEFAULT:
%                    50). Used only in SRLS.
%   opt.eta        - Factor reflecting a priori knowledge of how close the
%                    estimated filter weights at t=0 are to their optimal
%                    value at that time instant (DEFAULT: 5e-3). Used in
%                    both
%                    H_inf algorithms.
%   opt.rho        - Factor reflecting a priori knowledge of how rapidly
%                    the filter coefficients vary with time (DEFAULT: 1e-5). 
%                    Used in both H_inf algorithms.
%   opt.eps        - Positive constant described in [4] (DEFAULT: 1.5). Used in both
%                    H_inf algorithms.
%
%   You can copy and paste the following code in order to make the exact
%   opt structure that you need, regarding the regression algorithm you
%   want to employ. Note that the only thing that you have to complete is
%   the opt.EOG with the reference EOG signals as described above. 
%
%   opt structure for LMS. 
%   ----------------------
%   opt.M=3;
%   opt.mu=1e-6;
%
%   opt structure for CRLS. 
%   ---------------------------
%   opt.M=3;
%   opt.lambda=0.9999;
%   opt.sigma=0.01;
%
%   opt structure for SRLS. 
%   ---------------------------
%   opt.M=3;
%   opt.lambda=0.9999;
%   opt.sigma=0.01;
%   opt.prec=50;
%  
%   opt structure for H_inf_TV. 
%   ---------------------------
%   opt.M=3;
%   opt.eta=5e-3;
%   opt.rho=1e-5;
%   opt.eps=1.5;
%
%   opt structure for H_inf_EW. 
%   ---------------------------
%   opt.M=3;
%   opt.lambda=0.99;
%   opt.eta=5e-3;
%   opt.rho=1e-5;
%   opt.eps=1.5;
%
%
% Note:
% 1) Infomax (runica, binica) is the ICA algorithm we use most. It is based 
%    on Tony Bell's infomax algorithm as implemented for automated use by 
%    Scott Makeig et al. using the natural gradient of Amari et al. It can 
%    also extract sub-Gaussian sources using the (recommended) 'extended' option 
%    of Lee and Girolami. Function runica() is the all-Matlab version; function 
%    binica() calls the (1.5x faster) binary version (a separate download) 
%    translated into C from runica() by Sigurd Enghoff.
% 2) jader() calls the JADE algorithm of Jean-Francois Cardoso. This is 
%    included in the EEGLAB toolbox by his permission. See >> help jader
% 3) To run fastica(), download the fastICA toolbox from its website,
%    http://www.cis.hut.fi/projects/ica/fastica/, and make it available 
%    in your Matlab path. According to its authors, default parameters
%    are not optimal: Try args 'approach', 'sym' to estimate components 
%    in parallel.
% 4) The function LMS implements the Least Mean Square Algorithm described in [2].
% 5) CRLS performs automatic EOG artifact correction using
%    multiple adaptive regression. The adaptation is made using the Conventional
%    Recursive Least Squares Algorithm (CRLS) [1,2]. A forgetting factor can be
%    used for dealing with time-varying scenarios. 
% 6) SRLS is a stable (but slower) version of the CRLS algorithm.The precision of the
%    computations is set so that stability is guaranteed [3]. If stability is
%    not an issue CRLS does the same job much faster.
% 7) The functions H_inf_EW and H_inf_TV implement the H infinity norm EW and 
%    H infinity norm TV algorithms  described in [4].  
%
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%   H   - Final filter weights (M*k x d)
%   Hh  - filter weights evolution (M*k x d x L)
%
% 
% Citations: 
% 
% 1) M.A.Klados, et al.,REG-ICA:A hybrid methodology combining Blind
% Source Separation and regression techniques for the rejection of ocular artifacts, 
% Biomed. Signal Process. Control (2011), doi:10.1016/j.bspc.2011.02.001.
% 
% 2) Klados, M.A. Papadelis, C.L. Bamidis, P.D.,"REG-ICA: A new hybrid
% method for EOG Artifact Rejection",  9th International Conference on 
% Information Technology and Applications in Biomedicine, 2009. ITAB 2009.doi:10.1109/ITAB.2009.5394295 
%
%
% Author: Manousos Klados 
%         Group of Applied Neurosciences
%         Lab of Medical Informatics
%         Medical School
%         Aristotle University of Thessaloniki
%         P.O. Box 323, 54124 Thessaloniki
%         Greece
% 
%         e-mail: mklados (at) med.auth.gr
%         Τel: +30-2310-999332
%         Fax: +30-2310-999263
%         Personal Page: http://lomiweb.med.auth.gr/gan/mklados/
% 
%
% Aknowledgments
% This software package includes parts from the EEGLAB toolbox, as well as 
% from the Automatic Artifact Removal v 1.3 (AAR1.3) plugin for EEGLAB developed 
% by German Gomez-Herrero.  
% 
% 
% 
% Copyright (C) 2011 Manousos Klados, Aristotle University of Thessaloniki, Greece 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Some parts of the toolbox are owned by others and therefore they might have
% their own licenses. Please visit the related urls for license information
% and the latest versions of the code:
%
%
% References & Credits:
% [1] P. He, G. Wilson, C. Russell, Removal of ocular artifacts from 
%     electroencephalogram by adaptive filtering, Med. Biol. Eng. Comput. 42 (2004) 407-412.
% [2] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
% [3] A. P. Liavas and P. A. Regalia, "On the numerical stability and accuracy of the conventional 
%     RLS algorithm," IEEE Trans. Signal Processing, vol. 47, no. 1, pp. 88-96, January 1999
% [4] S. Puthusserypady and T. Ratnarajah. 2006. Robust adaptive techniques for minimization 
%     of EOG artefacts from EEG signals. Signal Process. 86, 9 (September 2006), 2351-2363. 
%
% See also: hinfew_regression(), lms_regression(),
% scrls_regression(),crls_regression(),
% Schlogl_algo(),pop_regica(),runica(),binica(),sobi(),acsobiro(),jader()

function  [cl_data, cl_ics, regica_weights1, regica_sphere1]=regica(X,opt)

if nargin < 1, help regica; return; end

if ~exist('opt','var'),
    opt = def_regica;
else
    opt = def_regica(opt);
end
if isempty(opt.EOG),
    error('regica(): I need at least one EOG signal!');
end

bss = lower(opt.icatype);
reg=lower(opt.regtype);

if isfield(opt, 'extica')
    opt.extica=lower(opt.extica);
end;



if(strcmp(bss,'runica')==1)
    
    if(strcmp(opt.extica,'yes')==1)
        disp('regica(): Please wait...ICA tries to decompose your signals to ICs...');
        [regica_weights,regica_sphere,compvars,bias,signs,lrates,activations]=runica(X,'extended',1);
    else
        disp('regica(): Please wait...Extended ICA tries to decompose your signals to ICs...');
        [regica_weights,regica_sphere,compvars,bias,signs,lrates,activations]=runica(X);
    end;
end;

if(strcmp(bss,'binica')==1)
    if(strcmp(opt.ext,'yes')==1)
        disp('regica(): Please wait...BINICA tries to decompose your signals to ICs...');
        [regica_weights regica_sphere]=runica(X,'extended',1);
        activations=regica_weights*regica_sphere*X;
    else
        disp('regica(): Please wait...Extended BINICA tries to decompose your signals to ICs...');
        [regica_weights regica_sphere]=runica(X);
        activations=regica_weights*regica_sphere*X;
    end;
end;

if(strcmp(bss,'sobi')==1)
    disp('regica(): Please wait...SOBI tries to decompose your signals to ICs...');
   [winv  activations]=sobi(X);
   regica_weights=pinv(winv);
   regica_sphere=eye(min(size(X))); 
end;

if(strcmp(bss,'acsobiro')==1)
    disp('regica(): Please wait...ACSOBIRO tries to decompose your signals to ICs...');
   [winv activations]=acsobiro(X);
   regica_weights=pinv(winv);
   regica_sphere=eye(min(size(X))); 
end;
if(strcmp(bss,'jader')==1)
    disp('regica(): Please wait...JADE tries to decompose your signals to ICs...');
    regica_weights = jader(X);
    regica_sphere=eye(min(size(X))); 
    activations=regica_weights*regica_sphere*X;
end;



opt.refdata=opt.EOG;

contICs=[];

   switch lower(opt.crittype)
    case 'fractal dimension'
          contICs=regica_fd(activations);
          
          
    case 'correlation'
          opt1.eogref=opt.EOG;
          opt1.corrcoef=opt.corthr;
          contICs=regica_corr(activations,opt1);
          
    end;
regica_Y1=activations;

if(strcmp(reg,'lms')==1)
    disp('regica(): Filtering the EOG artifacts from the ICs using LMS...');
    [regica_Y1(contICs,:)]=lms_regression(activations(contICs,:), opt);
    regica_Y=inv(regica_weights*regica_sphere)*regica_Y1;
    disp('regica(): EOG Artifacts are rejected.');
end;

if(strcmp(reg,'crls')==1)
    disp('regica(): Filtering the EOG artifacts from the ICs using CRLS...');
    [regica_Y1(contICs,:)]=crls_regression(activations(contICs,:), opt);
    regica_Y=inv(regica_weights*regica_sphere)*regica_Y1;
    disp('regica(): EOG Artifacts are rejected.');
end;

if(strcmp(reg,'srls')==1)
    disp('regica(): Filtering the EOG artifacts from the ICs using SRLS...');
    [regica_Y1(contICs,:)]=scrls_regression(activations(contICs,:), opt);
    regica_Y=inv(regica_weights*regica_sphere)*regica_Y1;
    disp('regica(): EOG Artifacts are rejected.');
end;

if(strcmp(reg,'h_inf_tv')==1)
    disp('regica(): Filtering the EOG artifacts from the ICs using Hinf-TV...');
    [regica_Y1(contICs,:)]=hinftv_regression(activations(contICs,:), opt);
    regica_Y=inv(regica_weights*regica_sphere)*regica_Y1;
    disp('regica(): EOG Artifacts are rejected.');
end;
if(strcmp(reg,'h_inf_ew')==1)
    disp('regica(): Filtering the EOG artifacts from the ICs using Hinf-EW...');
    [regica_Y1(contICs,:)]=hinfew_regression(activations(contICs,:), opt);
    regica_Y=pinv(regica_weights*regica_sphere)*regica_Y1;
    disp('regica(): EOG Artifacts are rejected.');
end;
if(strcmp(reg,'schlogl_algo')==1)
    disp('regica(): Filtering the EOG artifacts from the ICs using Schlogl''s Algorithm...');
    [regica_Y1(contICs,:)]=Schlogl_algo(activations(contICs,:), opt.EOG);
    regica_Y=pinv(regica_weights*regica_sphere)*regica_Y1;
    disp('regica(): EOG Artifacts are rejected.');
end;


cl_data=regica_Y; cl_ics=regica_Y1; regica_weights1=regica_weights; regica_sphere1=regica_sphere;

function [opt] = def_regica(opt)
if ~exist('opt','var') || isempty(opt) || ~isfield(opt,'EOG'),
    opt.EOG = [];
end
if ~exist('opt','var') || ~isfield(opt, 'icatype'),
    opt.icatype = 'runica';
end
if ~exist('opt','var') || (((strcmp(opt.icatype, 'runica')==1)||(strcmp(opt.icatype, 'binica')==1))&&(~isfield(opt, 'extica')))
    opt.extica = 'yes';
end
if ~isfield(opt, 'regtype') || ~isfield(opt, 'regtype')
    opt.regtype = 'SRLS';
end
if ~isfield(opt,'sigma') || isempty(opt.sigma),
    opt.sigma = 0.01;
end
if ~isfield(opt, 'lambda') || isempty(opt.lambda),
    opt.lambda = 0.999;
end
if ~isfield(opt, 'M') || isempty(opt.M),
    opt.M = 3;
end
if ~isfield(opt, 'prec') || isempty(opt.prec),
    opt.prec = 50;
end

