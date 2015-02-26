% pop_regica() - Reject ocular artifacts from EEG signals using the reg-ica
%                methodology
%  
%
% Usage:
%   >> OUT_EEG = pop_regica( EEG ); % pops-up a data entry window
%
%
% Graphic interface:
%   "BSS algorithm to use" - [list box] The BSS algorithm to use for 
%                 Blind Source Separation. Command line equivalent:'icatype'
%   "Commandline options" - [edit box] Command line options to forward
%                 to the BSS algorithm. Command line equivalent: 'options'
%   "EOG channel indices" - [edit box] indices of the EOG reference
%                 channels [e.g 62 63 ...etc] Command line equivalent:
%                 'EOGindex'
%   "Regression algorithm to use" - [list box] The regression algorithm to 
%                 use for the rejection of ocular artifacts. Command line
%                 equivalent:'regtype'
%   'M'         - Order of the adaptive filter (DEFAULT: 3). Used on
%                 LMS,CRLS,SRLS,H_inf_TV,H_inf_EW
%   'mu'        - Learning rate (DEFAULT: 1e-6). Used only in LMS
%   'lambda'    - Forgetting factor (DEFAULT: 0.999 for CRLS,SRLS and 0.9 for 
%                 H_inf_EW). Used only in CRLS,SRLS and H_inf_EW.
%   'sigma'     - Parameter described in [1] (DEFAULT: 0.01). Used only in
%                 CRLS and SRLS
%   'prec'      - Precision (in bits) to use for the computations (DEFAULT:
%                 50). Used only in SRLS.
%   'eta'       - Factor reflecting a priori knowledge of how close the
%                 estimated filter weights at t=0 are to their optimal
%                 value at that time instant (DEFAULT: 5e-3). Used in both
%                 H_inf algorithms.
%   'rho'       - Factor reflecting a priori knowledge of how rapidly
%                 the filter coefficients vary with time (DEFAULT: 1e-5). 
%                 Used in both H_inf algorithms.
%   'eps'       - Positive constant described in [4] (DEFAULT: 1.5). Used in both
%                 H_inf algorithms.
% 
% Inputs:
%   EEG         - input EEG dataset or array of datasets
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
% 4) The function LMS implements the Least Mean Square Algorithm described
% in [2].
% 5) CRLS performs automatic EOG artifact correction using
%    multiple adaptive regression. The adaptation is made using the Conventional
%    Recursive Least Squares Algorithm (CRLS) [1,2]. A forgetting factor can be
%    used for dealing with time-varying scenarios. 
% 6) SRLS is a stable (but slower) version of the CRLS algorithm.The precision of the
%    computations is set so that stability is guaranteed [3]. If stability is
%    not an issue CRLS does the same job much faster.
% 7) The functions H_inf_EW and H_inf_TV implement the H infinity norm EW and 
%    H infinity norm TV algorithms  described in [4].  
% 8) [optional] Install, i.e. include in Matlab path, the following Blind 
%    Source Separation (BSS) algorithms:
%    8.1) JADE: 
%         (c) Jean-Francois Cardoso
%         WWW: http://www.tsi.enst.fr/~cardoso/Algo/Jade/jadeR.m
%    8.2) FastICA:
%	      (c) H. Gδvert, J. Hurri, J. Sδrelδ and A. Hyvδrinen
%	      WWW: http://www.cis.hut.fi/projects/ica/fastica/

% Outputs:
%   OUT_EEG = The input EEGLAB dataset with new fields icaweights, icasphere 
%             and icachansind (channel indices). 
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





function [ALLEEG, com] = pop_regica( ALLEEG, varargin )

com = '';
if nargin < 1   
    help pop_runica;
    return;
end;

%Params for results
parica=1;
paricapar=2;
parcrit=3;
parthr=4;
parreg=5;
pareog=6;
parm=7;
parl=8;
parff=9;
pars=10;
parprec=11;
pareta=12;
parrho=13;
pareps=14;
parchan=15;

% runica
% find available algorithms
% -------------------------
allalgs   = { 'runica' 'binica' 'jader' 'jadeop' 'jade_td_p' 'MatlabshibbsR' 'fastica' ...
              'tica' 'erica' 'simbec' 'unica' 'amuse' 'fobi' 'evd' 'evd24' 'sons' 'sobi' 'ng_ol' ...
              'acsobiro' 'acrsobibpf' 'pearson_ica' 'egld_ica' 'eeA' 'tfbss' 'icaML' 'icaMS' }; % do not use egld_ica => too slow
allregs   = { 'LMS' 'CRLS' 'SRLS' 'H_inf_TV' 'H_inf_EW'};          
allcrit  = { 'Fractal Dimension' 'Correlation' };  
selectalg = {};
linenb    = 1;
count     = 1;
for index = length(allalgs):-1:1
    if exist(allalgs{index}) ~= 2 && exist(allalgs{index}) ~= 6
        allalgs(index) = [];
    end;
end;

% popup window parameters
% -----------------------
fig = [];
if nargin < 2

    cb_ica = [ 'if get(gcbo, ''value'') < 3, ' ...
               '     set(findobj(gcbf, ''tag'', ''params''), ''string'', ''''''extended'''', 1'');' ...
               'else set(findobj(gcbf, ''tag'', ''params''), ''string'', '''');' ...
               'end;' ];

    cb_reg11 = [ 'if get(gcbo, ''value'') ==1 , ' ...
               '     set(findobj(gcbf, ''tag'', ''fo2''), ''string'', ''3'', ''enable'', ''on'' );' ...
               '     set(findobj(gcbf, ''tag'', ''lr2''), ''string'', ''1e-6'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''ff''), ''string'', '''', ''enable'', ''off'' );' ...
               '     set(findobj(gcbf, ''tag'', ''sigma''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''prec''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eta''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''rho''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eps''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''evanchanss''), ''string'', '''', ''enable'', ''on'');' ...
               'end;'...
               '     if get(gcbo, ''value'') ==2 , ' ...
               '     set(findobj(gcbf, ''tag'', ''fo2''), ''string'', ''3'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''lr2''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''ff''), ''string'', ''0.9999'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''sigma''), ''string'', ''0.01'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''prec''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eta''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''rho''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eps''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''evanchanss''), ''string'', '''', ''enable'', ''on'');' ...
               'end;'...
               '     if get(gcbo, ''value'') ==3 , ' ...
               '     set(findobj(gcbf, ''tag'', ''fo2''), ''string'', ''3'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''lr2''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''ff''), ''string'', ''0.9999'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''sigma''), ''string'', ''0.01'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''prec''), ''string'', ''50'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''eta''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''rho''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eps''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''evanchanss''), ''string'', '''', ''enable'', ''on'');' ...
               'end;'...
               '     if get(gcbo, ''value'') ==4 , ' ...
               '     set(findobj(gcbf, ''tag'', ''fo2''), ''string'', ''3'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''lr2''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''ff''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''sigma''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''prec''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eta''), ''string'', ''5e-3'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''rho''), ''string'', ''1e-5'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''eps''), ''string'', ''1.5'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''evanchanss''), ''string'', '''', ''enable'', ''on'');' ...
               'end;'...
               '     if get(gcbo, ''value'') ==5 , ' ...
               '     set(findobj(gcbf, ''tag'', ''fo2''), ''string'', ''3'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''lr2''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''ff''), ''string'', ''0.99'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''sigma''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''prec''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eta''), ''string'', ''5e-3'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''rho''), ''string'', ''1e-5'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''eps''), ''string'', ''1.5'', ''enable'', ''on'');' ...
               '     set(findobj(gcbf, ''tag'', ''evanchanss''), ''string'', '''', ''enable'', ''on'');' ...
               'end;'...
               '     if get(gcbo, ''value'') ==6 , ' ...
               '     set(findobj(gcbf, ''tag'', ''fo2''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''lr2''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''ff''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''sigma''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''prec''), ''string'', '''', ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eta''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''rho''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''eps''), ''string'', '''',  ''enable'', ''off'');' ...
               '     set(findobj(gcbf, ''tag'', ''evanchanss''), ''string'', '''', ''enable'', ''off'');' ...
               'end;'];
    cb_crit = [ 'if get(gcbo, ''value'') ==1 , ' ...
                   'set(findobj(gcbf, ''tag'', ''cor_thres''), ''enable'', ''off'', ''String'', '' '' );'...
                   'else set(findobj(gcbf, ''tag'', ''cor_thres''), ''enable'', ''on'', ''String'', ''0.25'' );'...
                   'end;'];
               


        
   
    


         commandchans = [ 'tmpchans = get(gcbf, ''userdata'');' ...
                     'tmpchans = tmpchans{1};' ...
                     'set(findobj(gcbf, ''tag'', ''chantype''), ''string'', ' ...
                     '       int2str(pop_chansel( tmpchans )));' ...
                     'clear tmpchans;' ];
      

    promptstr    = { { 'style' 'text'       'string' 'BSS Selection' 'fontweight', 'bold'} ...
                     { 'style' 'text'       'string' 'BSS algorithm to use (click to select)' } ...
                     { 'style' 'listbox'    'string' strvcat(allalgs{:}) 'callback', cb_ica } ...
                     { 'style' 'text'       'string' 'Commandline options (See help messages)' } ...
                     { 'style' 'edit'       'string' '''extended'', 1' 'tag' 'params' } ...
                     { 'style' 'text'       'string' 'Choose Criterion for Automatic Identification' 'tag' 'crit_txt' }...
                     { 'style' 'listbox'    'string' strvcat(allcrit{:}) 'tag' 'crit_cho' 'callback' cb_crit } ...
                     { 'style' 'text'       'string' 'Correlation Threshold (0<...<1)'  } ...
                     { 'style' 'edit'       'string' '' 'tag' 'cor_thres' 'enable' 'off'  } ...
                     {}...
                     { 'style' 'text'       'string' 'Regression Selection' 'fontweight', 'bold'} ...
                     { 'style' 'text'       'string' 'Regression algorithm to use (click to select)' } ...
                     { 'style' 'listbox'    'string' strvcat(allregs{:}) 'callback', cb_reg11,} ...
                     { 'style' 'text'       'string' 'EOG channel indices:'} ...
                     { 'style' 'edit'       'string' ''} ...
                     { 'style' 'text'       'string' 'Filter order (M):' 'tag' 'fo1' } ...
                     { 'style' 'edit'       'string' '' 'tag' 'fo2' 'string' '3'  } ...
                     { 'style' 'text'       'string' 'Learning rate (mu):' 'tag' 'lr1' } ...
                     { 'style' 'edit'       'string' '' 'tag' 'lr2' 'string' '1e-6'} ...
                     { 'style' 'text'       'string' 'Forgetting factor (lambda):'} ...
                     { 'style' 'edit'       'string' '' 'tag' 'ff' 'enable' 'off'} ...
                     { 'style' 'text'       'string' 'Sigma:'} ...
                     { 'style' 'edit'       'string' '' 'tag' 'sigma' 'enable' 'off'} ...
                     { 'style' 'text'       'string' 'Precision (in bits) to use for the computations'} ...
                     { 'style' 'edit'       'string' '' 'tag' 'prec' 'enable' 'off'} ...
                     { 'style'  'text'      'string'    'Distance at t=0 to optimal solution (eta):'} ...
                     { 'style'  'edit'      'string'  ''  'tag' 'eta' 'enable' 'off'} ...        
                     { 'style'  'text'      'string'     'Speed of variation of filter coefficients (rho):'} ...
                     { 'style'  'edit'      'string'  ''  'tag' 'rho' 'enable' 'off'} ...
                     { 'style'  'text'      'string'     'Positive constant epsilon:'} ...
                     { 'style'  'edit'      'string'  ''  'tag' 'eps' 'enable' 'off'} ...
                     { 'style' 'text'       'string' 'Store filter weights for channels:'} ...
                     { 'style' 'edit'       'string' '' 'tag' 'evanchanss' 'enable' 'off'} ...
                     };
                geometry = { [1] [2 1] [2 1]  [2 1] [2 1] [1] [1] [2 1] [2 1] [2 1] [2 1] [2 1] [2 1] [2 1]  [2 1] [2 1] [2 1] [2 1] };
                 
    % channel types
    % -------------
    if isfield(ALLEEG(1).chanlocs, 'type'), alltypes = { ALLEEG(1).chanlocs.type };
                                            indempty = cellfun('isempty', alltypes);
                                            alltypes(indempty) = '';
                                            alltypes = unique(alltypes);
    else                                    alltypes = '';
    end;
    
    % channel labels
    % --------------
    if ~isempty(ALLEEG(1).chanlocs)
        alllabels = { ALLEEG(1).chanlocs.labels };
    else
        for index = 1:ALLEEG(1).nbchan
            alllabels{index} = int2str(index);
        end;
    end;
    
    % gui
    % ---
    result       = inputgui( 'geometry', geometry, 'uilist', promptstr, ...
                             'helpcom', 'pophelp(''pop_regica'')', ...
                             'title', 'Artifact Rejection using REG-ICA Methodology -- pop_regica()', 'userdata', { alllabels alltypes } );
    if length(result) == 0 return; end;        
    options = { 'crittype' allcrit{result{parcrit}} 'regtype' allregs{result{parreg}} 'icatype' allalgs{result{parica}} 'dataset' [1:length(ALLEEG)] 'options' eval( [ '{' result{paricapar} '}' ]) };
        if ~isempty(result{pareog})
                if ~isempty(str2num(result{pareog})), options = { options{:} 'chanind' str2num(result{pareog}) };
                else                             options = { options{:} 'chanind' eeg_chantype(ALLEEG(1).chanlocs, parsetxt(result{pareog})) }; 
                end;
        end;

        else 
            if ~strcmpi(varargin{parica}, 'icatype')
                options = { 'icatype' varargin{1:end} };
            else
                options = varargin;
            end;
            if ~strcmpi(varargin{parreg}, 'regtype')
                options = { 'regtype' varargin{1:end} };
            else
                options = varargin;
            end;
            if ~strcmpi(varargin{parcrit}, 'crittype')
                options = { 'crittype' varargin{1:end} };
            else
                options = varargin;
            end;
end;
% decode input arguments
% ----------------------
[ g options ] = finputcheck( options, { 'crittype' 'string' allcrit ''; ...
                            'regtype' 'string' allregs 'LMS'; ...
                            'icatype'        'string'  allalgs   'runica'; ...
                            'dataset'        'integer' []        [1:length(ALLEEG)];
                            'options'        'cell'    []        {};}, ...
                            'pop_regica', 'ignore');
if ischar(g), error(g); end;
if isempty(g.options), g.options = options; end;

% select datasets, create new big dataset if necessary
% ----------------------------------------------------
if length(g.dataset) == 1
    EEG = ALLEEG(g.dataset);
else
    error('pop_regica:REG-ICA can be performed only in one dataset')
end;    

% Store and then remove current EEG ICA weights and sphere
% ---------------------------------------------------
fprintf('\n');
if ~isempty(EEG.icaweights)
    fprintf('Saving current ICA decomposition in "EEG.etc.oldicaweights" (etc.).\n');
    if ~isfield(EEG,'etc'), EEG.etc = []; end;
    if ~isfield(EEG.etc,'oldicaweights')
        EEG.etc.oldicaweights = {};
        EEG.etc.oldicasphere = {};
        EEG.etc.oldicachansind = {};
    end;
    EEG.etc.oldicaweights = { EEG.icaweights  EEG.etc.oldicaweights{:} };
    EEG.etc.oldicasphere  = { EEG.icasphere   EEG.etc.oldicasphere{:}  };
    EEG.etc.oldicachansind  = { EEG.icachansind EEG.etc.oldicachansind{:}  };
    fprintf('               Decomposition saved as entry %d.\n',length(EEG.etc.oldicaweights));
end
EEG.icaweights = [];
EEG.icasphere  = [];
EEG.icawinv    = [];
EEG.icaact     = [];

% select sub_channels
% -------------------
if ~isempty(g.chanind)
    
    EOGindex = eval(['[' result{pareog} ']']);
    if (min(EOGindex)>EEG.nbchan)
        error('pop_regica(): There are not such EOG signals');
        return;
    else
    tmp_difa= 1:EEG.nbchan;
    g.chanind=setdiff(tmp_difa,EOGindex);
    [tmp_m tmp_n]=size(g.chanind);
    g.nbchan=max(tmp_m, tmp_n);
    end;
end;
tmp_dataa=EEG.data(g.chanind,:,:);
EOG=EEG.data(EOGindex,:,:);
EEG.data=zeros(size(tmp_dataa));
EEG.data=tmp_dataa;
EEG.icachansind = g.chanind;
EEG.nbchan=g.nbchan;

% is pca already an option?
% -------------------------
pca_opt = 0;
for i = 1:length(g.options)
    if isstr(g.options{i})
        if strcmpi(g.options{i}, 'pca')
            pca_opt = 1;
        end;
    end;
end;

%------------------------------
% compute ICA on a definite set
% -----------------------------
tmpdata = reshape( EEG.data(g.chanind,:,:), length(g.chanind), EEG.pnts*EEG.trials);
tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean 
if ~strcmpi(lower(g.icatype), 'binica')
    try,
        disp('Attempting to convert data matrix to double precision (more accurate ICA results)')
        tmpdata = double(tmpdata);
        tmpdata2 = tmpdata+1; % check for more memory
        clear tmpdata2;
    catch,
        disp('*************************************************************')
        disp('Not enougth memory to convert data Matrix to double precision')
        disp('All computation will be done in single precision. Matlab 7.x')
        disp('(under 64-bit Linux and others) is imprecise in this mode')
        disp('We advise that you use "binica" instead of "runica"')
        disp('*************************************************************')
    end;
end;
switch lower(g.icatype)
    case 'runica' 
        if nargin < 2
            fig = figure('visible', 'off');
            supergui( fig, {1 1}, [], {'style' 'text' 'string' 'Press button to interrupt runica()' }, ...
                      {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'figure(gcbf); set(gcbf, ''tag'', ''stop'');' } );
            drawnow;
        end;
        tmprank = rank(double(tmpdata(:,1:min(3000, size(tmpdata,2)))));
        if tmprank == size(tmpdata,1) | pca_opt
            [EEG.icaweights,EEG.icasphere,compvars,bias,signs,lrates,activations] = runica( tmpdata, 'lrate', 0.001, g.options{:} );
            EEG.icaact=activations;
        else 
            disp(['Data rank (' int2str(tmprank) ') less than the number of channels (' int2str(size(tmpdata,1)) ').']);
            [EEG.icaweights,EEG.icasphere, compvars,bias,signs,lrates,activations] = runica( tmpdata, 'lrate', 0.001, 'pca', tmprank, g.options{:} );
            EEG.icaact=activations;
        end;
     case 'binica'
        icadefs;
        fprintf(['Warning: If the binary ICA function does not work, check that you have added the\n' ...
                 'binary file location (in the EEGLAB directory) to your Unix /bin directory (.cshrc file)\n']);
        if exist(ICABINARY) ~= 2
            error('pop_regica(): binary ICA executable not found. Edit icadefs.m file to specify the ICABINARY location');
        end;
        tmprank = rank(double(tmpdata(:,1:min(3000, size(tmpdata,2)))));
        if tmprank == size(tmpdata,1) | pca_opt
            [EEG.icaweights,EEG.icasphere,compvars,bias,signs,lrates,activations] = binica( tmpdata, 'lrate', 0.001, g.options{:} );
            EEG.icaact=activations;
        else 
            disp(['Data rank (' int2str(tmprank) ') is less than the number of channels (' int2str(size(tmpdata,1)) ').']);
            [EEG.icaweights,EEG.icasphere,compvars,bias,signs,lrates,activations] = binica( tmpdata, 'lrate', 0.001, 'pca', tmprank, g.options{:} );
            EEG.icaact=activations;
        end;
     case 'pearson_ica' 
        if isempty(g.options)
            disp('Warning: EEG default for pearson ICA changed to 1000 iterations and epsilon=0.0005');
            [tmp EEG.icaweights] = pearson_ica( tmpdata, 'maxNumIterations', 1000,'epsilon',0.0005);
        else    
            [tmp EEG.icaweights] = pearson_ica( tmpdata, g.options{:});
        end;
     case 'egld_ica', disp('Warning: this algorithm is very slow!!!');
                      [tmp EEG.icaweights] = egld_ica( tmpdata, g.options{:} );
     case 'tfbss' 
        if  isempty(g.options)
             [tmp EEG.icaweights] = tfbss( tmpdata, size(tmpdata,1), 8, 512 );
        else    
             [tmp EEG.icaweights] = tfbss( tmpdata, g.options{:} );
        end;
     case 'jader',         [EEG.icaweights] = jader( tmpdata, g.options{:} );
     case 'matlabshibbsr', [EEG.icaweights] = MatlabshibbsR( tmpdata, g.options{:} );
     case 'eea',           [EEG.icaweights] = eeA( tmpdata, g.options{:} );
     case 'icaml',         [tmp EEG.icawinv] = icaML( tmpdata, g.options{:} );
     case 'icams',         [tmp EEG.icawinv] = icaMS( tmpdata, g.options{:} );
     case 'fastica',       [ ICAcomp, EEG.icawinv, EEG.icaweights] = fastica( tmpdata, 'displayMode', 'off', g.options{:} );
     case { 'tica' 'erica' 'simbec' 'unica' 'amuse' 'fobi' 'evd' 'sons' ...
            'jadeop' 'jade_td_p' 'evd24' 'sobi' 'ng_ol' 'acsobiro' 'acrsobibpf' } 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        
        if isempty(g.options), g.options = { size(tmpdata,1) }; end;
        switch lower(g.icatype)
         case 'tica',     EEG.icaweights = tica( tmpdata, g.options{:} );
         case 'erica',    EEG.icaweights = erica( tmpdata, g.options{:} );
         case 'simbec',   EEG.icaweights = simbec( tmpdata, g.options{:} );
         case 'unica',    EEG.icaweights = unica( tmpdata, g.options{:} );
         case 'amuse',    EEG.icaweights = amuse( tmpdata );
         case 'fobi',     [tmp EEG.icaweights] = fobi( tmpdata, g.options{:} );
         case 'evd',      EEG.icaweights = evd( tmpdata, g.options{:} );
         case 'sons',     EEG.icaweights = sons( tmpdata, g.options{:} );
         case 'jadeop',   EEG.icaweights = jadeop( tmpdata, g.options{:} );
         case 'jade_td_p',EEG.icaweights = jade_td_p( tmpdata, g.options{:} );
         case 'evd24',    EEG.icaweights = evd24( tmpdata, g.options{:} );
         case 'sobi',     EEG.icawinv = sobi( EEG.data );
         case 'ng_ol',    [tmp EEG.icaweights] = ng_ol( tmpdata, g.options{:} );
         case 'acsobiro',  [EEG.icawinv] = acsobiro(EEG.data);
         case 'acrsobibpf', EEG.icawinv = acrsobibpf( tmpdata, g.options{:} );
        end;
        clear tmp;
        
        close(fig);
        
     otherwise, error('pop_regica(): unrecognized algorithm');
end;

% update weight and inverse matrices etc...
% -----------------------------------------
if ~isempty(fig), try, close(fig); catch, end; end;
if isempty(EEG.icaweights)
    EEG.icaweights = pinv(EEG.icawinv);
end;
if isempty(EEG.icasphere)
    EEG.icasphere  = eye(size(EEG.icaweights,2));
end;
if isempty(EEG.icawinv)
    EEG.icawinv    = pinv(EEG.icaweights*EEG.icasphere);
end;



if isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data;
end;



   switch lower(g.crittype)
    case 'fractal dimension'
          EEG.contICs=regica_fd(EEG.icaact);
          
    case 'correlation'
        
        
        opt.thr=eval(['[' result{parthr} ']']);
        [m1 n1]=size(EEG.icaact);
        [m2 n2]=size(EOG);
        
        xcorr12=zeros(m1,m2);
        for i=1:m1
            for j=1:m2
            xcorr12(i,j)=max(xcorr(EEG.icaact(i,:),EOG(j,:),'coeff'));
            end;
        end;
        j=1;
        for i=1:m1
            if any(xcorr12(i,:)>opt.thr)
                EEG.contICs(j,1)=i;
                j=j+1;
            end;
        end;
       

     
    end;
    
    if isempty(EEG.contICs)
        error('Artifactual ICs cannot be identified...Try another algorithm')
    end;

switch lower(g.regtype)
         
     case 'lms'
            disp('pop_regica(): LMS algorithm will be used for the rejection of ocular artifacts'); 
            if isempty(EEG.icaact)
                disp('pop_regica(): cannot clean an empty dataset'); 
                return;
            end;


         
        % reading params
        % -------------------
        EOGindex = EOG;
        M = eval(['[' result{parm} ']']);
        mu = eval(['[' result{parl} ']']);    
        
        
        
        
        % build the options structure
        % -----------------------------
        opt.refdata = EOG;
        

        opt.M = M;
        opt.mu = mu;
       

        % run the EOG correction
        % ----------------------

            tmp_in = reshape(EEG.icaact,EEG.nbchan,EEG.pnts*EEG.trials);

            [tmp,H,Hh] = lms_regica(tmp_in(EEG.contICs,:), opt);
            EEG.icaact(EEG.contICs,:)=tmp;

            EEG.Hh = Hh;
           
            EEG.data=EEG.icawinv*EEG.icaact;
            
     case 'crls' 
                     disp('pop_regica(): cRLS algorithm will be used for the rejection of ocular artifacts'); 
         
         if isempty(EEG.icaact)
            disp('pop_regica(): cannot clean an empty dataset'); 
            return;
         end;

        % try to guess which are the EOG channels
        % --------------------------------------------------
        if isempty(EOGindex),
            for i = 1:length(EEG.chanlocs),
                labels = EEG.chanlocs(i).labels;
                if ~isempty(strfind(lower(labels),'eog')),
                    EOGindex = [EOGindex i];
                end
            end
        end
        
        
        % reading params
        % -------------------
        EOGindex = eval(['[' result{pareog} ']']);
        M = eval(['[' result{parm} ']']);
        lambda = eval(['[' result{parff} ']']); 
        sigma = eval(['[' result{pars} ']']); 
        evchans = eval(['[' result{parchan} ']']);
        
        % build the options structure
        % -----------------------------
        opt.refdata = reshape(EOG,length(EOGindex),EEG.pnts*EEG.trials);
        opt.M = M;
        opt.lambda = lambda;
        opt.sigma = sigma;
        EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

        % run the EOG correction
        % ----------------------
        

        tmp_in = reshape(EEG.icaact,EEG.nbchan,EEG.pnts*EEG.trials);
        [tmp,H,Hh] = crls_regica(tmp_in(EEG.contICs,:), opt);
        EEG.data=EEG.icawinv*EEG.icaact;
        
        
     case 'srls' 
                     disp('pop_regica(): sRLS algorithm will be used for the rejection of ocular artifacts'); 
        if isempty(EEG.icaact)

            disp('pop_regica(): cannot clean an empty dataset'); 
            return;
        end;

        % try to guess which are the EOG channels
        % --------------------------------------------------
        if isempty(EOGindex),
            for i = 1:length(EEG.chanlocs),
                labels = EEG.chanlocs(i).labels;
                if ~isempty(strfind(lower(labels),'eog')),
                    EOGindex = [EOGindex i];
                end
            end
        end
        
        % reading params
        % -------------------
        EOGindex = eval(['[' result{pareog} ']']);
        M = eval(['[' result{parm} ']']);
        lambda = eval(['[' result{parff} ']']); 
        sigma = eval(['[' result{pars} ']']); 
        prec = eval(['[' result{12} ']']); 
        evchans = eval(['[' result{parchan} ']']);
        
        % build the options structure
        % -----------------------------
        opt.refdata = reshape(EOG,length(EOGindex),EEG.pnts*EEG.trials);
        opt.M = M;
        opt.lambda = lambda;
        opt.sigma = sigma;
        opt.prec = prec;
        EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

        % run the EOG correction
        % ----------------------

        tmp_in = reshape(EEG.icaact,EEG.nbchan,EEG.pnts*EEG.trials);
        [tmp,H,Hh] = scrls_regica(tmp_in(EEG.contICs,:), opt);
        EEG.icaact(EEG.contICs,:,:) = reshape(tmp,[length(EEG.contICs),EEG.pnts,EEG.trials]);    
        EEG.icaact(EEG.contICs,:)=tmp;
        EEG.data=EEG.icawinv*EEG.icaact;
        
     case 'h_inf_tv' 
         disp('pop_regica(): H_inf_tv algorithm will be used for the rejection of ocular artifacts'); 
         if isempty(EEG.icaact)
            disp('pop_regica(): cannot clean an empty dataset'); 
            return;
        end;

        % try to guess which are the EOG channels
        % --------------------------------------------------
        if isempty(EOGindex),
            for i = 1:length(EEG.chanlocs),
                labels = EEG.chanlocs(i).labels;
                if ~isempty(strfind(lower(labels),'eog')),
                    EOGindex = [EOGindex i];
                end
            end
        end
         
        % reading params
        % -------------------
        EOGindex = eval(['[' result{pareog} ']']);
        M = eval(['[' result{parm} ']']);
        eta = eval(['[' result{pareta} ']']);    
        rho = eval(['[' result{parrho} ']']);
        eps = eval(['[' result{pareps} ']']);
        evchans = eval(['[' result{parchan} ']']);
         
        % build the options structure
        % -----------------------------
        opt.refdata = reshape(EOG,length(EOGindex),EEG.pnts*EEG.trials);
        opt.M = M;
        opt.eta = eta;
        opt.rho = rho;
        opt.eps = eps;
        EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

        % run the EOG correction
        % ----------------------
        tmp_in = reshape(EEG.icaact,EEG.nbchan,EEG.pnts*EEG.trials);
        tmp = hinftv_regica(tmp_in(EEG.contICs,:), opt);
        EEG.data=EEG.icawinv*EEG.icaact;
         
         
     case 'h_inf_ew' 
            disp('pop_regica(): H_inf_ew algorithm will be used for the rejection of ocular artifacts'); 
            if isempty(EEG.icaact)
                disp('pop_regica(): cannot clean an empty dataset'); 
                return;
            end;
            EOGindex = [];
            % try to guess which are the EOG channels
            % --------------------------------------------------
            if isempty(EOGindex),
                for i = 1:length(EEG.chanlocs),
                    labels = EEG.chanlocs(i).labels;
                    if ~isempty(strfind(lower(labels),'eog')),
                        EOGindex = [EOGindex i];
                    end
                end
            end
            
            % reading params
            % -------------------
            EOGindex = eval(['[' result{pareog} ']']);
            M = eval(['[' result{parm} ']']);
            eta = eval(['[' result{pareta} ']']);    
            rho = eval(['[' result{parrho} ']']);
            eps = eval(['[' result{pareps} ']']);
            lambda = eval(['[' result{parff} ']']);
            evchans = eval(['[' result{parchan} ']']);
            
            % build the options structure
            % -----------------------------
            opt.refdata = reshape(EOG,length(EOGindex),EEG.pnts*EEG.trials);
            opt.M = M;
            opt.eta = eta;
            opt.rho = rho;
            opt.eps = eps;
            opt.lambda = lambda;
            EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

            % run the EOG correction
            % ----------------------
            
            tmp_in = reshape(EEG.icaact,EEG.nbchan,EEG.pnts*EEG.trials);
            tmp = hinfew_regica(tmp_in(EEG.contICs,:), opt);
            EEG.data=EEG.icawinv*EEG.icaact;
            

end;

if nargin < 2
    com = sprintf('%s = pop_regica(%s, %s);', inputname(1),inputname(1), ...
                  vararg2str({ 'regtype' g.regtype 'icatype' g.icatype 'dataset' g.dataset 'options' g.options }) );
end;

 



            if length(g.dataset) > 1
            for i = g.dataset
                ALLEEG(i).icaweights  = EEG.icaweights;
                ALLEEG(i).icasphere   = EEG.icasphere;
                ALLEEG(i).icawinv     = EEG.icawinv;
                ALLEEG(i).icaact      = EEG.icaact;
                ALLEEG(i).data      = EEG.data;
                ALLEEG(i).icachansind = g.chanind;
            end;            
                ALLEEG = eeg_checkset(ALLEEG);
            else
                EEG = eeg_checkset(EEG);
                ALLEEG = eeg_store(ALLEEG, EEG, g.dataset);
            end;
clear g;

return;
