% eegplugin_regica() - EEGLAB plugin for Automatic EOG Artifact Rejection
% using REGICA v1.0
% 
% 
%
% Usage:
%   >> eegplugin_aar(fig, trystrs, catchstrs)
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks. 
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
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
% Citations: 
% 
% 1) M.A.Klados, et al.,REG-ICA:A hybrid methodology combining Blind
% Source Separation and regression techniques for the rejection of ocular artifacts, 
% Biomed. Signal Process. Control (2011), doi:10.1016/j.bspc.2011.02.001.
% 
% 2) Klados, M.A. Papadelis, C.L. Bamidis, P.D.,"REG-ICA: A new hybrid
% method for EOG Artifact Rejection",  9th International Conference on 
% Information Technology and Applications in Biomedicine, 2009. ITAB
% 2009.doi:10.1109/ITAB.2009.5394295
%
% Copyright (C) 2011 Manousos Klados, Aristotle University of Thessaloniki,
% Greece 
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

function vers = eegplugin_regica(fig, trystrs, catchstrs)

vers = 'REGICA v1.0';
if nargin < 3
    error('eegplugin_regica requires 3 arguments');
end;


% add plugin folder to path
% -----------------------
if exist('pop_regica.m','file')
    p = which('eegplugin_regica');
    p = p(1:findstr(p,'eegplugin_regica.m')-1);
    addpath(p);    
end;

% find tools menu
% ---------------------
menu = findobj(fig, 'tag', 'tools');




% menu callbacks
% --------------
regica_cback = [ trystrs.no_check '[EEG LASTCOM] = pop_regica(EEG);' catchstrs.new_and_hist ];



% create menus if necessary
% -------------------------
submenu = uimenu( menu, 'Label', 'Automatic EOG Artifact Rejection');
uimenu( submenu, 'Label', 'REGICA-Methodology', 'CallBack', regica_cback);



