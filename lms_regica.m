function [Y,H,Hout] = lms_regica(X, opt)
% lms_regica() - Performs automatic EOG artifact correction using
% multiple adaptive regression. The adaptation is made using the Least Mean
% Squares (LMS) algorithm [1].
%
% Usage:
%   >>  [Y,H,Hh] = lms_regica( X, opt)
%
% Inputs:
%   X               - Input data matrix (dxN)
%   opt             - Analysis options:
%   opt.refdata     - Reference signal (s) (dref x N) (default: [])
%   opt.M           - Order of the adaptive filter (default: 3)
%   opt.mu          - Learning rate (default: 1e-6)
% 
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%   H   - Final filter weights (M*dref x d)
%   Hh  - filter weights evolution (M*dref x d x N)
%
% References:
% [1] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
%
%
% Author: German Gomez-Herrero         
%         http://www.cs.tut.fi/~gomezher/index.htm
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% See also:
%   POP_lms_regica, EEGLAB
%

% Copyright (C) <2007>  <German Gomez-Herrero>
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

% OVERFLOW
OVERFLOW = 1E12;

if nargin < 1,
    help lms_regica;
    return;
end

if ~exist('opt','var'),
    opt = def_lms_regica;
else
    opt = def_lms_regica(opt);
end

if isempty(opt.refdata),
    error('(lms_regica) I need a reference signal!');
end

mu   = opt.mu;
M    = opt.M;
Xref = opt.refdata;
[deeg,Leeg] = size(X);
[dref,Lref] = size(Xref);
if Leeg~=Lref, 
    error('(lms_regica) Input data and reference signal must have the same length'); 
end

% Initialization of the adaptation loop
% ---------------------------------------------
H = zeros(dref*M,deeg); 
Y = zeros(deeg,Leeg);
if nargout > 2, 
    Hout = zeros(dref*M,deeg,Leeg);
    Hout(:,:,1:M-1) = repmat(H,[1,1,M-1]); 
end

% Adaptation loop
% ---------------------------------------------
if opt.verbose, fprintf('\n(lms_regica) '); end
for i = M:Leeg
    r = Xref(:,i:-1:(i-M+1));
    r = reshape(r', M*dref,1); 
    eupdate = r'*H;
    e0 = X(:,i)'-eupdate;    
    if ~isempty(find(abs(eupdate(:))>OVERFLOW, 1)),
        error('(lms_regica) Algorithm became unstable');
    end
    H = H+mu*r*e0;
    if nargout > 2, Hout(:,:,i) = H; end
    Y(:,i) = e0; 
    if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
end
if opt.verbose,fprintf('[OK]\n');end


% sub-function to initialize the default values
% ----------------------------------------------
function [opt] = def_lms_regica(opt)

if ~exist('opt','var') || isempty(opt) || ~isfield(opt,'refdata'),
    opt.refdata = [];
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt, 'M') || isempty(opt.M),
    opt.M = 3;
end
if ~isfield(opt, 'mu') || isempty(opt.mu),
    opt.mu = 1e-6;
end

