function [Y,H,Hh] = hinfew_regica(X, opt)
% hinfew_regica() - Performs automatic EOG artifact correction using
% multiple adaptive regression. The adaptation is made using the H infinity
% exponentially weighted (EW) algorithm [1].
%
% Usage:
%   >>  [Y,H,Hh] = hinfew_regica(X, opt)
%
% Inputs:
%   X               - Input data matrix
%   opt             - Analysis options (see below)
%   opt.refdata     - Reference signal (s) (dref x N) (default: [])
%   opt.M           - Order of the adaptive filter (default: 5)
%   opt.eta         - Factor reflecting a priori knowledge of how close the
%                     estimated filter weights at t=0 are to their optimal
%                     value at that time instant (default: 5e-3)                       
%   opt.rho         - Factor reflecting a priori knowledge of how rapidly
%                     the filter coefficients vary with time
%                     (default: 1e-5)
%   opt.eps         - Positive constant described in [1] (default: 1.5)
%   opt.lambda      - Forgetting factor (default:0.9)
% 
% Outputs:
%   Y   - Output data matrix (dxN) (artifact corrected)
%   H   - Final filter weights (M*dref x d)
%   Hh  - filter weights evolution (M*dref x d x N)
%
% Notes:
%   1) This function implements the H infinity norm EW algorithm described
%      in [1].
%
% References:
% [1] S. Puthusserypady and T. Ratnarajah, IEEE Signal Processing Letters
% 12, 816-819
%
%
% Author: German Gomez-Herrero
%         http://www.cs.tut.fi/~gomezher/index.htm         
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% See also:
%   POP_hinfew_regica, EEGLAB
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
    help hinfew_regica;
    return;
end

if ~exist('opt','var'),
    opt = def_hinfew_regica;
else
    opt = def_hinfew_regica(opt);
end

if isempty(opt.refdata),
    error('(hinfew_regica) I need a reference signal!');
end

eta = opt.eta;
M = opt.M;
eps = opt.eps;
lambda = opt.lambda;
Xref = opt.refdata;
[deeg,Leeg] = size(X);
[dref,Lref] = size(Xref);
if Leeg~=Lref, 
    error('(hinfew_regica) Input data and reference signal must have the same length'); 
end

% initialization of the adaptation loop
% -----------------------------------------------
H = zeros(dref*M,deeg); 
Y = zeros(deeg,Leeg);
r = Xref(:,(M+1):-1:2);
r = reshape(r', M*dref,1);
invP_hat = (1/eta)*eye(M*dref)-eps^(-2)*r*r';

if nargout > 2, 
    Hh = zeros(dref*M,deeg,Leeg);
    Hh(:,:,1:M-1) = repmat(H,[1,1,M-1]); 
end

% adaptation loop
% -----------------------------------------------
if opt.verbose, fprintf('\n(hinfew_regica) '); end
for i = M:Leeg-1
    r = Xref(:,i:-1:(i-M+1));
    r = reshape(r', M*dref,1);
    r_next = Xref(:,i+1:-1:(i-M+2));
    r_next = reshape(r_next', M*dref,1);
    invP_hat = lambda*invP_hat+lambda*r*r'-eps^(-2)*r_next*r_next';
    P_hat = inv(invP_hat);
    g = (P_hat*r)./(1+r'*P_hat*r);    
    update = r'*H;
    if ~isempty(find(abs(update(:))>OVERFLOW,1)),
        error('(hinfew_regica) Algorithm became unstable');
    end
    Y(:,i) = (X(:,i)'-update)';
    e0 = X(:,i)'-r'*H; 
    H_prev = H;    
    H = H_prev+g*e0; 
    if nargout > 2, Hh(:,:,i) = H; end
    if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
end
if opt.verbose,fprintf('[OK]\n');end


function [opt] = def_hinfew_regica(opt)

if ~exist('opt','var') || isempty(opt) || ~isfield(opt,'refdata'),
    opt.refdata = [];
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt,'eta') || isempty(opt.eta),
    opt.eta = 5e-3;
end
if ~isfield(opt,'rho') || isempty(opt.rho),
    opt.rho = 1e-5;
end
if ~isfield(opt, 'eps') || isempty(opt.eps),
    opt.eps = 1.5;
end
if ~isfield(opt, 'M') || isempty(opt.M),
    opt.M = 3;
end
if ~isfield(opt, 'lambda') || isempty(opt.lambda),
    opt.lambda = 0.99;
end

